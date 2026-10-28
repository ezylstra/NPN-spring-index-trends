# Exploration of SI-x station data

library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(sf)
library(gamm4)
library(sdmTMB)
library(fmesher)
library(ggplot2)
library(rnaturalearth)
library(RColorBrewer)

# Using Claude to help sketch out workflow (including coordinate projection,
# model specification, extraction of site-level trends, and estimation of
# hemisphere-level trend)

# Load, format data -----------------------------------------------------------#
 
dfload <- read.csv("SI-x_NH_1950_2024_25.csv", na.strings = c("", 999, "999"))
# GHCN = 3-digit country code
# Stat = (up to 7 digit) station number within each country  
# Need combined 10-digit number to uniquely ID a station
# 4 indices: leaf, bloom, last freeze, damage (= diff between leaf, last freeze)

# Make all column names lowercase
colnames(dfload) <- str_to_lower(colnames(dfload))

# Create station ID, then remove unnecessary columns and format
df <- dfload %>%
  mutate(station = paste0(ghcn, "-", sprintf("%07d", stat))) %>%
  rename(leaf = six_fleaf,
         leaf_norm = six_fleaf_norm,
         bloom = six_fbloom,
         bloom_norm = six_fbloom_norm,
         frz = lstfrz28,
         frz_norm = lstfrz28_norm,
         lat = latitude,
         lon = longitude,
         elev = elevation) %>%
  select(ghcn, station, lat, lon, elev, year, leaf, bloom, frz, damage, 
         leaf_norm, bloom_norm, frz_norm, damage_norm)

# Summarize data by station ---------------------------------------------------#

stations <- df %>%
  group_by(ghcn, station, lat, lon, elev) %>%
  summarize(nyrs = n_distinct(year),
            n_leaf = sum(!is.na(leaf)),
            n_bloom = sum(!is.na(bloom)),
            n_frz = sum(!is.na(frz)),
            n_damage = sum(!is.na(damage)),
            .groups = "drop") %>%
  mutate(prop_frz = n_frz/nyrs) %>%
  data.frame()

dim(stations) # 4080 stations
summary(stations) # max of 75 years
all.equal(stations$nyrs, stations$n_leaf)
all.equal(stations$nyrs, stations$n_bloom)
# There's a leaf & bloom value for every year & station
# There's a freeze damage value for most years & stations (a few stations have
# no freeze values)

# Project site coordinates ----------------------------------------------------#

# Want site coords to be projected (not raw lon/lat) because sites span a wide
# area and gaussian process covariance assumes a roughly Euclidean distance 
# metric. Will use pole-centered Lambert azimuthal equal-area (polar LAEA)
# (EPSG:6931, "EASE-Grid 2.0 North")

# We'll also want units in km, not m, given the calculations needed for the
# gaussian surface

sites_sf <- st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
# sites_ease <- st_transform(sites_sf, crs = 6931)
sites_ease <- st_transform(sites_sf, "+proj=laea +lat_0=90 +lon_0=0 +units=km")
# plot(st_geometry(sites_ease))

df <- df %>%
  mutate(x = st_coordinates(sites_ease)[, 1],
         y = st_coordinates(sites_ease)[, 2])

# Create subset of data for initial model exploration -------------------------#
# For now, using raw index values, not deviations from 30-year norms

# Data filters
index <- "leaf"
years <- 1991:2024

# Create new years variables; make sure station is a factor
dat <- df %>%
  filter(year %in% years) %>%
  mutate(yr0 = year - min(year)) %>%
  mutate(yrc = (year - mean(year))/10) %>%  # centered, scaled to decade
  mutate(station = factor(station)) %>%
  select(ghcn, station, lat, lon, x, y, year, yr0, yrc, contains(index)) %>%
  arrange(station, year)

# Check for temporal autocorrelation ------------------------------------------#

# Check whether we need to include an AR1 component in models (to account for
# temporal autocorrelation within sites). To do this, fit a model without an AR1 
# term, look at the lag-1 autocorrelation of the residuals within sites (and 
# then could use this estimate as rho)

# Using gamm4 and not bam for this. gamm4 fits the random effects through lme4, 
# which uses sparse matrices, so having lots of sites is cheap. The GP smooth is 
# associated with the intercepts and handled like mgcv. 

start <- Sys.time()
m0 <- gamm4(leaf ~ yrc + s(x, y, bs = "gp", m = 2, k = 100),
            random = ~ (1 | station) + (0 + yrc | station),
            data = dat)
end <- Sys.time()
end - start # Run time = 14.5 min

# check lag-1 correlation
dat$r <- residuals(m0$mer)
dat <- dat[order(dat$station, dat$yrc), ]
lag1 <- with(dat, ave(r, station, FUN = function(z) c(NA, head(z, -1))))
cor(dat$r, lag1, use = "complete.obs")
# Value = -0.02 (and values less than ~0.1 suggest no meaningful correlation), 
# so really no evidence of temporally correlated errors
# Leaving out AR1 terms in subsequent models.

# Model trends while accounting for spatial autocorrelation -------------------#

# We have a few options for including spatial smooths: 
# mgcv::bam, gamm4::gamm4, sdmTMB::sdmTMB, maybe others?

# Using Template Model Builder (TMB) via sdmTMB might be the fastest option for
# models with lots of random effects (and Gaussian random fields), like we do. 
# sdmTMB can include spatial and spatiotemporal latent Gaussian Markov random
# fields using the stochastic partial differential equation (SPDE) approach.
# SPDE matrices are constructed with fmesher, and estimation is conducted via 
# maximum marginal likelihood with TMB. See: 
# https://www.jstatsoft.org/article/view/v115i02
# https://sdmtmb.github.io/sdmTMB/

# Create mesh for Gaussian random field ---------------------------------------#

# For initial exploration, we want ~800-2000 vertices. The interior of the mesh 
# should contain most of the vertices/triangles, with an outer buffer area that 
# has larger triangles (created to avoid edge effects)

inla_mesh <- fm_mesh_2d_inla(
  loc      = as.matrix(dat[, c("x", "y")]),
  max.edge = c(250, 1000),   # max triangle edge (km): inside data area, outer buffer
  offset   = c(250, 1000),   # inner and outer buffer widths (km)
  cutoff   = 75              # merge vertices closer than this (km)
)
mesh <- make_mesh(dat, c("x", "y"), mesh = inla_mesh)
plot(mesh)
mesh$mesh$n   # 5618 vertices 
# This is probably overkill. Would like to reduce this so fits don't take too long

# Look at a range of options.
loc <- as.matrix(dat[, c("x", "y")])
settings <- expand.grid(inner = c(300, 400, 500, 700), cut_frac = c(0.3, 0.5))
settings$n <- apply(settings, 1, function(s) {
  m <- fm_mesh_2d_inla(
    loc      = loc,
    max.edge = c(s["inner"], s["inner"] * 4),
    offset   = c(s["inner"], s["inner"] * 3),
    cutoff   = s["inner"] * s["cut_frac"]
  )
  m$n
})
settings
# Selecting option 7:
inner <- settings$inner[7]      
cut_frac <- settings$cut_frac[7]

inla_mesh <- fm_mesh_2d_inla(
  loc      = loc,
  max.edge = c(inner, inner * 4), # max.edge = c(500, 2000)
  offset   = c(inner, inner * 3), # offset = c(500, 1500)
  cutoff   = inner * cut_frac     # cutoff = 250
)
mesh <- make_mesh(dat, c("x", "y"), mesh = inla_mesh)
mesh$mesh$n # 1495 vertices
plot(mesh)
# Looks reasonable

