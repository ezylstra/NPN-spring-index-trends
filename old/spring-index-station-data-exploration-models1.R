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

# Packages below may not be needed
# library(mgcv)

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

# Fit trend models ------------------------------------------------------------#

# Just linear trends (at least for now)

# sdmTMB notes:
# Intercepts: can keep spatial = "on", a smooth spatial intercept, and 
  # (1 | station) for site-level departures from it, or drop either. Keeping a 
  # spatial intercept is usually wise. Otherwise, broad spatial differences in 
  # mean leaf timing can leak into the trend field.
# spatial_varying = ~ 0 + yrc: the field is a zero-mean deviation from the fixed 
  # yrz coefficient. The local trend at any location is the fixed slope plus the 
  # field value there. Wrote 0 + explicitly so there's no ambiguity about an 
  # intercept field being included twice.
# Site-level trend deviations: Can allow trends at each site to deviate
  # from the smooth trend surface with random effect terms: (1 + yrc | station) 
  # or (1 | station) + (0 + yrc | station). If those terms aren't included 
  # then trends at each site are determined by the overall mean (yrc) and 
  # the spatial smooth. If we include both, can assess how much trend variation
  # is spatially structure or not. 

# Start with a simpler model:
start <- Sys.time()
fit1 <- sdmTMB(
  leaf ~ yrc + (1 | station), 
  data = dat,
  mesh = mesh,
  spatial = "on",                  # spatial intercept field
  spatial_varying = ~ 0 + yrc,     # spatial field on the trend
  family = gaussian()
)
end <- Sys.time()
end - start # 2.8 min

# Check for problems
sanity(fit1)
# 2 warnings about high gradient values. Not sure this is a big deal, but will 
# run extra optimization in case.

fit1_update <- run_extra_optimization(fit1)
sanity(fit1_update)
# looks good

fit1_update
# Range = distance at which two locations are effectively independent (about 
  # 13% correlated). Here, 5500 km. Ideally, this value should be about 5-10x
  # the inner max.edge of mesh (500).
# Spatial SD (sigma_O, for “omega”): how much the spatial random field varies. 
# Spatially varying coefficient SD: Larger values mean more spatial variation in
  # trends among sites.
# See https://sdmtmb.github.io/sdmTMB/articles/basic-intro.html for helpful
  # info on interpreting results

# Can extract coefficients via tidy():
tidy(fit1_update, conf.int = TRUE)
tidy(fit1_update, effects = "ran_par", conf.int = TRUE)


# Try a model that includes random slope effects for each site (so slopes aren't
# solely determined by the spatial mesh). Keeping a complex model for intercepts
# (a smooth spatial intercept and (1 | station) for site-level departures from 
# it) for now.
start <- Sys.time()
fit2 <- sdmTMB(
  leaf ~ yrc + (1 | station) + (0 + yrc | station),  
  data = dat,
  mesh = mesh,
  spatial = "on",                  # spatial intercept field
  spatial_varying = ~ 0 + yrc,     # spatial field on the trend
  family = gaussian()
) 
end <- Sys.time()
end - start # 9.9 min

# Check for problems
sanity(fit2)
# 3 warnings about high gradient values. 
fit2_update <- run_extra_optimization(fit2)
sanity(fit2_update)
# just a warning about RE SDs

fit2_update
tidy(fit2_update, effects = "ran_par", conf.int = TRUE)
# The SD associated with the random trend effect is ~0. All other effects
# look almost identical to that from the model that didn't include that random
# effect. This suggests that most of the trend variation is spatially 
# structured. 

# Note also that the range is big (5500 km) relative to the distance between
# stations (esp in North America and Europe), so the random trend effect should
# have been free to capture local departures from the spatial field. 

# Also comparing sigma_Z to the station-level trend RE from the gamm4 model
# that didn't have a spatial field associated with trends. 
m0$mer
# Trend RE from m0 had a SE of 0.91, which is on the same order as sigma_Z in
# the sdmTMB model (2.61), which suggests that all the spatial variation in 
# trends can be captured by the gaussian smooth. 

# finally, model comparisons:
logLik(fit2_update); logLik(fit1_update)
# The log-likelihoods are the same. So no need to allow for station-level trend 
# variation beyond the spatial field.

# Will use fit1_update for inferences

# Map predictions -------------------------------------------------------------#

# Create a prediction grid:
# Unique site locations in projected km (no CRS needed; just planar coords)
sites <- st_as_sf(dat[!duplicated(dat$station), c("x", "y")], 
                  coords = c("x", "y"))

# Domain = everything within 200 km of a site
domain <- st_union(st_buffer(sites, 300))

# 50 km grid of cell centers, kept only inside the domain (can adjust this)
g <- st_make_grid(domain, cellsize = 50, what = "centers")
g <- g[lengths(st_intersects(g, domain)) > 0]

xy <- st_coordinates(g)
pred_grid <- data.frame(
  x = xy[, 1],
  y = xy[, 2],
  yrc = 0,
  # Pick a station ID (doesn't matter which one since we'll ignore station
  # level REs when making preds)
  station = factor(levels(dat$station)[1], levels = levels(dat$station))
)
nrow(pred_grid)

# Predict
p <- predict(fit1_update, newdata = pred_grid, re_form_iid = NA)
# Prediction is composed of multiple parts:
  # est_non_rf: the part from the fixed effects (here, with yrc = 0, its the Int)
  # omega_s: the spatial random field associated with the intercept
  # zeta_s_yrc: the spatially varying coefficient (yrc)
  # est: the full prediction

# To get local trends, we just need the spatially-varying coefficient predictions
p$local_trend <- fixef(fit1_update)["yrc"] + p$zeta_s_yrc   # change per decade

# Getting layer with country boundaries
laea_km <- "+proj=laea +lat_0=90 +lon_0=0 +units=km"
world <- ne_countries(scale = "medium", returnclass = "sf")

# Crop in lon/lat first; s2 off avoids spherical-geometry complaints on the crop
sf_use_s2(FALSE)
north <- st_crop(world, xmin = -180, xmax = 180, ymin = 20, ymax = 90)
sf_use_s2(TRUE)

north <- st_transform(north, laea_km)

# Calculate to make sure color ranmp is symmetric
lim <- max(abs(quantile(p$local_trend, c(0.02, 0.98))))

# Plot (Note that here we're "squishing" the extremes which sets the color 
# ramp limits at 2nd and 98th percentile and pins anything beyond that to 
# the end colors. Would need to mention this clampling in a figure caption)
ggplot() +
  geom_tile(data = p, aes(x, y, fill = local_trend), width = 50, height = 50) +
  geom_sf(data = north, fill = NA, color = "grey30", linewidth = 0.25) +
  scale_fill_distiller(palette = "RdBu", limits = c(-lim, lim),
                       oob = scales::squish, direction = 1, 
                       name = "Trend\n(per decade)") +
  # scale_fill_gradient2(name = "Trend\n(per decade)") +
  coord_sf(crs = laea_km,
           xlim = range(p$x), ylim = range(p$y), expand = FALSE) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

### TO DO: Clip prediction grid to land areas...
# Reducing the domain buffer doesn't work well since that creates a bunch of 
# isolated prediction islands in Russia.

# Assess local trend signficance ----------------------------------------------#

# From Claude:
# sdmTMB can simulate draws from the joint uncertainty of the fixed trend and 
# the spatial trend field. From those draws you get a credible interval for the 
# local trend at every grid cell and can flag cells where it excludes zero.

# The local trend is the fixed slope plus the value of the spatial field at that 
# location. A clean way to get draws of exactly that is to predict each cell at 
# yrc = 0 and yrc = 1 in the same call, so both rows use the same parameter 
# draws, and then subtract:

n <- nrow(pred_grid)
nd <- rbind(transform(pred_grid, yrc = 0),
            transform(pred_grid, yrc = 1))

set.seed(1)
sims <- predict(fit1_update, newdata = nd, nsim = 500, re_form_iid = NA)

tr <- sims[(n + 1):(2 * n), ] - sims[1:n, ]   # trend draws: cells x 500

pred_grid$trend <- rowMeans(tr)
pred_grid$lwr   <- apply(tr, 1, quantile, 0.025)
pred_grid$upr   <- apply(tr, 1, quantile, 0.975)
pred_grid$p_pos <- rowMeans(tr > 0)
pred_grid$sig   <- pred_grid$lwr > 0 | pred_grid$upr < 0

# Check that pred_grid$trend is very simliar to p$local_trend
summary(pred_grid$trend - p$local_trend)
# Confirming that differences are just simulation noise...
point <- fixef(fit1)["yrc"] + p$zeta_s_yrc
mc_se <- apply(tr, 1, sd) / sqrt(ncol(tr))
z <- (pred_grid$trend - point) / mc_se
summary(z)
mean(abs(z) > 2)

#### TODO: These are pointwise intervals, so we have the same multiple 
# comparison issues that have come up before. Figure out if/how we want to deal
# with that. 

# In map, should see that uncertainty grows as you move away from stations.

# Map just showing where trends are significantly different than zero
ggplot() +
   #geom_tile(data = p, aes(x, y, fill = local_trend), width = 50, height = 50) +
  geom_tile(data = pred_grid, aes(x, y, fill = trend, alpha = sig),
            width = 50, height = 50) +
  geom_sf(data = north, fill = NA, color = "grey30", linewidth = 0.25) +
  geom_contour(data = pred_grid, aes(x, y, z = as.numeric(sig)),
               breaks = 0.5, color = "black", linewidth = 0.4) +
  scale_alpha_manual(values = c("FALSE" = 0, "TRUE" = 1), guide = "none") +
  scale_fill_distiller(palette = "RdBu", limits = c(-lim, lim),
                       oob = scales::squish, direction = 1, 
                       name = "Trend\n(per decade)") +
  coord_sf(crs = laea_km,
           xlim = range(p$x), ylim = range(p$y), expand = FALSE) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

# TODO: Try different meshes to see if that effects estimates much
# TODO: Figure out how to calculate hemispheric trends

