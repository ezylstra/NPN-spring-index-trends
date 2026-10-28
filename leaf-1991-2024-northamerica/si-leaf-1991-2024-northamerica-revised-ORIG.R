# Exploration of SI-x station data

library(dplyr)
library(stringr)
library(tidyr)
library(sf)
library(sdmTMB)
library(fmesher)
library(ggplot2)
library(rnaturalearth)

# Load, format data -----------------------------------------------------------#

dfload <- read.csv("SI-x_NH_1950_2024_25.csv", na.strings = c("", 999, "999"))

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

# Subset data to North America ------------------------------------------------#

# ghcns <- df %>%
#   group_by(ghcn) %>%
#   summarize(nstations = n_distinct(station), 
#             minlat = min(lat),
#             maxlat = max(lat),
#             minlon = min(lon),
#             maxlon = max(lon)) %>%
#   data.frame()
# ghcns
# Looks like US and Canada stations all have ghcn = 403 or 425

dfna <- df %>%
  filter(ghcn %in% c(403, 425))

# Project site coordinates ----------------------------------------------------#

# Want site coords to be projected (not raw lon/lat) because sites span a wide
# area and gaussian process covariance assumes a roughly Euclidean distance 
# metric. Will use a Lambert azimuthal equal-area projection centred on the
# North America

# We'll also want units in km, not m, given the calculations needed for the
# gaussian surface

laea_na <- "+proj=laea +lat_0=50 +lon_0=-100 +datum=WGS84 +units=km"

sites_sf <- st_as_sf(dfna, coords = c("lon", "lat"), crs = 4326)
sites_ease <- st_transform(sites_sf, laea_na)
# plot(st_geometry(sites_ease))

dfna <- dfna %>%
  mutate(x = st_coordinates(sites_ease)[, 1],
         y = st_coordinates(sites_ease)[, 2])

# Create subset of data for initial model exploration -------------------------#
# For now, using raw index values, not deviations from 30-year norms

# Data filters
index <- "leaf"
years <- 1991:2024

# Create new years variables; make sure station is a factor
dat <- dfna %>%
  filter(year %in% years) %>%
  mutate(yrc = (year - mean(year))/10) %>%  # centered, scaled to decade
  mutate(station = factor(station)) %>%
  select(ghcn, station, lat, lon, x, y, year, yrc, contains(index)) %>%
  arrange(station, year)

yr_mean <- mean(dat$year)    

# Create mesh for Gaussian random field ---------------------------------------#

# For initial exploration, we want ~800-2000 vertices. The interior of the mesh 
# should contain most of the vertices/triangles, with an outer buffer area that 
# has larger triangles (created to avoid edge effects)

# Look at a range of options.
loc <- as.matrix(dat[, c("x", "y")])
# settings <- expand.grid(inner = c(300, 400, 500, 700), cut_frac = c(0.3, 0.5))
# settings$n <- apply(settings, 1, function(s) {
#   m <- fm_mesh_2d_inla(
#     loc      = loc,
#     max.edge = c(s["inner"], s["inner"] * 4),
#     offset   = c(s["inner"], s["inner"] * 3),
#     cutoff   = s["inner"] * s["cut_frac"]
#   )
#   m$n
# })
# settings

inner <- 300   
cut_frac <- 0.5

inla_mesh <- fm_mesh_2d_inla(
  loc      = loc,
  max.edge = c(inner, inner * 4),
  offset   = c(inner, inner * 3),
  cutoff   = inner * cut_frac
)
mesh <- make_mesh(dat, c("x", "y"), mesh = inla_mesh)
# mesh$mesh$n # 1195 vertices
# plot(mesh)
# Looks reasonable

# Model trend model with yearly spatial anomaly fields ------------------------#

# using yrc (centered year, scaled to decade)
# time = "year" + spatiotemporal = "iid" adds a zero-mean GRF for each year.

fit1 <- sdmTMB(
  leaf ~ yrc + (1 | station),
  data = dat,
  mesh = mesh,
  time = "year",
  spatial = "on",
  spatiotemporal = "iid",
  spatial_varying = ~ 0 + yrc,
  family = gaussian()
)
sanity(fit1)
# tidy(fit1, effects = "fixed", conf.int = TRUE) 
# tidy(fit1, effects = "ran_par", conf.int = TRUE) 

# Create a prediction grid ----------------------------------------------------#

# Unique site locations in projected km (no CRS needed; just planar coords)
sites <- st_as_sf(dat[!duplicated(dat$station), c("x", "y")], 
                  coords = c("x", "y"))

# Domain = everything within 300 km of a site
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

# Getting layer with country boundaries
world <- ne_countries(scale = "medium", returnclass = "sf")

# Crop in lon/lat first; s2 off avoids spherical-geometry complaints on the crop
sf_use_s2(FALSE)
north <- st_crop(world, xmin = -180, xmax = -50, ymin = 20, ymax = 75)
sf_use_s2(TRUE)
north <- st_transform(north, laea_na)

# Land-only prediction grid
# Keep the 300-km buffer around stations and intersect it with land. This
# restricts the grid to land that has nearby stations
sf_use_s2(FALSE)
land <- st_union(st_make_valid(north))
sf_use_s2(TRUE)
pg_sf <- st_as_sf(pred_grid, coords = c("x", "y"), crs = laea_na, remove = FALSE)
pred_grid$land <- lengths(st_intersects(pg_sf, land)) > 0

lg <- pred_grid[pred_grid$land, c("x", "y", "station")]
n_land <- nrow(lg)
n_land # 5839 cells

# check: 
# ggplot() +
#   geom_sf(data = north, fill = NA, color = "grey30", linewidth = 0.25) +
#   geom_point(data = lg, aes(x = x, y = y))

# Settings ---------------------------------------------------------------------#

nsim <- 1000  # consider 2000 for lfsr tails
mem  <- n_land * length(years) * nsim * 8 / 1e9
coarsen <- mem > 6  # set relative to your HPC memory request

if (coarsen) {
  keep <- (round((lg$x - min(lg$x)) / 50) %% 2 == 0) &
    (round((lg$y - min(lg$y)) / 50) %% 2 == 0)
  lg_c <- lg[keep, c("x", "y", "station")]
} else {
  lg_c <- lg[, c("x", "y", "station")]
}
n_c <- nrow(lg_c)

# Build one newdata: [trend rows at yrc = 0] [trend rows at yrc = 1] [cell x year]
nd <- rbind(transform(lg, yrc = 0, year = years[1]),
            transform(lg, yrc = 1, year = years[1]))
grid_yrs_c <- as.data.frame(tidyr::expand_grid(lg_c, year = years)) %>%
  mutate(yrc = (year - yr_mean) / 10)
newdat <- rbind(nd, grid_yrs_c[, names(nd)])

i0 <- seq_len(n_land)
i1 <- n_land + seq_len(n_land)
ig <- 2 * n_land + seq_len(nrow(grid_yrs_c))

# One joint simulation call ---------------------------------------------------#
set.seed(1234)
sims <- predict(fit1, newdata = newdat, nsim = nsim, re_form_iid = NA)

tr <- sims[i1, ] - sims[i0, ]  # n_land x nsim
ann_draws <- rowsum(sims[ig, , drop = FALSE], grid_yrs_c$year) / n_c  # years x nsim
rm(sims); gc()

# Point-estimate local trend (needs the model, so do it here) -----------------#
p1 <- predict(fit1, newdata = transform(lg, year = years[1], yrc = 0),
              re_form_iid = NA)
lg$trend_pt <- fixef(fit1)["yrc"] + p1$zeta_s_yrc

# Hemisphere-wide summaries ---------------------------------------------------#
underlying_trend <- colMeans(tr)  # per draw, days/decade

ann_sum <- data.frame(
  year = as.numeric(rownames(ann_draws)),
  est  = rowMeans(ann_draws),
  lwr  = apply(ann_draws, 1, quantile, 0.025),
  upr  = apply(ann_draws, 1, quantile, 0.975)
)
xc <- (ann_sum$year - yr_mean) / 10
w  <- (xc - mean(xc)) / sum((xc - mean(xc))^2)
realized_trend <- as.vector(crossprod(w, ann_draws))  # per draw, days/decade

# Station info for maps and methods --------------------------------------------#
stations <- dat %>% distinct(station, lat, lon, x, y)
n_by_year <- dat %>% filter(!is.na(leaf)) %>% count(year, name = "n_stations")

# Save -------------------------------------------------------------------------#
saveRDS(fit1, "output/si-leaf-1991-2024-northamerica-fit.rds")

saveRDS(list(lg = lg, # cell coords + trend_pt; rows match tr
             tr = tr,
             ann_draws = ann_draws,
             ann_sum = ann_sum,
             underlying_trend = underlying_trend,
             realized_trend = realized_trend,
             stations = stations,
             n_by_year = n_by_year,
             fit_sanity = fit_sanity,
             fit_fixed = fit_fixed,
             fit_ranpar = fit_ranpar,
             years = years, yr_mean = yr_mean, laea_na = laea_na,
             n_land = n_land, n_cells_annual = n_c,
             cellsize_km_annual = ifelse(coarsen, 100, 50),
             nsim = nsim, seed = 1234,
             session = sessionInfo()),
        "output/si-leaf-1991-2024-northamerica-results.rds")
