# Exploration of SI-x station data

library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(sf)
library(sdmTMB)
library(fmesher)
library(ggplot2)
library(rnaturalearth)
library(RColorBrewer)
library(lme4)
library(scales)
library(here)

# Location of saved files
outputloc <- "leaf-1991-2024-northamerica/job24387720-output/"

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

# Subset data to North America ------------------------------------------------#

ghcns <- df %>%
  group_by(ghcn) %>%
  summarize(nstations = n_distinct(station), 
            minlat = min(lat),
            maxlat = max(lat),
            minlon = min(lon),
            maxlon = max(lon)) %>%
  data.frame()
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
  mutate(yr0 = year - min(year)) %>%
  mutate(yrc = (year - mean(year))/10) %>%  # centered, scaled to decade
  mutate(station = factor(station)) %>%
  select(ghcn, station, lat, lon, x, y, year, yr0, yrc, contains(index)) %>%
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

# fit1 <- sdmTMB(
#   leaf ~ yrc + (1 | station),
#   data = dat,
#   mesh = mesh,
#   time = "year",
#   spatial = "on",
#   spatiotemporal = "iid",
#   spatial_varying = ~ 0 + yrc,
#   family = gaussian()
# )

# Save to file
# saveRDS(fit1, "output/leaf-na-1991-2024-model.rds")

fit1 <- readRDS(paste0(outputloc, "leaf-na-1991-2024-model.rds"))

sanity(fit1)
tidy(fit1, effects = "fixed", conf.int = TRUE) 
tidy(fit1, effects = "ran_par", conf.int = TRUE) 

# Create a prediction grid ----------------------------------------------------#

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

# 3. Local trend draws (land cells) -------------------------------------------#
# Predict each cell at yrc = 0 and yrc = 1 in the same year, using the same
# draws, and take the difference. The yearly field and the intercept field
# cancel, so each draw is exactly b + zeta(s).
# Use more draws than 500. The lfsr/FDR step depends on the tails.

# nsim <- 1000 # 2000 is probably better ########
# nd <- rbind(transform(lg, yrc = 0, year = years[1]),
#             transform(lg, yrc = 1, year = years[1]))
# set.seed(1)
# sims <- predict(fit1, newdata = nd, nsim = nsim, re_form_iid = NA)
# tr <- sims[(n_land + 1):(2 * n_land), ] - sims[1:n_land, ]   # n_land x nsim
# rm(sims)

# Save tr
# saveRDS(tr, "output/si-leaf-1991-2024-northamerica-tr.rds")

tr <- readRDS(paste0(outputloc, "si-leaf-1991-2024-northamerica-tr.rds"))

lg$trend <- rowMeans(tr)
lg$sd    <- apply(tr, 1, sd)
lg$lwr   <- apply(tr, 1, quantile, 0.025)
lg$upr   <- apply(tr, 1, quantile, 0.975)
lg$p_pos <- rowMeans(tr > 0)

# Create a map with local trends ----------------------------------------------#
# Not worrying about estimating signficance for the time being

# One row per land cell. year must be one of the fitted years (fit1 has
# time = "year"), but which year you pick doesn't affect zeta_s_yrc.
lg_pred <- transform(lg[, c("x", "y", "station")], year = years[1], yrc = 0)

p1 <- predict(fit1, newdata = lg_pred, re_form_iid = NA)
p1$local_trend <- fixef(fit1)["yrc"] + p1$zeta_s_yrc   # change per decade

lim1 <- max(abs(quantile(p1$local_trend, c(0.02, 0.98))))

trendmap <- ggplot() +
  geom_tile(data = p1, aes(x, y, fill = local_trend), width = 50, height = 50) +
  geom_sf(data = north, fill = NA, color = "grey30", linewidth = 0.25) +
  scale_fill_distiller(palette = "RdBu", limits = c(-lim1, lim1),
                       oob = scales::squish, direction = 1,
                       name = "Trend\n(per decade)") +
  coord_sf(crs = laea_na, xlim = range(p1$x), ylim = range(p1$y),
           expand = FALSE) +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank())

# Same idea, but using mean values from tr
trendmap2 <- ggplot() +
  geom_tile(data = lg, aes(x, y, fill = trend), width = 50, height = 50) +
  geom_sf(data = north, fill = NA, color = "grey30", linewidth = 0.25) +
  scale_fill_distiller(palette = "RdBu", limits = c(-lim1, lim1),
                       oob = scales::squish, direction = 1,
                       name = "Trend\n(per decade)") +
  coord_sf(crs = laea_na, xlim = range(p1$x), ylim = range(p1$y),
           expand = FALSE) +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank())

trendmap
trendmap2
# Very similar but not identical. 
# They estimate the same quantity, b + ζ(s), in two ways:
  # trendmap (p1) uses the model’s point estimate: the fixed effect plus the 
    # conditional mode of ζ. It is exact and has no simulation noise.
  # trendmap2 (lg$trend) is the mean of 1,000 draws from the model’s approximate 
    # joint distribution. That distribution is centered on the same estimates, 
    # so the mean converges to trendmap as nsim grows. The remaining differences 
    # are Monte Carlo noise, roughly sd / sqrt(1000) per cell.
# Use trendmap for point estimates, and use the draws for uncertainty and 
# significance.

# 3a. Bayesian FDR on the sign (local false sign rate) ------------------------#
# Here’s the logic, step by step:
  # p_pos is the posterior probability that a cell’s trend is above 0. It’s the 
    # fraction of the 1,000 draws that are positive.
  # lfsr (local false sign rate) is the posterior probability that the sign of 
    # the trend is wrong and is min(p_pos, 1 − p_pos). A cell with p_pos = 0.97 
    # has lfsr = 0.03.
  # Rank the cells by lfsr from most to least certain, and keep adding cells 
    # while the average lfsr of the flagged set stays ≤ 0.05.
# All cells have equal area (equal-area projection), so this is an area-based
# rate.

q <- 0.05
lg$lfsr <- pmin(lg$p_pos, 1 - lg$p_pos)
o  <- order(lg$lfsr)
cm <- cumsum(lg$lfsr[o]) / seq_along(o)
k  <- max(c(0, which(cm <= q)))
lg$sig_fdr <- FALSE
lg$sig_fdr[o[seq_len(k)]] <- TRUE
mean(lg$sig_fdr)                 # proportion of land area flagged
mean(lg$lfsr[lg$sig_fdr])        # realized average lfsr of the flagged set

# Map: FDR-flagged cells
lim <- max(abs(quantile(p1$local_trend, c(0.02, 0.98))))
ggplot() +
  geom_tile(data = lg, aes(x, y, fill = trend, alpha = sig_fdr),
            width = 50, height = 50) +
  geom_sf(data = north, fill = NA, color = "grey30", linewidth = 0.25) +
  scale_alpha_manual(values = c("FALSE" = 0.15, "TRUE" = 1), guide = "none") +
  scale_fill_distiller(palette = "RdBu", limits = c(-lim, lim),
                       oob = scales::squish, direction = 1,
                       name = "Trend\n(per decade)") +
  coord_sf(crs = laea_na, xlim = range(lg$x), ylim = range(lg$y),
           expand = FALSE) +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank())

# 3b. Simultaneous 95% band (stricter; all cells jointly) --------------------
# Max-standardized-deviation method: the band covers the whole trend surface
# with ~95% probability.
zmax <- apply(abs((tr - lg$trend) / lg$sd), 2, max)
c95  <- quantile(zmax, 0.95)
c95 # compare with 1.96 for pointwise intervals
lg$sig_simul <- (lg$trend - c95 * lg$sd > 0) | (lg$trend + c95 * lg$sd < 0)
mean(lg$sig_simul)

# Comparing 3a and 3b:
# They control different error rates, much like Benjamini–Hochberg versus 
# Bonferroni:
  # 3b (simultaneous band) controls the chance of making any wrong call anywhere 
    # on the map at about 5%. It’s very strict and flags fewer cells.
  # 3a (FDR) allows some wrong calls but limits them to about 5% of the flagged 
    # area. It’s less strict and flags more cells.
# A reasonable way to report them: use the FDR map to show where trends appear
  # and if you want, , and point to the simultaneous-band cells as the areas 
  # you’re most confident about.

# 4. Hemisphere-wide trend ----------------------------------------------------#
# Area-weighted mean of the local trend over land cells, computed in each draw.
# This is not the same as fixef(yrc). The SVC field is centered over the whole
# mesh, not over land, and the station network is dominated by the US.
hemi_trend <- colMeans(tr)
c(est = mean(hemi_trend), quantile(hemi_trend, c(0.025, 0.975)))
# Units: index units (days) per decade

# Averaged over land, it's the model’s underlying trend: -0.154 (-0.99, 0.65).
# This does NOT take the yearly anomaly fields into account.

# 5. Hemisphere-wide annual values -------------------------------------------#
# Predict every land cell in every year, then take the area-weighted mean.
# area = 1/n_land turns get_index()'s area-weighted sum into a mean. The SEs
# come from TMB's generalized delta method. get_index() builds its CIs from
# log_est and se on the log scale. That is fine for positive day-of-year
# values. For symmetric intervals, use est +/- 1.96 * se * est (delta method).
# If get_index() is slow, coarsen the grid (e.g., 100 km). Equal-area cells mean
# coarsening does not bias the mean.

grid_yrs <- tidyr::expand_grid(lg[, c("x", "y", "station")], year = years) %>%
  mutate(yrc = (year - yr_mean) / 10)
# hemi <- get_index(obj = fit1, 
#                   newdata = grid_yrs, 
#                   predict_args = list(re_form_iid = NA), 
#                   area = rep(1 / n_land, nrow(grid_yrs)),
#                   bias_correct = FALSE)

# Save hemi
# saveRDS(hemi, "output/si-leaf-1991-2024-northamerica-hemi.rds")
hemi <- readRDS(paste0(outputloc, "si-leaf-1991-2024-northamerica-hemi.rds"))

# Create plot
ggplot(hemi, aes(x = year)) +
  geom_segment(aes(xend = year, y = lwr, yend = upr), color = "gray") +
  geom_point(aes(y = est)) +
  geom_smooth(aes(y = est), method = "lm")
summary(lm(est ~ I((year - yr_mean)/10), data = hemi)) 
confint(lm(est ~ I((year - yr_mean)/10), data = hemi))
# Trend = 0.0789 (-1.05, 1.21)

# Why is this different than hemi_trend from section 4?

# First, this DOES take the annual anomalies into account (re_form_iid = NA 
  # ignores the (1 |station) effect but not the annual GRFs. Spatial and 
  # spatiotemporal fields are controlled by the re_form argument)

# Predict on grid, taking everything but random intercepts into account:
pg <- predict(fit1, newdata = grid_yrs, re_form_iid = NA)
# Get average doy for each year and the average yearly anomaly
ann <- pg %>%
  group_by(year) %>%
  summarize(total = mean(est), eps = mean(epsilon_st), .groups = "drop") %>%
  mutate(yrc = (year - yr_mean) / 10) %>%
  data.frame()

# (a) Do these annual means match get_index()?
summary(ann$total - hemi$est[match(ann$year, hemi$year)])   # should be ~0

# (b) The three pieces
slope_total <- coef(lm(total ~ yrc, data = ann))[2]
slope_eps   <- coef(lm(eps   ~ yrc, data = ann))[2]
trend_point <- fixef(fit1)["yrc"] + mean(pg$zeta_s_yrc[pg$year == years[1]])

c(slope_total = slope_total, trend_point = trend_point, slope_eps = slope_eps,
  sum = trend_point + slope_eps)   # slope_total should equal sum
mean(hemi_trend)                   # should be close to trend_point


# Notes: 
# pg (from predict()) has one row per cell per year. Each row holds the point 
# prediction (est) plus its components (omega_s, zeta_s_yrc, epsilon_st, etc). 
# It has no uncertainty: no SEs and no intervals. Averaging est over cells 
# within a year gives you the annual means, but nothing about how precise they 
# are.

# hemi (from get_index()) has one row per year. Internally, get_index() runs the 
# same prediction on grid_yrs. It then does the area-weighted sum inside TMB, 
# the model’s estimation engine, which lets it compute an SE for each annual 
# value with the generalized delta method. That SE accounts for joint 
# uncertainty in all the fixed effects and random fields. The trade-off is that 
# the cell-level detail and the components are gone.

# So the relationship is:
# hemi$est = the land average of pg$est, plus a properly propagated SE for each year.
