# SI-x leaf index, North America 1991-2024: PRODUCTION RUN
#
# Assumes si-leaf-1991-2024-northamerica-validate.R has been run and all
# checks passed. Fits the model (or reuses the saved fit), then builds all
# products from joint parameter draws projected with a sparse A matrix:
#   tr               n_land x nsim   cell trend draws (days/decade)
#   lg$trend_pt      point-estimate local trend
#   ann_draws        years x nsim    land-averaged annual index draws
#   ann_sum          annual means and 95% intervals
#   underlying_trend, realized_trend (per draw, same draws)
# Peak memory during simulation is set by `batch`, not by cells x years x nsim.

library(dplyr)
library(stringr)
library(tidyr)
library(sf)
library(sdmTMB)
library(fmesher)
library(Matrix)
library(rnaturalearth)

dir.create("output", showWarnings = FALSE)

# Settings --------------------------------------------------------------------#

use_saved_fit <- TRUE   # reuse the fit from the validation run if it exists
fit_file <- "output/si-leaf-1991-2024-northamerica-fit.rds"
out_file <- "output/si-leaf-1991-2024-northamerica-results.rds"

nsim  <- 1000   # consider 2000 for lfsr tails
batch <- 100    # draws held in memory at once; lower if memory is tight
seed  <- 1234

# Load, format data -----------------------------------------------------------#

dfload <- read.csv("SI-x_NH_1950_2024_25.csv", na.strings = c("", 999, "999"))
colnames(dfload) <- str_to_lower(colnames(dfload))

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

# US and Canada stations have ghcn = 403 or 425
dfna <- df %>% filter(ghcn %in% c(403, 425))

# Project site coordinates (km) -----------------------------------------------#

laea_na <- "+proj=laea +lat_0=50 +lon_0=-100 +datum=WGS84 +units=km"

sites_sf   <- st_as_sf(dfna, coords = c("lon", "lat"), crs = 4326)
sites_proj <- st_transform(sites_sf, laea_na)
dfna <- dfna %>%
  mutate(x = st_coordinates(sites_proj)[, 1],
         y = st_coordinates(sites_proj)[, 2])

# Analysis data ---------------------------------------------------------------#

index <- "leaf"
years <- 1991:2024

dat <- dfna %>%
  filter(year %in% years) %>%
  mutate(yrc = (year - mean(year)) / 10) %>%   # centered, per decade
  mutate(station = factor(station)) %>%
  select(ghcn, station, lat, lon, x, y, year, yrc, contains(index)) %>%
  arrange(station, year)

yr_mean <- mean(dat$year)

# Mesh ------------------------------------------------------------------------#

inner    <- 300
cut_frac <- 0.5
loc <- as.matrix(dat[, c("x", "y")])

inla_mesh <- fm_mesh_2d_inla(
  loc      = loc,
  max.edge = c(inner, inner * 4),
  offset   = c(inner, inner * 3),
  cutoff   = inner * cut_frac
)
mesh <- make_mesh(dat, c("x", "y"), mesh = inla_mesh)

# Fit (or reuse) --------------------------------------------------------------#

if (use_saved_fit && file.exists(fit_file)) {
  fit1 <- readRDS(fit_file)
  # Guard against a stale fit: mesh and data must match this script's
  stopifnot(fit1$spde$mesh$n == mesh$mesh$n, nrow(fit1$data) == nrow(dat))
} else {
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
  saveRDS(fit1, fit_file)
}
fit_sanity <- sanity(fit1)
fit_fixed  <- tidy(fit1, effects = "fixed", conf.int = TRUE)
fit_ranpar <- tidy(fit1, effects = "ran_par", conf.int = TRUE)

# Land prediction grid (50 km, within 300 km of a station) --------------------#

sites  <- st_as_sf(dat[!duplicated(dat$station), c("x", "y")],
                   coords = c("x", "y"))
domain <- st_union(st_buffer(sites, 300))
g <- st_make_grid(domain, cellsize = 50, what = "centers")
g <- g[lengths(st_intersects(g, domain)) > 0]
xy <- st_coordinates(g)

pred_grid <- data.frame(x = xy[, 1], y = xy[, 2])

world <- ne_countries(scale = "medium", returnclass = "sf")
sf_use_s2(FALSE)
north <- st_crop(world, xmin = -180, xmax = -50, ymin = 20, ymax = 75)
north <- st_transform(north, laea_na)
land  <- st_union(st_make_valid(north))
sf_use_s2(TRUE)

pg_sf <- st_as_sf(pred_grid, coords = c("x", "y"), crs = laea_na, remove = FALSE)
lg <- pred_grid[lengths(st_intersects(pg_sf, land)) > 0, c("x", "y")]
n_land <- nrow(lg)

# Parameters, joint precision, projection -------------------------------------#

par_hat <- fit1$tmb_obj$env$last.par.best
pn      <- names(par_hat)
Q       <- fit1$sd_report$jointPrecision

n_v   <- mesh$mesh$n
n_t   <- length(years)
yrc_t <- (years - yr_mean) / 10

ib <- which(pn == "b_j")         # (Intercept), yrc
io <- which(pn == "omega_s")     # spatial intercept field
iz <- which(pn == "zeta_s")      # spatially varying yrc field
ie <- which(pn == "epsilon_st")  # yearly anomaly fields (vertex fastest)

# Cheap structural guards (full checks are in the validation script)
stopifnot(!is.null(Q), nrow(Q) == length(par_hat),
          length(ib) == 2, length(io) == n_v, length(iz) == n_v,
          length(ie) == n_v * n_t)

A    <- fmesher::fm_basis(mesh$mesh, loc = as.matrix(lg[, c("x", "y")]))
abar <- Matrix::colMeans(A)                                        # n_v
Kt   <- kronecker(Diagonal(n_t), Matrix(abar, nrow = 1, sparse = TRUE))

# Point-estimate local trend --------------------------------------------------#

lg$trend_pt <- par_hat[ib[2]] + as.vector(A %*% par_hat[iz])

# Joint draws, in batches ------------------------------------------------------#

L <- Cholesky(Q, super = TRUE)
rm(Q); gc()

draw_batch <- function(n) {
  z <- matrix(rnorm(length(par_hat) * n), ncol = n)
  z <- solve(L, z, system = "Lt")
  z <- solve(L, z, system = "Pt")
  par_hat + as.matrix(z)
}

tr        <- matrix(NA_real_, n_land, nsim)
ann_draws <- matrix(NA_real_, n_t, nsim, dimnames = list(years, NULL))

set.seed(seed)
for (k in seq(1, nsim, by = batch)) {
  j <- k:min(k + batch - 1, nsim)
  d <- draw_batch(length(j))

  b0   <- d[ib[1], ]
  b1   <- d[ib[2], ]
  zeta <- d[iz, , drop = FALSE]

  # Cell trends: b1 + zeta(s)
  tr[, j] <- sweep(as.matrix(A %*% zeta), 2, b1, "+")

  # Land-averaged annual values
  om_bar  <- as.vector(abar %*% d[io, , drop = FALSE])
  ze_bar  <- as.vector(abar %*% zeta)
  eps_bar <- as.matrix(Kt %*% d[ie, , drop = FALSE])   # n_t x batch

  ann_draws[, j] <- outer(rep(1, n_t), b0 + om_bar) +
    outer(yrc_t, b1 + ze_bar) + eps_bar

  rm(d, zeta, eps_bar); gc()
}
rm(L); gc()

# Hemisphere-wide summaries ---------------------------------------------------#

underlying_trend <- colMeans(tr)   # per draw, days/decade

ann_sum <- data.frame(
  year = years,
  est  = rowMeans(ann_draws),
  lwr  = apply(ann_draws, 1, quantile, 0.025),
  upr  = apply(ann_draws, 1, quantile, 0.975)
)
w <- (yrc_t - mean(yrc_t)) / sum((yrc_t - mean(yrc_t))^2)
realized_trend <- as.vector(crossprod(w, ann_draws))   # per draw, days/decade

# Station info for maps and methods --------------------------------------------#

stations  <- dat %>% distinct(station, lat, lon, x, y)
n_by_year <- dat %>% filter(!is.na(leaf)) %>% count(year, name = "n_stations")

# Save -------------------------------------------------------------------------#

saveRDS(list(lg = lg,   # cell coords + trend_pt; rows match tr
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
             n_land = n_land, cellsize_km = 50, n_vertices = n_v,
             nsim = nsim, batch = batch, seed = seed,
             session = sessionInfo()),
        out_file)
