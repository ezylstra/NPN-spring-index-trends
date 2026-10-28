# SI-x leaf index, North America 1991-2024: VALIDATION RUN
#
# Fits the model, saves it, and checks that the low-memory draw approach
# (joint parameter draws projected with a sparse A matrix) reproduces sdmTMB's
# own predictions. Run this once; if every check passes, use
# si-leaf-1991-2024-northamerica-run.R for the production analysis.
#
# Checks:
#   1. Model years match `years`, and parameter blocks have the expected sizes
#   2. Point estimates: manual reconstruction vs predict() for several years,
#      overall and component by component (omega_s, zeta_s_yrc, epsilon_st)
#   3. Point-estimate local trend vs predict()'s zeta_s_yrc + fixed slope
#   4. Uncertainty: per-row mean and SD of manual draws vs predict(nsim = )
#      on a small subset of cells (should agree within Monte Carlo error)
#
# Outputs:
#   output/si-leaf-1991-2024-northamerica-fit.rds          (reusable fit)
#   output/si-leaf-1991-2024-northamerica-validation.rds   (all check results)
#   output/si-leaf-1991-2024-northamerica-validation.txt   (readable log)

library(dplyr)
library(stringr)
library(tidyr)
library(sf)
library(sdmTMB)
library(fmesher)
library(Matrix)
library(rnaturalearth)
library(lme4)

dir.create("output", showWarnings = FALSE)

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
cat("Mesh vertices:", mesh$mesh$n, "\n")

# Fit -------------------------------------------------------------------------#

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
fit1 <- readRDS("output/si-leaf-1991-2024-northamerica-fit.rds")
fit_sanity <- sanity(fit1)
# saveRDS(fit1, "output/si-leaf-1991-2024-northamerica-fit.rds")

# Land prediction grid (50 km, within 300 km of a station) --------------------#

sites  <- st_as_sf(dat[!duplicated(dat$station), c("x", "y")],
                   coords = c("x", "y"))
domain <- st_union(st_buffer(sites, 300))
g <- st_make_grid(domain, cellsize = 50, what = "centers")
g <- g[lengths(st_intersects(g, domain)) > 0]
xy <- st_coordinates(g)

pred_grid <- data.frame(
  x = xy[, 1],
  y = xy[, 2],
  station = factor(levels(dat$station)[1], levels = levels(dat$station))
)

world <- ne_countries(scale = "medium", returnclass = "sf")
sf_use_s2(FALSE)
north <- st_crop(world, xmin = -180, xmax = -50, ymin = 20, ymax = 75)
north <- st_transform(north, laea_na)
land  <- st_union(st_make_valid(north))
sf_use_s2(TRUE)

pg_sf <- st_as_sf(pred_grid, coords = c("x", "y"), crs = laea_na, remove = FALSE)
pred_grid$land <- lengths(st_intersects(pg_sf, land)) > 0
lg <- pred_grid[pred_grid$land, c("x", "y", "station")]
n_land <- nrow(lg)
cat("Land cells:", n_land, "\n")

# Parameter bookkeeping -------------------------------------------------------#

par_hat <- fit1$tmb_obj$env$last.par.best
pn      <- names(par_hat)
Q       <- fit1$sd_report$jointPrecision

n_v   <- mesh$mesh$n
n_t   <- length(years)
yrc_t <- (years - yr_mean) / 10

ib <- which(pn == "b_j")
io <- which(pn == "omega_s")
iz <- which(pn == "zeta_s")
ie <- which(pn == "epsilon_st")

A <- fmesher::fm_basis(mesh$mesh, loc = as.matrix(lg[, c("x", "y")]))

# CHECK 1: structure ----------------------------------------------------------#

time_lu_years <- tryCatch(as.numeric(fit1$time_lu$time_from_data),
                          error = function(e) NA)
check1 <- list(
  joint_precision_present = !is.null(Q),
  Q_dim_matches_par       = !is.null(Q) && nrow(Q) == length(par_hat),
  n_b_j                   = length(ib),
  b_j_names_fixef         = names(fixef(fit1)),
  n_omega_s               = length(io),
  n_zeta_s                = length(iz),
  n_epsilon_st            = length(ie),
  expected_epsilon_st     = n_v * n_t,
  model_years             = time_lu_years,
  years_match             = isTRUE(all.equal(time_lu_years, as.numeric(years))),
  A_dim                   = dim(A)
)
check1$pass <- with(check1,
  joint_precision_present && Q_dim_matches_par && n_b_j == 2 &&
  n_omega_s == n_v && n_zeta_s == n_v &&
  n_epsilon_st == expected_epsilon_st && years_match &&
  A_dim[1] == n_land && A_dim[2] == n_v)

# CHECK 2: point estimates for several years, by component --------------------#

chk_years <- years[c(1, ceiling(n_t / 2), n_t)]

eps_hat <- function(t) par_hat[ie[(t - 1) * n_v + seq_len(n_v)]]

check2 <- lapply(chk_years, function(yy) {
  t  <- match(yy, years)
  yc <- yrc_t[t]
  p  <- predict(fit1, newdata = transform(lg, year = yy, yrc = yc),
                re_form_iid = NA)
  m_omega <- as.vector(A %*% par_hat[io])
  m_zeta  <- as.vector(yc * (A %*% par_hat[iz]))
  m_eps   <- as.vector(A %*% eps_hat(t))
  m_est   <- par_hat[ib[1]] + par_hat[ib[2]] * yc + m_omega + m_zeta + m_eps
  data.frame(
    year       = yy,
    est        = max(abs(p$est - m_est)),
    omega_s    = if ("omega_s" %in% names(p)) max(abs(p$omega_s - m_omega)) else NA,
    zeta_s_yrc = if ("zeta_s_yrc" %in% names(p)) max(abs(p$zeta_s_yrc - m_zeta)) else NA,
    epsilon_st = if ("epsilon_st" %in% names(p)) max(abs(p$epsilon_st - m_eps)) else NA
  )
})
check2 <- do.call(rbind, check2)
check2_pass <- all(check2$est < 1e-6)

# CHECK 3: point-estimate local trend -----------------------------------------#

p1 <- predict(fit1, newdata = transform(lg, year = years[1], yrc = 1),
              re_form_iid = NA)
# predict's zeta_s_yrc = yrc * zeta(s); at yrc = 1 it is zeta(s)
trend_pred   <- fixef(fit1)["yrc"] + p1$zeta_s_yrc
trend_manual <- par_hat[ib[2]] + as.vector(A %*% par_hat[iz])
check3 <- max(abs(trend_pred - trend_manual))
check3_pass <- check3 < 1e-6

# CHECK 4: uncertainty vs predict(nsim) on a cell subset ----------------------#

n_sub  <- 200
nsim_v <- 500
set.seed(99)
sub <- sort(sample(n_land, n_sub))
nd_sub <- as.data.frame(tidyr::expand_grid(lg[sub, ], year = years)) %>%
  mutate(yrc = (year - yr_mean) / 10)

set.seed(1)
sim_pred <- predict(fit1, newdata = nd_sub, nsim = nsim_v, re_form_iid = NA)

L <- Cholesky(Q, super = TRUE)
draw_batch <- function(n) {
  z <- matrix(rnorm(length(par_hat) * n), ncol = n)
  z <- solve(L, z, system = "Lt")
  z <- solve(L, z, system = "Pt")
  par_hat + as.matrix(z)
}
set.seed(2)
d <- draw_batch(nsim_v)

A_sub <- A[sub, , drop = FALSE]
sim_man <- matrix(NA_real_, nrow(nd_sub), nsim_v)
for (t in seq_len(n_t)) {
  rows <- which(nd_sub$year == years[t])     # same cell order as sub
  yc   <- yrc_t[t]
  eta  <- A_sub %*% d[io, , drop = FALSE] +
    yc * (A_sub %*% d[iz, , drop = FALSE]) +
    A_sub %*% d[ie[(t - 1) * n_v + seq_len(n_v)], , drop = FALSE]
  eta  <- sweep(as.matrix(eta), 2, d[ib[1], ] + d[ib[2], ] * yc, "+")
  sim_man[rows, ] <- eta
}

sd_pred <- apply(sim_pred, 1, sd)
sd_man  <- apply(sim_man, 1, sd)
mean_pred <- rowMeans(sim_pred)
mean_man  <- rowMeans(sim_man)

check4 <- list(
  sd_ratio_median   = median(sd_man / sd_pred),
  sd_ratio_range    = range(sd_man / sd_pred),
  sd_correlation    = cor(sd_man, sd_pred),
  mean_diff_over_sd = summary(abs(mean_man - mean_pred) / sd_pred)
)
# With 500 draws, SD Monte Carlo error is ~3%; medians near 1 and a high
# correlation indicate the same distribution
check4_pass <- abs(check4$sd_ratio_median - 1) < 0.05 &&
  check4$sd_correlation > 0.95
rm(d, sim_pred, sim_man, L); gc()

# Report ----------------------------------------------------------------------#

all_pass <- check1$pass && check2_pass && check3_pass && check4_pass

validation <- list(
  all_pass = all_pass,
  check1_structure = check1,
  check2_point_estimates = check2, check2_pass = check2_pass,
  check3_trend_max_abs_diff = check3, check3_pass = check3_pass,
  check4_uncertainty = check4, check4_pass = check4_pass,
  fit_sanity = fit_sanity,
  fit_fixed  = tidy(fit1, effects = "fixed", conf.int = TRUE),
  fit_ranpar = tidy(fit1, effects = "ran_par", conf.int = TRUE),
  n_vertices = n_v, n_land = n_land,
  session = sessionInfo()
)
saveRDS(validation, "output/si-leaf-1991-2024-northamerica-validation.rds")

sink("output/si-leaf-1991-2024-northamerica-validation.txt")
cat("ALL CHECKS PASS:", all_pass, "\n\n")
cat("Check 1 (structure):", check1$pass, "\n"); str(check1)
cat("\nCheck 2 (point estimates, max abs diff):", check2_pass, "\n"); print(check2)
cat("\nCheck 3 (local trend, max abs diff):", check3_pass, "\n"); print(check3)
cat("\nCheck 4 (uncertainty vs predict(nsim)):", check4_pass, "\n"); str(check4)
cat("\nSanity:\n"); print(fit_sanity)
sink()

cat("\nALL CHECKS PASS:", all_pass, "\n")
