# SI-x indices, Northern Hemisphere 1951-2024: STEP 3, PRODUCTION RUN (NO-SVC MODEL)
#
# Usage: Rscript si-3-run-noSVC.R <index>     (after si-2-validate-noSVC.R passes)
#
# Same products, names and file structure as si-3-run.R, so
# visualizing-results.R works unchanged (apart from the file name). The
# difference: this model has no spatially varying trend, so the UNDERLYING
# local trend is the fixed slope b1 in every cell:
#   tr               n_land x nsim   underlying trend draws (each column = b1,
#                                    the same in every cell)
#   lg$trend_pt      point-estimate underlying trend (= b1 in every cell)
#   ann_draws        years x nsim    land-averaged annual index draws
#   ann_sum          annual means and 95% intervals
#   underlying_trend, realized_trend (per draw, same draws)
#   tr_real          list of n_land x nsim matrices, one per period in
#                    `periods`: REALIZED local trend draws (days/decade), i.e.
#                    the OLS slope through each cell's modeled annual values
#                    (including the yearly anomaly fields) over that period.
#                    These still vary across space.
#   realized_by_period  list of per-draw land-averaged realized trends
# Peak memory during simulation is set by `batch`, not by cells x years x nsim.
#
# Output: output/si-<index>-1951-2024-northhemis-noSVC-results.rds

source("si-setup-noSVC.R")   # dat, years, yr_mean, mesh, laea_km, file names, helpers

# Settings --------------------------------------------------------------------#

nsim  <- 1000   # consider 2000 for lfsr tails
batch <- 100    # draws held in memory at once; lower if memory is tight
seed  <- 1234

# Periods for realized trends (local and land-averaged). Each must be a run of
# years inside `years`. Adding a period does not change the draws themselves.
periods <- list(full = 1951:2024,
                recent = 1991:2024,
                recent30 = 1991:2020)

# Load fit --------------------------------------------------------------------#

fit1 <- readRDS(fit_file)
# Guard against a stale fit: mesh and data must match the setup script
stopifnot(fit1$spde$mesh$n == mesh$mesh$n, nrow(fit1$data) == nrow(dat))
fit_sanity <- sanity(fit1)
fit_fixed  <- tidy(fit1, effects = "fixed", conf.int = TRUE)
fit_ranpar <- tidy(fit1, effects = "ran_par", conf.int = TRUE)
stamp("Fit loaded:", fit_file)

# Land prediction grid (same function as validation) --------------------------#

lg <- make_land_grid(dat)
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
iz <- which(pn == "zeta_s")      # spatially varying yrc field: absent here
ie <- which(pn == "epsilon_st")  # yearly anomaly fields (vertex fastest)

# Cheap structural guards (full checks are in the validation script)
stopifnot(!is.null(Q), nrow(Q) == length(par_hat),
          length(ib) == 2, length(io) == n_v, length(iz) == 0,
          length(ie) == n_v * n_t)

A    <- fmesher::fm_basis(mesh$mesh, loc = as.matrix(lg[, c("x", "y")]))
abar <- Matrix::colMeans(A)                                        # n_v
Kt   <- kronecker(Diagonal(n_t), Matrix(abar, nrow = 1, sparse = TRUE))

rm(fit1); gc()   # everything needed below has been extracted

# Point-estimate underlying trend: the fixed slope, same in every cell --------#

lg$trend_pt <- rep(unname(par_hat[ib[2]]), n_land)

# Joint draws, in batches ------------------------------------------------------#

stamp("Cholesky of joint precision")
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

# Realized trends per period: OLS weights on decade-scaled time. Because the
# weights sum to 0, the intercept and omega fields drop out, and because
# sum(w * yrc) = 1, a cell's realized slope is
#   b1 + sum_t w_t * epsilon_t(s)  = tr + A %*% (weighted eps)
stopifnot(all(unlist(periods) %in% years))
per_idx <- lapply(periods, function(p) match(p, years))
per_w   <- lapply(per_idx, function(i) {
  x <- yrc_t[i] - mean(yrc_t[i])
  x / sum(x^2)
})
tr_real <- lapply(periods, function(p) matrix(NA_real_, n_land, nsim))

set.seed(seed)
for (k in seq(1, nsim, by = batch)) {
  j <- k:min(k + batch - 1, nsim)
  d <- draw_batch(length(j))

  b0 <- d[ib[1], ]
  b1 <- d[ib[2], ]

  # Underlying trend: b1, the same in every cell
  tr[, j] <- matrix(b1, n_land, length(j), byrow = TRUE)

  # Realized local trends: add the trend in the yearly anomaly fields
  for (p in names(periods)) {
    ew <- matrix(0, n_v, length(j))
    for (m in seq_along(per_idx[[p]])) {
      t  <- per_idx[[p]][m]
      ew <- ew + per_w[[p]][m] * d[ie[(t - 1) * n_v + seq_len(n_v)], , drop = FALSE]
    }
    tr_real[[p]][, j] <- tr[, j] + as.matrix(A %*% ew)
  }

  # Land-averaged annual values
  om_bar  <- as.vector(abar %*% d[io, , drop = FALSE])
  eps_bar <- as.matrix(Kt %*% d[ie, , drop = FALSE])   # n_t x batch

  ann_draws[, j] <- outer(rep(1, n_t), b0 + om_bar) +
    outer(yrc_t, b1) + eps_bar

  rm(d, eps_bar); gc()
  stamp("Draws", max(j), "of", nsim, "| peak mem",
        round(mem_gb()[["peak"]], 1), "GB")
}
rm(L); gc()

# Hemisphere-wide summaries ---------------------------------------------------#

underlying_trend <- colMeans(tr)   # per draw, days/decade (= b1)

ann_sum <- data.frame(
  year = years,
  est  = rowMeans(ann_draws),
  lwr  = apply(ann_draws, 1, quantile, 0.025),
  upr  = apply(ann_draws, 1, quantile, 0.975)
)

realized_by_period <- lapply(names(periods), function(p) {
  as.vector(crossprod(per_w[[p]], ann_draws[per_idx[[p]], , drop = FALSE]))
})
names(realized_by_period) <- names(periods)
realized_trend <- realized_by_period$full   # per draw, days/decade

# Consistency check: the land average of the realized local trends must equal
# the realized trend of the land-averaged annual values (same draws)
for (p in names(periods)) {
  dev <- max(abs(colMeans(tr_real[[p]]) - realized_by_period[[p]]))
  stamp("Check realized trend,", p, "| max abs diff:", signif(dev, 3))
  stopifnot(dev < 1e-6)
}

# Station info for maps and methods --------------------------------------------#

stations  <- dat %>% distinct(station, lat, lon, x, y)
n_by_year <- dat %>% count(year, name = "n_stations")

# Save -------------------------------------------------------------------------#

saveRDS(list(index = index,
             model = "noSVC",   # no spatially varying trend
             lg = lg,   # cell coords + trend_pt; rows match tr
             tr = tr,
             ann_draws = ann_draws,
             ann_sum = ann_sum,
             underlying_trend = underlying_trend,
             realized_trend = realized_trend,
             periods = periods,
             tr_real = tr_real,   # list by period; rows match lg
             realized_by_period = realized_by_period,
             stations = stations,
             n_by_year = n_by_year,
             fit_sanity = fit_sanity,
             fit_fixed = fit_fixed,
             fit_ranpar = fit_ranpar,
             years = years, yr_mean = yr_mean, laea_km = laea_km,
             n_land = n_land, cellsize_km = 50, n_vertices = n_v,
             nsim = nsim, batch = batch, seed = seed,
             session = sessionInfo()),
        results_file)
stamp("Results saved:", results_file,
      sprintf("(%.0f MB)", file.size(results_file) / 1024^2))
