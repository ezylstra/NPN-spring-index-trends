# SI-x indices, Northern Hemisphere 1991-2024: STEP 3, PRODUCTION RUN
#
# Usage: Rscript si-3-run-1991.R <index>     (after si-2-validate-1991.R passes)
#
# Loads the saved fit and builds all products from joint parameter draws
# projected with a sparse A matrix:
#   tr               n_land x nsim   underlying local trend draws (days/decade)
#   lg$trend_pt      point-estimate local trend
#   ann_draws        years x nsim    land-averaged annual index draws
#   ann_sum          annual means and 95% intervals
#   underlying_trend, realized_trend (per draw, same draws)
#   tr_real          n_land x nsim   REALIZED local trend draws (days/decade),
#                    i.e. the OLS slope through each cell's modeled annual values
#                    (including the yearly anomaly fields) over 1991-2024
# Realized trends are for the full period only, so (unlike si-3-run-1955.R)
# there is no list of periods and tr_real is a single matrix.
# Peak memory during simulation is set by `batch`, not by cells x years x nsim.
#
# Output: output/si-<index>-1991-2024-northhemis-results.rds

source("si-setup-1991.R")   # dat, years, yr_mean, mesh, laea_km, file names, helpers

# Settings --------------------------------------------------------------------#

nsim  <- 1000   # consider 2000 for lfsr tails
batch <- 100    # draws held in memory at once; lower if memory is tight
seed  <- 1234

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
iz <- which(pn == "zeta_s")      # spatially varying yrc field
ie <- which(pn == "epsilon_st")  # yearly anomaly fields (vertex fastest)

# Cheap structural guards (full checks are in the validation script)
stopifnot(!is.null(Q), nrow(Q) == length(par_hat),
          length(ib) == 2, length(io) == n_v, length(iz) == n_v,
          length(ie) == n_v * n_t)

A    <- fmesher::fm_basis(mesh$mesh, loc = as.matrix(lg[, c("x", "y")]))
abar <- Matrix::colMeans(A)                                        # n_v
Kt   <- kronecker(Diagonal(n_t), Matrix(abar, nrow = 1, sparse = TRUE))

rm(fit1); gc()   # everything needed below has been extracted

# Point-estimate local trend --------------------------------------------------#

lg$trend_pt <- par_hat[ib[2]] + as.vector(A %*% par_hat[iz])

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
tr_real   <- matrix(NA_real_, n_land, nsim)
ann_draws <- matrix(NA_real_, n_t, nsim, dimnames = list(years, NULL))

# Realized trends: OLS weights on decade-scaled time. Because the weights sum
# to 0, the intercept and omega fields drop out, and because
# sum(w * yrc) = 1, a cell's realized slope is
#   b1 + zeta(s) + sum_t w_t * epsilon_t(s)  = tr + A %*% (weighted eps)
# (yrc_t is re-centered because yr_mean is the mean over rows, not over years)
x_t <- yrc_t - mean(yrc_t)
w_t <- x_t / sum(x_t^2)

set.seed(seed)
for (k in seq(1, nsim, by = batch)) {
  j <- k:min(k + batch - 1, nsim)
  d <- draw_batch(length(j))

  b0   <- d[ib[1], ]
  b1   <- d[ib[2], ]
  zeta <- d[iz, , drop = FALSE]

  # Underlying local trends: b1 + zeta(s)
  tr[, j] <- sweep(as.matrix(A %*% zeta), 2, b1, "+")

  # Realized local trends: add the trend in the yearly anomaly fields
  ew <- matrix(0, n_v, length(j))
  for (t in seq_len(n_t)) {
    ew <- ew + w_t[t] * d[ie[(t - 1) * n_v + seq_len(n_v)], , drop = FALSE]
  }
  tr_real[, j] <- tr[, j] + as.matrix(A %*% ew)

  # Land-averaged annual values
  om_bar  <- as.vector(abar %*% d[io, , drop = FALSE])
  ze_bar  <- as.vector(abar %*% zeta)
  eps_bar <- as.matrix(Kt %*% d[ie, , drop = FALSE])   # n_t x batch

  ann_draws[, j] <- outer(rep(1, n_t), b0 + om_bar) +
    outer(yrc_t, b1 + ze_bar) + eps_bar

  rm(d, zeta, eps_bar); gc()
  stamp("Draws", max(j), "of", nsim, "| peak mem",
        round(mem_gb()[["peak"]], 1), "GB")
}
rm(L); gc()

# Hemisphere-wide summaries ---------------------------------------------------#

underlying_trend <- colMeans(tr)                              # per draw, days/decade
realized_trend   <- as.vector(crossprod(w_t, ann_draws))      # per draw, days/decade

ann_sum <- data.frame(
  year = years,
  est  = rowMeans(ann_draws),
  lwr  = apply(ann_draws, 1, quantile, 0.025),
  upr  = apply(ann_draws, 1, quantile, 0.975)
)

# Consistency check: the land average of the realized local trends must equal
# the realized trend of the land-averaged annual values (same draws)
dev <- max(abs(colMeans(tr_real) - realized_trend))
stamp("Check realized trend | max abs diff:", signif(dev, 3))
stopifnot(dev < 1e-6)

# Station info for maps and methods --------------------------------------------#

stations  <- dat %>% distinct(station, lat, lon, x, y)
n_by_year <- dat %>% count(year, name = "n_stations")

# Save -------------------------------------------------------------------------#

saveRDS(list(index = index,
             lg = lg,   # cell coords + trend_pt; rows match tr and tr_real
             tr = tr,
             ann_draws = ann_draws,
             ann_sum = ann_sum,
             underlying_trend = underlying_trend,
             realized_trend = realized_trend,
             tr_real = tr_real,
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
