# SI-x indices, Northern Hemisphere 1951-2024: STEP 2, VALIDATE (NO-SVC MODEL)
#
# Usage: Rscript si-2-validate-noSVC.R <index>     (after si-1-fit-noSVC.R)
#
# Checks that the low-memory draw approach used in si-3-run-noSVC.R (joint
# parameter draws projected with a sparse A matrix) reproduces sdmTMB's own
# predictions for this fit:
#   1. Model years match `years`, parameter blocks have the expected sizes,
#      and there is no zeta_s (spatially varying trend) block
#   2. Point estimates: manual reconstruction vs predict() for several years,
#      overall and component by component (omega_s, epsilon_st)
#   3. Point-estimate underlying trend: predict() at yrc = 1 minus yrc = 0
#      (same year, so the yearly field cancels) equals the fixed slope in
#      every cell
#   4. Uncertainty: per-row mean and SD of manual draws vs predict(nsim = )
#      on a small subset of cells (should agree within Monte Carlo error)
#
# Output: output/si-<index>-1951-2024-northhemis-noSVC-validation.txt
# Exits with an error if any check fails, so a run job submitted with
# --dependency=afterok on this job will not start.

source("si-setup-noSVC.R")   # dat, years, yr_mean, mesh, file names, helpers
library(lme4)                # fixef()

# Load the saved fit ----------------------------------------------------------#

fit1 <- readRDS(fit_file)
stopifnot(fit1$spde$mesh$n == mesh$mesh$n,   # same mesh as the fit
          nrow(fit1$data) == nrow(dat))      # same data as the fit
stamp("Fit loaded:", fit_file)

fit_sanity <- sanity(fit1)

# Land prediction grid --------------------------------------------------------#
# predict() needs every model variable, so add a placeholder station
# (ignored: re_form_iid = NA drops the station intercepts)

lg <- make_land_grid(dat)
lg$station <- factor(levels(dat$station)[1], levels = levels(dat$station))
n_land <- nrow(lg)

# Parameter bookkeeping -------------------------------------------------------#

par_hat <- fit1$tmb_obj$env$last.par.best
pn      <- names(par_hat)
Q       <- fit1$sd_report$jointPrecision

n_v   <- mesh$mesh$n
n_t   <- length(years)
yrc_t <- (years - yr_mean) / 10

ib <- which(pn == "b_j")
io <- which(pn == "omega_s")
iz <- which(pn == "zeta_s")       # should be empty in this model
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
  n_omega_s == n_v && n_zeta_s == 0 &&
  n_epsilon_st == expected_epsilon_st && years_match &&
  A_dim[1] == n_land && A_dim[2] == n_v)
stamp("Check 1 done:", check1$pass)

# CHECK 2: point estimates for several years, by component --------------------#

chk_years <- years[c(1, ceiling(n_t / 2), n_t)]

eps_hat <- function(t) par_hat[ie[(t - 1) * n_v + seq_len(n_v)]]

check2 <- lapply(chk_years, function(yy) {
  t  <- match(yy, years)
  yc <- yrc_t[t]
  p  <- predict(fit1, newdata = transform(lg, year = yy, yrc = yc),
                re_form_iid = NA)
  m_omega <- as.vector(A %*% par_hat[io])
  m_eps   <- as.vector(A %*% eps_hat(t))
  m_est   <- par_hat[ib[1]] + par_hat[ib[2]] * yc + m_omega + m_eps
  data.frame(
    year       = yy,
    est        = max(abs(p$est - m_est)),
    omega_s    = if ("omega_s" %in% names(p)) max(abs(p$omega_s - m_omega)) else NA,
    epsilon_st = if ("epsilon_st" %in% names(p)) max(abs(p$epsilon_st - m_eps)) else NA
  )
})
check2 <- do.call(rbind, check2)
check2_pass <- all(check2$est < 1e-6)
stamp("Check 2 done:", check2_pass)

# CHECK 3: point-estimate underlying trend ------------------------------------#
# Same year (same yearly field), yrc = 1 vs yrc = 0: the difference is the
# trend per decade, which in this model is the fixed slope everywhere

p1 <- predict(fit1, newdata = transform(lg, year = years[1], yrc = 1),
              re_form_iid = NA)
p0 <- predict(fit1, newdata = transform(lg, year = years[1], yrc = 0),
              re_form_iid = NA)
trend_pred   <- p1$est - p0$est
trend_manual <- par_hat[ib[2]]
check3 <- max(abs(trend_pred - trend_manual),
              abs(fixef(fit1)[["yrc"]] - trend_manual))
check3_pass <- check3 < 1e-6
rm(p0, p1)
stamp("Check 3 done:", check3_pass)

# CHECK 4: uncertainty vs predict(nsim) on a cell subset ----------------------#

n_sub  <- 200
nsim_v <- 500
set.seed(99)
sub <- sort(sample(n_land, n_sub))
nd_sub <- as.data.frame(tidyr::expand_grid(lg[sub, ], year = years)) %>%
  mutate(yrc = (year - yr_mean) / 10)

stamp("Check 4: predict(nsim =", nsim_v, ") on", nrow(nd_sub), "rows")
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
    A_sub %*% d[ie[(t - 1) * n_v + seq_len(n_v)], , drop = FALSE]
  eta  <- sweep(as.matrix(eta), 2, d[ib[1], ] + d[ib[2], ] * yc, "+")
  sim_man[rows, ] <- eta
}

sd_pred   <- apply(sim_pred, 1, sd)
sd_man    <- apply(sim_man, 1, sd)
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
stamp("Check 4 done:", check4_pass, "| peak mem",
      round(mem_gb()[["peak"]], 1), "GB")
rm(d, sim_pred, sim_man, L); gc()

# Report ----------------------------------------------------------------------#

all_pass <- check1$pass && check2_pass && check3_pass && check4_pass

sink(validation_file)
cat("Run:", run_name, "\n")
cat("ALL CHECKS PASS:", all_pass, "\n\n")
cat("Check 1 (structure):", check1$pass, "\n"); str(check1)
cat("\nCheck 2 (point estimates, max abs diff):", check2_pass, "\n"); print(check2)
cat("\nCheck 3 (underlying trend = fixed slope, max abs diff):", check3_pass, "\n"); print(check3)
cat("\nCheck 4 (uncertainty vs predict(nsim)):", check4_pass, "\n"); str(check4)
cat("\nSanity:\n"); print(fit_sanity)
sink()

stamp("ALL CHECKS PASS:", all_pass, "| report:", validation_file)
if (!all_pass) quit(save = "no", status = 1)   # blocks the dependent run job
