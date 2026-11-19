# SI-x indices, Northern Hemisphere 1951-2024: STEP 1, FIT (NO-SVC MODEL)
#
# Usage: Rscript si-1-fit-noSVC.R <index>     (leaf, bloom, frz or damage)
# Output: output/si-<index>-1951-2024-northhemis-noSVC-fit.rds
#
# Same model as si-1-fit.R (fixed trend, station intercepts, spatial field,
# iid yearly anomaly fields) but WITHOUT the spatially varying trend
# (spatial_varying = ~ 0 + yrc). The underlying trend is therefore the same
# everywhere (the fixed slope); spatial variation in trends comes only from
# the realized trends (spatial field + yearly anomaly fields).
#
# Use for indices where sigma_Z collapses to ~0 and the joint precision
# matrix is not positive definite (leaf, 1951-2024).
#
# Exits with an error if the Hessian is not positive definite or any
# standard error is NA, so the dependent validation job will not start.
#
# To watch optimizer progress in the log (e.g. if a fit seems stuck), add
# silent = FALSE to sdmTMB(); each "outer mgc" line is one evaluation.

source("si-setup-noSVC.R")

form <- as.formula(paste(index, "~ yrc + (1 | station)"))

stamp("Fitting started (no SVC):", deparse(form))
t_fit <- system.time(
  fit1 <- sdmTMB(
    form,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "on",
    spatiotemporal = "iid",
    family = gaussian()
  )
)
stamp(sprintf("Fitting done in %.1f hours | peak mem %.1f GB",
              t_fit[["elapsed"]] / 3600, mem_gb()[["peak"]]))

# Save immediately, before anything else can fail
saveRDS(fit1, fit_file)
stamp("Fit saved:", fit_file)

s <- sanity(fit1)
print(s)
print(fit1)

# Stop the chain here if the fit can't support joint draws
if (!isTRUE(s$hessian_ok) || !isTRUE(s$se_na_ok)) {
  stamp("Hessian not positive definite or NA standard errors: stopping")
  quit(save = "no", status = 1)
}
stamp("Done")
