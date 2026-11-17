# SI-x indices, Northern Hemisphere 1951-2024: STEP 1, FIT
#
# Usage: Rscript si-1-fit.R <index>     (leaf, bloom, frz or damage)
# Output: output/si-<index>-1951-2024-northhemis-fit.rds
#
# To watch optimizer progress in the log (e.g. if a fit seems stuck), add
# silent = FALSE to sdmTMB(); each "outer mgc" line is one evaluation.

source("si-setup.R")

form <- as.formula(paste(index, "~ yrc + (1 | station)"))

stamp("Fitting started:", deparse(form))
t_fit <- system.time(
  fit1 <- sdmTMB(
    form,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "on",
    spatiotemporal = "iid",
    spatial_varying = ~ 0 + yrc,
    family = gaussian()
  )
)
stamp(sprintf("Fitting done in %.1f hours | peak mem %.1f GB",
              t_fit[["elapsed"]] / 3600, mem_gb()[["peak"]]))

# Save immediately, before anything else can fail
saveRDS(fit1, fit_file)
stamp("Fit saved:", fit_file)

print(sanity(fit1))
print(fit1)
stamp("Done")
