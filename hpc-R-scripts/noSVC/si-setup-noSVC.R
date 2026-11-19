# SI-x indices, Northern Hemisphere 1951-2024: SHARED SETUP, NO-SVC MODEL
#
# Sourced by si-1-fit-noSVC.R, si-2-validate-noSVC.R and si-3-run-noSVC.R.
# Runs the standard setup (si-setup.R: same data, mesh, land grid, helpers)
# and then renames the output files with "-noSVC", so these runs never
# overwrite the files from the model with the spatially varying trend.
#
# Creates everything si-setup.R creates; run_name, fit_file, validation_file
# and results_file point to the -noSVC versions, e.g.
#   output/si-leaf-1951-2024-northhemis-noSVC-fit.rds

source("si-setup.R")

run_name        <- paste0(run_name, "-noSVC")
fit_file        <- file.path("output", paste0(run_name, "-fit.rds"))
validation_file <- file.path("output", paste0(run_name, "-validation.txt"))
results_file    <- file.path("output", paste0(run_name, "-results.rds"))

stamp("No-SVC model | output files:", run_name)
