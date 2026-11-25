# SI-x indices, Northern Hemisphere 1991-2024: SHARED SETUP
#
# Sourced by si-1-fit-1991.R, si-2-validate-1991.R and si-3-run-1991.R so all three use the
# same data, mesh, land grid and file names. The index comes from the first
# command-line argument, e.g.  Rscript si-1-fit-1991.R leaf  (default: bloom).
#
# Creates: index, years, yr_mean, dat, mesh, laea_km, run_name, fit_file,
#          validation_file, results_file, stamp(), mem_gb(), make_land_grid()

library(dplyr)
library(stringr)
library(tidyr)
library(sf)
library(sdmTMB)
library(fmesher)
library(Matrix)
library(rnaturalearth)

dir.create("output", showWarnings = FALSE)

# Helpers ---------------------------------------------------------------------#

# Timestamped progress line for the Slurm log
stamp <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
  flush.console()
}

# Current and peak memory of this R process (GB), read from Linux
mem_gb <- function() {
  s <- readLines("/proc/self/status")
  get <- function(key) {
    v <- grep(paste0("^", key, ":"), s, value = TRUE)
    as.numeric(gsub("[^0-9]", "", v)) / 1024^2   # kB -> GB
  }
  c(current = get("VmRSS"), peak = get("VmHWM"))
}

# Index, years, file names ----------------------------------------------------#

args    <- commandArgs(trailingOnly = TRUE)
index   <- if (length(args) >= 1) args[1] else "bloom"
indices <- c("leaf", "bloom", "frz", "damage")
if (!index %in% indices) {
  stop("Index must be one of: ", paste(indices, collapse = ", "))
}
years <- 1991:2024

run_name        <- sprintf("si-%s-%d-%d-northhemis", index, min(years), max(years))
fit_file        <- file.path("output", paste0(run_name, "-fit.rds"))
validation_file <- file.path("output", paste0(run_name, "-validation.txt"))
results_file    <- file.path("output", paste0(run_name, "-results.rds"))

stamp("Setup started:", run_name, "| node:", Sys.info()[["nodename"]])

# OpenBLAS check: a 2000 x 2000 cross-product takes ~0.2 s with OpenBLAS and
# ~4-5 s with R's reference BLAS. Uses fixed values so the RNG is not touched.
blas_sec <- system.time(
  crossprod(matrix(seq_len(4e6) / 4e6, 2000))
)[["elapsed"]]
stamp(sprintf("BLAS speed test: %.2f s", blas_sec))
if (blas_sec > 1.5) {
  warning("BLAS looks slow; OpenBLAS may not be active (check LD_PRELOAD). ",
          "The fit will run, but ~14x slower.", immediate. = TRUE)
}

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

# Project site coordinates (km) -----------------------------------------------#

laea_km <- "+proj=laea +lat_0=90 +lon_0=0 +units=km"

sites_sf   <- st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
sites_proj <- st_transform(sites_sf, laea_km)
df <- df %>%
  mutate(x = st_coordinates(sites_proj)[, 1],
         y = st_coordinates(sites_proj)[, 2])

# Analysis data ---------------------------------------------------------------#
# Rows with a missing index value are dropped (none for bloom)

dat <- df %>%
  filter(year %in% years, !is.na(.data[[index]])) %>%
  mutate(yrc = (year - mean(year)) / 10) %>%   # centered, per decade
  mutate(station = factor(station)) %>%
  select(ghcn, station, lat, lon, x, y, year, yrc, contains(index)) %>%
  arrange(station, year)

yr_mean <- mean(dat$year)

# Mesh ------------------------------------------------------------------------#

inner    <- 500
cut_frac <- 0.3
loc <- as.matrix(dat[, c("x", "y")])

inla_mesh <- fm_mesh_2d_inla(
  loc      = loc,
  max.edge = c(inner, inner * 4),
  offset   = c(inner, inner * 3),
  cutoff   = inner * cut_frac
)
mesh <- make_mesh(dat, c("x", "y"), mesh = inla_mesh)

stamp("Setup done | rows:", nrow(dat),
      "| stations:", nlevels(dat$station),
      "| mesh vertices:", mesh$mesh$n,
      "| random effects (approx):", mesh$mesh$n * (2 + length(years)) +
        nlevels(dat$station))

# Land prediction grid --------------------------------------------------------#
# Cell centers on a regular grid within buffer_km of a station, kept if they
# fall on land. Used by both validation and the production run.

make_land_grid <- function(dat, buffer_km = 300, cellsize_km = 50) {
  sites  <- st_as_sf(dat[!duplicated(dat$station), c("x", "y")],
                     coords = c("x", "y"))
  domain <- st_union(st_buffer(sites, buffer_km))
  g  <- st_make_grid(domain, cellsize = cellsize_km, what = "centers")
  g  <- g[lengths(st_intersects(g, domain)) > 0]
  xy <- st_coordinates(g)
  pred_grid <- data.frame(x = xy[, 1], y = xy[, 2])

  world <- ne_countries(scale = "medium", returnclass = "sf")
  sf_use_s2(FALSE)
  north <- st_crop(world, xmin = -180, xmax = 180, ymin = 20, ymax = 90)
  north <- st_transform(north, laea_km)
  land  <- st_union(st_make_valid(north))
  sf_use_s2(TRUE)

  pg_sf <- st_as_sf(pred_grid, coords = c("x", "y"), crs = laea_km,
                    remove = FALSE)
  lg <- pred_grid[lengths(st_intersects(pg_sf, land)) > 0, c("x", "y")]
  rownames(lg) <- NULL
  stamp("Land cells:", nrow(lg))
  lg
}
