# SI-x indices, Northern Hemisphere 1955-2024: VISUALIZE RESULTS
#
# Evaluate trends in the length of time between when bloom and leaf indices
# are reached.

library(dplyr)
library(tidyr)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(terra)

# Steps
# 1) Load results from bloom & leaf models
# 2) For each period of interest (1955-2024, 1991-2020, and 1991-2024):
  # 2a) Generate rasters for each index with local trend estimates (realized trends)
  # 2b) Difference the rasters (bloom - leaf) to estimate whether the gap between
  #     the two indexes is compressing or lengthening over time
#
# The gap is bloom DOY minus leaf DOY, so the trend in the gap equals the bloom
# trend minus the leaf trend. Positive = gap lengthening; negative = gap
# shrinking (days/decade).
#
# The two models were fit separately, so draw i of one is independent of
# draw i of the other. Differencing draws pairwise gives valid draws of the
# gap trend (uncertainty from both models is carried through).

# Settings --------------------------------------------------------------------#

indices <- c("bloom", "leaf")
run_names <- setNames(paste0("si-", indices, "-1955-2024-northhemis"), indices)
q <- 0.05   # error rate for significance (sign) calls

# TRUE: one color scale shared by maps of all periods. FALSE: each period gets
# its own scale (see visualizing-results.R).
shared_color_ramp <- FALSE

# 1) Load results -------------------------------------------------------------#

res <- lapply(run_names, function(r) {
  readRDS(paste0("model-fits-results/", r, "-results.rds"))
})

# Settings that must agree between the two runs
for (k in c("periods", "nsim", "laea_km", "cellsize_km", "years")) {
  if (!identical(res$bloom[[k]], res$leaf[[k]])) {
    stop("'", k, "' differs between the bloom and leaf results")
  }
}
for (ix in indices) {
  for (p in names(res[[ix]]$periods)) {
    stopifnot(nrow(res[[ix]]$lg) == nrow(res[[ix]]$tr_real[[p]]))
  }
}

periods     <- res$bloom$periods
years       <- res$bloom$years
laea_km     <- res$bloom$laea_km
cellsize_km <- res$bloom$cellsize_km
period_lab  <- sapply(periods, function(p) paste0(min(p), "-", max(p)))

# Align prediction grids ------------------------------------------------------#
# Each model's land grid is built around its own stations, so the grids can
# differ at the margins. Keep cells present in both.

grid_b <- res$bloom$lg %>%
  transmute(x = round(x, 6), y = round(y, 6), ib = row_number())
grid_l <- res$leaf$lg %>%
  transmute(x = round(x, 6), y = round(y, 6), il = row_number())
cells <- inner_join(grid_b, grid_l, by = c("x", "y"))
ib <- cells$ib
il <- cells$il
cells <- cells[, c("x", "y")]

cat("Bloom cells:", nrow(grid_b), "| leaf cells:", nrow(grid_l),
    "| shared:", nrow(cells), "\n")
cat("Dropped (bloom only):", nrow(grid_b) - nrow(cells),
    "| dropped (leaf only):", nrow(grid_l) - nrow(cells), "\n")

# Helpers (from visualizing-results.R) ----------------------------------------#

# Per-cell summaries of a cells x draws matrix of trend draws
#   lfsr: local false sign rate = share of draws with the less likely sign
#   sig_cell: cell's own lfsr <= q (at least 1 - q of draws agree on sign)
#   sig_simul: simultaneous band excludes 0 (all cells jointly, ~95%)
summarize_trend_draws <- function(m, q = 0.05) {
  est   <- rowMeans(m)
  sd    <- apply(m, 1, sd)
  p_pos <- rowMeans(m > 0)
  lfsr  <- pmin(p_pos, 1 - p_pos)
  zmax  <- apply(abs((m - est) / sd), 2, max)
  c95   <- unname(quantile(zmax, 1 - q))
  data.frame(
    est = est, sd = sd,
    lwr = apply(m, 1, quantile, q / 2),
    upr = apply(m, 1, quantile, 1 - q / 2),
    p_pos = p_pos, lfsr = lfsr,
    sig_cell  = lfsr <= q,
    sig_simul = (est - c95 * sd > 0) | (est + c95 * sd < 0),
    c95 = c95
  )
}

# Base map layers
world <- ne_countries(scale = "medium", returnclass = "sf")
sf_use_s2(FALSE)
north <- st_crop(world, xmin = -180, xmax = 180, ymin = 20, ymax = 90)
north <- st_transform(north, laea_km)
land  <- st_union(st_make_valid(north))
sf_use_s2(TRUE)

# Color-scale limit for a set of values: largest absolute 2nd/98th percentile
scale_lim <- function(v) max(abs(quantile(v, c(0.02, 0.98))))

# Limits for a named list of value vectors (one per period)
period_lims <- function(vals) {
  lims <- sapply(vals, scale_lim)
  if (shared_color_ramp) lims[] <- max(lims)
  lims
}

# Map of a trend surface; optional `sig` fades non-significant cells.
# RdBu with direction = 1: red = negative (gap shrinking), blue = positive
# (gap lengthening).
trend_map <- function(cells, fill, sig = NULL, lim = NULL, title = NULL,
                      legend_name = "Change in gap\n(days/decade)") {
  cells$fill <- fill
  if (is.null(lim)) lim <- scale_lim(fill)
  cells$sig <- if (is.null(sig)) TRUE else sig
  ggplot() +
    geom_sf(data = land, fill = "gray", color = "grey30", linewidth = 0.25) +
    geom_tile(data = cells, aes(x, y, fill = fill, alpha = sig),
              width = cellsize_km, height = cellsize_km) +
    geom_sf(data = north, fill = NA, color = "grey30", linewidth = 0.25) +
    scale_alpha_manual(values = c("FALSE" = 0.15, "TRUE" = 1), guide = "none") +
    scale_fill_distiller(palette = "RdBu", limits = c(-lim, lim),
                         oob = scales::squish, direction = 1,
                         name = legend_name) +
    coord_sf(crs = laea_km, xlim = range(cells$x), ylim = range(cells$y),
             expand = FALSE) +
    labs(title = title) +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white"),
          legend.position = "bottom")
}

# Raster from cell coordinates + numeric layer columns
cells_to_rast <- function(df, layers) {
  r <- rast(df[, c("x", "y", layers)], type = "xyz", crs = laea_km)
  names(r) <- layers
  r
}

# 2a) Realized local trends for each index, on shared cells -------------------#

real_idx <- lapply(indices, function(ix) {
  rows <- if (ix == "bloom") ib else il
  out <- lapply(names(periods), function(p) {
    bind_cols(cells,
              summarize_trend_draws(res[[ix]]$tr_real[[p]][rows, , drop = FALSE],
                                    q))
  })
  setNames(out, names(periods))
})
names(real_idx) <- indices

# 2b) Gap trend: bloom minus leaf, per draw -----------------------------------#

gap <- lapply(names(periods), function(p) {
  gd <- res$bloom$tr_real[[p]][ib, , drop = FALSE] -
    res$leaf$tr_real[[p]][il, , drop = FALSE]
  bind_cols(cells, summarize_trend_draws(gd, q))
})
names(gap) <- names(periods)

# Check: mean of differenced draws equals difference of means
for (p in names(periods)) {
  stopifnot(isTRUE(all.equal(gap[[p]]$est,
                             real_idx$bloom[[p]]$est - real_idx$leaf[[p]]$est)))
}

for (p in names(periods)) {
  s <- gap[[p]]
  cat("\n== Bloom-leaf gap trend,", period_lab[[p]], "==\n")
  cat("Mean gap trend (days/decade):", round(mean(s$est), 3), "\n")
  cat("Area with clear sign (lfsr <= q):", round(mean(s$sig_cell), 3), "\n")
  cat("  ...lengthening:", round(mean(s$sig_cell & s$est > 0), 3), "\n")
  cat("  ...shrinking:  ", round(mean(s$sig_cell & s$est < 0), 3), "\n")
  cat("Simultaneous band:", round(mean(s$sig_simul), 3), "\n")
}

# Maps ------------------------------------------------------------------------#

lim_gap <- period_lims(lapply(gap, function(s) s$est))

maps_gap <- lapply(names(periods), function(p) {
  trend_map(gap[[p]], gap[[p]]$est, sig = gap[[p]]$sig_cell,
            lim = lim_gap[[p]],
            title = paste("Trend in bloom-leaf gap,", period_lab[[p]],
                          "(faded: sign uncertain)"))
})
names(maps_gap) <- names(periods)
for (p in names(periods)) print(maps_gap[[p]])

maps_gap_all <- lapply(names(periods), function(p) {
  trend_map(gap[[p]], gap[[p]]$est, lim = lim_gap[[p]],
            title = paste("Trend in bloom-leaf gap,", period_lab[[p]]))
})
names(maps_gap_all) <- names(periods)
for (p in names(periods)) print(maps_gap_all[[p]])

# Hemisphere-wide gap trend ---------------------------------------------------#
# ann_draws (years x nsim) are land averages over each model's own grid, so
# the two averages cover slightly different areas if the grids differ.

ann_gap <- res$bloom$ann_draws - res$leaf$ann_draws

realized_slope <- function(per, draws) {
  i <- match(per, years)
  x <- (years[i] - mean(years[i])) / 10
  w <- (x - mean(x)) / sum((x - mean(x))^2)
  as.vector(crossprod(w, draws[i, , drop = FALSE]))
}

hemi_gap <- do.call(rbind, lapply(names(periods), function(p) {
  r <- realized_slope(periods[[p]], ann_gap)
  data.frame(period = period_lab[[p]],
             est = mean(r),
             lwr = unname(quantile(r, q / 2)),
             upr = unname(quantile(r, 1 - q / 2)),
             p_lengthening = mean(r > 0))
}))
hemi_gap   # days/decade

# Save rasters ----------------------------------------------------------------#

# gap_layers <- c("est", "sd", "lwr", "upr", "lfsr", "sig_cell")
# idx_layers <- c("est", "lfsr")
# 
# for (p in names(periods)) {
#   g <- gap[[p]]
#   g$sig_cell <- as.integer(g$sig_cell)
#   writeRaster(cells_to_rast(g, gap_layers),
#               paste0("output/bloom-leaf-gap-trend-", period_lab[[p]], ".tif"),
#               overwrite = TRUE)
#   for (ix in indices) {
#     writeRaster(cells_to_rast(real_idx[[ix]][[p]], idx_layers),
#                 paste0("output/", ix, "-trend-realized-", period_lab[[p]],
#                        ".tif"),
#                 overwrite = TRUE)
#   }
# }

# Save figures and table ------------------------------------------------------#

for (p in names(periods)) {
  # ggsave(paste0("output/bloom-leaf-gap-trendmap-", period_lab[[p]], ".png"),
  #        maps_gap[[p]], width = 6.5, height = 6.5, units = "in", dpi = 600)
  ggsave(paste0("output/bloom-leaf-gap-trendmap-", period_lab[[p]], "-all.png"),
         maps_gap_all[[p]], width = 6.5, height = 6.5, units = "in", dpi = 600)
}

# write.csv(hemi_gap, "output/bloom-leaf-gap-hemi-trends.csv", row.names = FALSE)
