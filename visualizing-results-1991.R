# SI-x indices, Northern Hemisphere 1991-2024: VISUALIZE RESULTS
#
# Reads the results list saved by si-3-run-1991.R (on the HPC) and makes maps
# of local trends and a plot of hemisphere-wide annual values and the realized
# trend, for models fitted to 1991-2024 only.
#
# Two kinds of trend appear below:
#   UNDERLYING local trend (tr, lg$trend_pt): the model's straight-line trend
#     at each location, b + zeta(s).
#   REALIZED trend (tr_real, realized_trend): the OLS slope through the
#     modeled annual values, INCLUDING each year's anomaly field. This is the
#     trend that actually happened over 1991-2024 according to the model.
#
# Note: these are not the same as the 1991-2024 realized trends in
# visualizing-results.R. Those come from a model fitted to 1955-2024, and
# because the yearly anomaly fields are estimated with shrinkage toward 0,
# they are pulled toward the 1955-2024 trend. Here the underlying trend is
# itself estimated from 1991-2024 only. Figures and tables are saved in their
# own folder, with "1991model" in the file names, so they don't overwrite
# those from the 1955-2024 models.

library(dplyr)
library(tidyr)
library(sf)
library(rnaturalearth)
library(ggplot2)

# Settings --------------------------------------------------------------------#

index_name <- "leaf"   # leaf, bloom, frz or damage
run_name   <- paste0("si-", index_name, "-1991-2024-northhemis")
region     <- "Northern Hemisphere"
q          <- 0.05   # error rate for significance (sign) calls

out_dir    <- "output/model-1991-2024"
out_prefix <- file.path(out_dir, paste0(index_name, "-1991model-"))

# Read in results -------------------------------------------------------------#

# The fit itself is not needed here (all products are in the results file)
results <- readRDS(paste0("model-fits-results/", run_name, "-results.rds"))
list2env(results, envir = .GlobalEnv)
rm(results)

period_lab <- paste0(min(years), "-", max(years))

# Run information -------------------------------------------------------------#

session      # session info from the HPC run
years        # years in the data
yr_mean      # mean year (centering for yrc)
n_by_year    # stations with data each year
head(stations)

laea_km      # projection
n_vertices   # mesh vertices
cellsize_km  # prediction cell size (km)
n_land       # land prediction cells
nsim; batch; seed   # simulation settings

fit_sanity
all(unlist(fit_sanity))   # TRUE if every check passed

fit_fixed
fit_ranpar

# Combine all estimates
fit_fixed_df <- fit_fixed %>%
  mutate(type = "fixed", .before = term) %>%
  mutate(group_name = NA, .after = type) %>%
  data.frame()
fit_ranpar_df <- fit_ranpar %>%
  mutate(type = "random", .before = group_name) %>%
  data.frame()
fit_est <- rbind(fit_fixed_df, fit_ranpar_df)

# Helpers ---------------------------------------------------------------------#

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

# Bayesian FDR (average lfsr of flagged set <= q). Reported for reference
# only: when most cells have lfsr ~ 0, it can flag every cell, including
# cells whose own sign is close to a coin flip.
fdr_flag <- function(lfsr, q = 0.05) {
  o  <- order(lfsr)
  cm <- cumsum(lfsr[o]) / seq_along(o)
  k  <- max(c(0, which(cm <= q)))
  flag <- rep(FALSE, length(lfsr))
  flag[o[seq_len(k)]] <- TRUE
  flag
}

# Significance summary table for one set of cell summaries
sig_table <- function(s) {
  s %>%
    mutate(direction = if_else(est > 0, "later", "earlier")) %>%
    group_by(direction, sig_cell) %>%
    summarize(n_cells   = n(),
              area_prop = n() / nrow(s),
              trend_min = min(est), trend_mean = mean(est),
              trend_max = max(est),
              lfsr_min  = min(lfsr), lfsr_mean = mean(lfsr),
              lfsr_max  = max(lfsr),
              .groups = "drop") %>%
    data.frame()
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

# Map of a trend surface; optional `sig` fades non-significant cells
trend_map <- function(cells, fill, sig = NULL, lim = NULL, title = NULL) {
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
                         name = "Trend\n(days/decade)") +
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

# Underlying local trends -----------------------------------------------------#
# Each draw is b + zeta(s); the intercept and yearly fields cancel in a
# per-decade difference.

und <- summarize_trend_draws(tr, q)
lg  <- bind_cols(lg[, c("x", "y", "trend_pt")], und)

cor(lg$trend_pt, lg$est)      # point estimate vs mean of draws (should be ~1)

# Significance
mean(lg$sig_cell)                    # land area with a clear sign (lfsr <= q)
mean(lg$lfsr[lg$sig_cell])           # realized average lfsr of that area
mean(lg$sig_simul)                   # simultaneous band (strictest)
und$c95[1]                           # band multiplier (vs 1.96 pointwise)
mean(fdr_flag(lg$lfsr, q))           # Bayesian FDR, for reference only
sig_table(lg)

map_und     <- trend_map(lg, lg$trend_pt,
                         title = paste0("Underlying local trend, ", period_lab))
map_und_sig <- trend_map(lg, lg$trend_pt, sig = lg$sig_cell,
                         title = paste0("Underlying local trend, ", period_lab,
                                        " (faded: sign uncertain)"))
map_und
map_und_sig

# Realized local trends -------------------------------------------------------#

real <- bind_cols(lg[, c("x", "y")], summarize_trend_draws(tr_real, q))

# Significance
mean(real$sig_cell)                  # land area with a clear sign (lfsr <= q)
mean(real$lfsr[real$sig_cell])       # average lfsr of that area
mean(real$sig_simul)                 # simultaneous band (strictest)
sig_table(real)

# Both maps use the same color scale
lim_real <- scale_lim(real$est)

map_real     <- trend_map(real, real$est, sig = real$sig_cell, lim = lim_real,
                          title = paste0("Realized local trend, ", period_lab,
                                         " (faded: sign uncertain)"))
map_real_all <- trend_map(real, real$est, lim = lim_real,
                          title = paste0("Realized local trend, ", period_lab))
map_real
map_real_all

# Hemisphere-wide annual values and realized trend ----------------------------#
# ann_draws: years x nsim, each value the land-average of one draw's annual
# surface (includes the yearly anomaly fields).
# realized_trend: OLS slope through the land-averaged annual values within
# each draw. The underlying trend ignores the anomalies.

summary(underlying_trend)               # per draw, days/decade

hemi_tab <- data.frame(
  period  = period_lab,
  est     = mean(realized_trend),
  lwr     = unname(quantile(realized_trend, q / 2)),
  upr     = unname(quantile(realized_trend, 1 - q / 2)),
  p_later = mean(realized_trend > 0)    # probability the trend is positive
)
hemi_tab   # days/decade

# Fitted line + band, from the same draws
x_t <- (years - mean(years)) / 10
line_draws <- sweep(outer(x_t, realized_trend), 2, colMeans(ann_draws), "+")
line_df <- data.frame(period = period_lab, year = years,
                      fit = rowMeans(line_draws),
                      fit_lwr = apply(line_draws, 1, quantile, q / 2),
                      fit_upr = apply(line_draws, 1, quantile, 1 - q / 2))

trend_label <- with(hemi_tab, setNames(
  sprintf("%s: %.2f days/decade (%.2f to %.2f)", period, est, lwr, upr),
  period))

# Same color as the 1991-2024 trend in the plots from visualizing-results.R
hemi_color <- setNames(RColorBrewer::brewer.pal(3, "Dark2")[2], period_lab)

hemiplot <- ggplot() +
  geom_ribbon(data = line_df,
              aes(x = year, ymin = fit_lwr, ymax = fit_upr, fill = period),
              alpha = 0.2) +
  geom_line(data = line_df, aes(x = year, y = fit, color = period),
            linetype = 2) +
  geom_linerange(data = ann_sum, aes(x = year, ymin = lwr, ymax = upr),
                 color = "gray") +
  geom_point(data = ann_sum, aes(x = year, y = est), size = 0.8) +
  scale_color_manual(values = hemi_color, labels = trend_label,
                     name = "Realized trend") +
  scale_fill_manual(values = hemi_color, labels = trend_label,
                    name = "Realized trend") +
  labs(y = paste("Mean", index_name, "index (DOY)"), x = NULL,
       title = paste("Land-averaged", index_name, "index,", region)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "vertical")
hemiplot

# Save figures ----------------------------------------------------------------#

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(paste0(out_prefix, "trendmap-underlying.png"), map_und,
       width = 6.5, height = 6.5, units = "in", dpi = 600)
ggsave(paste0(out_prefix, "trendmap-underlying-sig.png"), map_und_sig,
       width = 6.5, height = 6.5, units = "in", dpi = 600)
ggsave(paste0(out_prefix, "trendmap-realized-", period_lab, ".png"), map_real,
       width = 6.5, height = 6.5, units = "in", dpi = 600)
ggsave(paste0(out_prefix, "trendmap-realized-", period_lab, "-all.png"),
       map_real_all, width = 6.5, height = 6.5, units = "in", dpi = 600)
ggsave(paste0(out_prefix, "hemi-trend-", period_lab, ".png"), hemiplot,
       width = 6.5, height = 5, units = "in", dpi = 600)

# Save estimates --------------------------------------------------------------#

write.csv(fit_est,
          file = paste0(out_prefix, "model-ests-", period_lab, ".csv"),
          row.names = FALSE)
