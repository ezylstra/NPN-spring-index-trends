# SI-x indices, Northern Hemisphere 1951-2024: VISUALIZE RESULTS
#
# Reads the results list saved by si-3-run.R (on the HPC) and
# makes maps of local trends and plots of hemisphere-wide annual values and
# realized trends, for the full period and for any sub-periods saved in
# `periods` (default: full = 1951-2024, recent = 1991-2024).
#
# Two kinds of trend appear below:
#   UNDERLYING local trend (tr, lg$trend_pt): the model's straight-line trend
#     at each location, b + zeta(s). The model fits ONE line per location for
#     1951-2024, so this map is the same for any sub-period.
#   REALIZED trend (tr_real, realized_by_period): the OLS slope through the
#     modeled annual values over a period, INCLUDING each year's anomaly
#     field. This is the trend that actually happened over that period
#     according to the model, and it can differ between periods.
#
# Note: the yearly anomaly fields are estimated with shrinkage toward 0
# (more so where stations are sparse), so realized local trends for short
# periods are pulled toward the long-term trend. A model fitted only to
# 1991-2024 lets that period's trend vary freely; the two approaches answer
# slightly different questions.

library(dplyr)
library(tidyr)
library(sf)
library(rnaturalearth)
library(ggplot2)

# Settings --------------------------------------------------------------------#

index_name <- "bloom"   # leaf, bloom, frz or damage
run_name   <- paste0("si-", index_name, "-1951-2024-northhemis")
folder     <- paste0(index_name, "-1951-2024-northhemis/")   # local folder
region     <- "Northern Hemisphere"
q          <- 0.05   # error rate for significance (sign) calls

# Read in results -------------------------------------------------------------#

# The fit itself is not needed here (all products are in the results file)
results <- readRDS(paste0(folder, run_name, "-results.rds"))
list2env(results, envir = .GlobalEnv)
rm(results)

has_real <- exists("tr_real")   # FALSE for results from the older run script
if (!has_real) {
  message("Results have no realized local trends (tr_real); re-run the ",
          "updated run script on the HPC to get the period maps. ",
          "Hemisphere-wide period trends are computed from ann_draws below.")
  periods <- list(full = years, recent = 1991:2024)
}
period_lab <- sapply(periods, function(p) paste0(min(p), "-", max(p)))

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

# Map of a trend surface; optional `sig` fades non-significant cells
trend_map <- function(cells, fill, sig = NULL, lim = NULL, title = NULL) {
  cells$fill <- fill
  if (is.null(lim)) lim <- max(abs(quantile(fill, c(0.02, 0.98))))
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

# Underlying local trends (whole model, 1951-2024) ----------------------------#
# Each draw is b + zeta(s); the intercept and yearly fields cancel in a
# per-decade difference. Same map for any sub-period.

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
                         title = paste("Underlying local trend,",
                                       min(years), "-", max(years)))
map_und_sig <- trend_map(lg, lg$trend_pt, sig = lg$sig_cell,
                         title = paste("Underlying local trend,",
                                       min(years), "-", max(years),
                                       "(faded: sign uncertain)"))
map_und
map_und_sig

# Realized local trends by period ---------------------------------------------#

if (has_real) {
  real <- lapply(names(periods), function(p) {
    s <- summarize_trend_draws(tr_real[[p]], q)
    bind_cols(lg[, c("x", "y")], s)
  })
  names(real) <- names(periods)

  for (p in names(periods)) {
    cat("\n== Realized local trend,", period_lab[[p]], "==\n")
    s <- real[[p]]
    cat("Area with clear sign (lfsr <= q):", round(mean(s$sig_cell), 3), "\n")
    cat("Average lfsr of that area:", round(mean(s$lfsr[s$sig_cell]), 4), "\n")
    cat("Simultaneous band:", round(mean(s$sig_simul), 3), "\n")
    print(sig_table(s))
  }

  # Common color scale across periods so maps can be compared directly
  lim_real <- max(sapply(real, function(s)
    max(abs(quantile(s$est, c(0.02, 0.98))))))

  maps_real <- lapply(names(periods), function(p) {
    trend_map(real[[p]], real[[p]]$est, sig = real[[p]]$sig_cell,
              lim = lim_real,
              title = paste("Realized local trend,", period_lab[[p]],
                            "(faded: sign uncertain)"))
  })
  names(maps_real) <- names(periods)
  print(maps_real$full)     # print() needed: plots inside if {} don't auto-show
  print(maps_real$recent)

  # How different is the recent period from the long-term trend?
  diff_recent <- tr_real$recent - tr       # per draw, days/decade
  lg$diff_recent     <- rowMeans(diff_recent)
  lg$diff_recent_sig <- pmin(rowMeans(diff_recent > 0),
                             rowMeans(diff_recent < 0)) <= q
  map_diff <- trend_map(lg, lg$diff_recent, sig = lg$diff_recent_sig,
                        title = paste("Realized", period_lab[["recent"]],
                                      "minus underlying trend"))
  print(map_diff)
}

# Hemisphere-wide annual values and realized trends ---------------------------#
# ann_draws: years x nsim, each value the land-average of one draw's annual
# surface (includes the yearly anomaly fields).
# Realized trend: OLS slope through the land-averaged annual values within
# each draw, over each period. The underlying trend ignores the anomalies.

summary(underlying_trend)               # per draw, days/decade

realized_slope <- function(per) {
  i <- match(per, years)
  x <- (years[i] - mean(years[i])) / 10
  w <- (x - mean(x)) / sum((x - mean(x))^2)
  as.vector(crossprod(w, ann_draws[i, , drop = FALSE]))
}
real_hemi <- lapply(periods, realized_slope)

hemi_tab <- do.call(rbind, lapply(names(periods), function(p) {
  r <- real_hemi[[p]]
  data.frame(period = period_lab[[p]],
             est = mean(r),
             lwr = unname(quantile(r, q / 2)),
             upr = unname(quantile(r, 1 - q / 2)),
             p_later = mean(r > 0))      # probability the trend is positive
}))
hemi_tab   # days/decade

# Fitted line + band for each period, from the same draws
line_df <- do.call(rbind, lapply(names(periods), function(p) {
  i  <- match(periods[[p]], years)
  x  <- (years[i] - mean(years[i])) / 10
  ld <- sweep(outer(x, real_hemi[[p]]), 2,
              colMeans(ann_draws[i, , drop = FALSE]), "+")
  data.frame(period = period_lab[[p]], year = years[i],
             fit = rowMeans(ld),
             fit_lwr = apply(ld, 1, quantile, q / 2),
             fit_upr = apply(ld, 1, quantile, 1 - q / 2))
}))

trend_labels <- with(hemi_tab, setNames(
  sprintf("%s: %.2f days/decade (%.2f to %.2f)", period, est, lwr, upr),
  period))

hemiplot <- ggplot() +
  geom_ribbon(data = line_df,
              aes(x = year, ymin = fit_lwr, ymax = fit_upr, fill = period),
              alpha = 0.2) +
  geom_line(data = line_df, aes(x = year, y = fit, color = period),
            linetype = 2) +
  geom_linerange(data = ann_sum, aes(x = year, ymin = lwr, ymax = upr),
                 color = "gray") +
  geom_point(data = ann_sum, aes(x = year, y = est)) +
  scale_color_brewer(palette = "Dark2", labels = trend_labels,
                     name = "Realized trend") +
  scale_fill_brewer(palette = "Dark2", labels = trend_labels,
                    name = "Realized trend") +
  labs(y = paste("Mean", index_name, "index (DOY)"), x = NULL,
       title = paste("Land-averaged", index_name, "index,", region)) +
  theme(legend.position = "bottom", legend.direction = "vertical")
hemiplot

# Save figures ----------------------------------------------------------------#

ggsave(paste0(folder, "trendmap-underlying.png"), map_und,
       width = 6.5, height = 6.5, units = "in", dpi = 600)
ggsave(paste0(folder, "trendmap-underlying-sig.png"), map_und_sig,
       width = 6.5, height = 6.5, units = "in", dpi = 600)
if (has_real) {
  for (p in names(periods)) {
    ggsave(paste0(folder, "trendmap-realized-", period_lab[[p]], ".png"),
           maps_real[[p]], width = 6.5, height = 6.5, units = "in", dpi = 600)
  }
  ggsave(paste0(folder, "trendmap-realized-recent-minus-underlying.png"),
         map_diff, width = 6.5, height = 6.5, units = "in", dpi = 600)
}
ggsave(paste0(folder, "hemi-trends.png"), hemiplot,
       width = 6.5, height = 5, units = "in", dpi = 600)
