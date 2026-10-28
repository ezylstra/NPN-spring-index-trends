# Continues from spring-index-station-data-exploration.R
# Assumes these objects exist: dat, years, mesh, fit1_update, pred_grid, north,
# lim. The centering of yrc must match the centering used to fit the model.

yr_mean <- mean(dat$year)          # yrc = (year - yr_mean) / 10
laea_km <- "+proj=laea +lat_0=90 +lon_0=0 +units=km"

# 0. Diagnostic: are residuals correlated across stations within a year? ------#
# The AR1 check looked at correlation within a station across years. Leaf-out is
# driven by weather that is shared regionally, such as a warm spring across all
# of Europe. If that is happening, residuals are correlated between stations in
# the same year. fit1 treats those residuals as independent, which makes trend
# SEs too small.

p_fit <- predict(fit1_update)                      # includes station REs
dat$res1 <- dat$leaf - p_fit$est

# (a) Year-mean residuals. With independent errors, their SD across years should
#     be close to sd(res) / sqrt(stations per year).
yr_res <- dat %>% 
  group_by(year) %>% 
  summarize(m = mean(res1), n = n())
sd(yr_res$m)
sd(dat$res1) / sqrt(mean(yr_res$n))

# (b) Correlation of residuals between station pairs, binned by distance
rw <- dat %>%
  select(station, year, res1) %>%
  pivot_wider(names_from = station, values_from = res1) %>%
  arrange(year)
C <- cor(as.matrix(rw[, -1]), use = "pairwise.complete.obs")
st_xy <- dat %>% distinct(station, x, y)
st_xy <- st_xy[match(colnames(rw)[-1], as.character(st_xy$station)), ]
D <- as.matrix(dist(st_xy[, c("x", "y")]))
ut <- upper.tri(D)
data.frame(d = D[ut], r = C[ut]) %>%
  mutate(bin = cut(d, c(0, 100, 250, 500, 1000, 2000, 4000, Inf))) %>%
  group_by(bin) %>%
  summarize(mean_r = mean(r, na.rm = TRUE), n_pairs = n())
rm(C, D)
# Large positive mean_r at short distances means yearly spatial anomaly fields
# should be added to the model (fit3 below).

# 1. Fit model with yearly spatial anomaly fields -----------------------------#
# time = "year" + spatiotemporal = "iid" adds a zero-mean GRF for each year.
# The trend structure stays the same: spatial intercept field, fixed slope, SVC
# trend field, and station intercepts.
# start <- Sys.time()
# fit3 <- sdmTMB(
#   leaf ~ yrc + (1 | station),
#   data = dat,
#   mesh = mesh,
#   time = "year",
#   spatial = "on",
#   spatiotemporal = "iid",
#   spatial_varying = ~ 0 + yrc,
#   family = gaussian()
# )
# end <- Sys.time()
# end - start

fit3 <- readRDS("model-leaf-1991-2024.rds")

sanity(fit3) # no issues!

# Random effects
tidy(fit3, effects = "ran_par", conf.int = TRUE)  # sigma_E = SD of yearly fields
# Model comparison (vast improvement over simpler model)
AIC(fit1_update, fit3)
# Compare trend estimates from the 2 models
tidy(fit1_update, conf.int = TRUE); tidy(fit3, conf.int = TRUE)  # compare yrc SE

##### DO THESE CHECKS LATER (first requires running another, more complex model)

  # 1. Does a common year effect have non-zero variance?
  dat$year_f <- factor(dat$year)
  fit4 <- update(fit3, leaf ~ yrc + (1 | station) + (1 | year_f))
  tidy(fit4, "ran_par", conf.int = TRUE)   # SD of year_f near 0 -> fit3 is fine
  tidy(fit4, conf.int = TRUE)              # compare yrc SE with fit3
  
  # 2. Naive comparison: OLS slope through the 34 annual hemisphere means
  summary(lm(est ~ I((year - mean(year)) / 10), data = hemi))

fit_use <- fit3    # or fit1_update, if the diagnostics show nothing

# 2. Land-only prediction grid -----------------------------------------------#
# Keep the 300-km buffer around stations and intersect it with land. This
# restricts the grid to land that has nearby stations, which avoids both open
# ocean and the isolated islands in Russia.
sf_use_s2(FALSE)
land <- st_union(st_make_valid(north))             # north is already in laea_km
sf_use_s2(TRUE)
pg_sf <- st_as_sf(pred_grid, coords = c("x", "y"), crs = laea_km, remove = FALSE)
pred_grid$land <- lengths(st_intersects(pg_sf, land)) > 0
# Natural Earth country polygons include some large lakes, such as the Great
# Lakes. To drop them, subtract rnaturalearth lake polygons as well.

lg <- pred_grid[pred_grid$land, c("x", "y", "station")]
n_land <- nrow(lg)
n_land

# 3. Local trend draws (land cells) ------------------------------------------#
# # Predict each cell at yrc = 0 and yrc = 1 in the same year, using the same
# # draws, and take the difference. The yearly field and the intercept field
# # cancel, so each draw is exactly b + zeta(s).
# # Use more draws than 500. The lfsr/FDR step depends on the tails.
# # Memory use is about 2 * n_land * nsim doubles. Coarsen the grid if needed.
# nsim <- 1000 # 2000 is probably better ########
# nd <- rbind(transform(lg, yrc = 0, year = years[1]),
#             transform(lg, yrc = 1, year = years[1]))
# set.seed(1)
# sims <- predict(fit_use, newdata = nd, nsim = nsim, re_form_iid = NA)
# tr <- sims[(n_land + 1):(2 * n_land), ] - sims[1:n_land, ]   # n_land x nsim
# rm(sims)
# 
# lg$trend <- rowMeans(tr)
# lg$sd    <- apply(tr, 1, sd)
# lg$lwr   <- apply(tr, 1, quantile, 0.025)
# lg$upr   <- apply(tr, 1, quantile, 0.975)
# lg$p_pos <- rowMeans(tr > 0)
# 
# # 3a. Bayesian FDR on the sign (local false sign rate) ----------------------
# # lfsr = posterior probability that the sign of the trend is wrong.
# # Rank cells by lfsr, then take the largest set whose average lfsr is <= q.
# # The expected proportion of flagged area with the wrong sign is then <= q.
# # All cells have equal area (equal-area projection), so this is an area-based
# # rate.
# q <- 0.05
# lg$lfsr <- pmin(lg$p_pos, 1 - lg$p_pos)
# o  <- order(lg$lfsr)
# cm <- cumsum(lg$lfsr[o]) / seq_along(o)
# k  <- max(c(0, which(cm <= q)))
# lg$sig_fdr <- FALSE
# lg$sig_fdr[o[seq_len(k)]] <- TRUE
# mean(lg$sig_fdr)                 # proportion of land area flagged
# mean(lg$lfsr[lg$sig_fdr])        # realized average lfsr of the flagged set
# 
# # 3b. Simultaneous 95% band (stricter; all cells jointly) --------------------
# # Max-standardized-deviation method: the band covers the whole trend surface
# # with ~95% probability.
# zmax <- apply(abs((tr - lg$trend) / lg$sd), 2, max)
# c95  <- quantile(zmax, 0.95)
# c95                                  # compare with 1.96 for pointwise intervals
# lg$sig_simul <- (lg$trend - c95 * lg$sd > 0) | (lg$trend + c95 * lg$sd < 0)
# mean(lg$sig_simul)
# 
# # Map: FDR-flagged cells
# ggplot() +
#   geom_tile(data = lg, aes(x, y, fill = trend, alpha = sig_fdr),
#             width = 50, height = 50) +
#   geom_sf(data = north, fill = NA, color = "grey30", linewidth = 0.25) +
#   scale_alpha_manual(values = c("FALSE" = 0.15, "TRUE" = 1), guide = "none") +
#   scale_fill_distiller(palette = "RdBu", limits = c(-lim, lim),
#                        oob = scales::squish, direction = 1,
#                        name = "Trend\n(per decade)") +
#   coord_sf(crs = laea_km, xlim = range(lg$x), ylim = range(lg$y),
#            expand = FALSE) +
#   theme(axis.title = element_blank(), axis.text = element_blank(),
#         axis.ticks = element_blank())

# Avoiding draws and significance testing for now (takes a very long time)

# One row per land cell. year must be one of the fitted years (fit3 has
# time = "year"), but which year you pick doesn't affect zeta_s_yrc.
lg_pred <- transform(lg[, c("x", "y", "station")], year = years[1], yrc = 0)

p3 <- predict(fit3, newdata = lg_pred, re_form_iid = NA)
p3$local_trend <- fixef(fit3)["yrc"] + p3$zeta_s_yrc   # change per decade

lim3 <- max(abs(quantile(p3$local_trend, c(0.02, 0.98))))

ggplot() +
  geom_tile(data = p3, aes(x, y, fill = local_trend), width = 50, height = 50) +
  geom_sf(data = north, fill = NA, color = "grey30", linewidth = 0.25) +
  scale_fill_distiller(palette = "RdBu", limits = c(-lim3, lim3),
                       oob = scales::squish, direction = 1,
                       name = "Trend\n(per decade)") +
  coord_sf(crs = laea_km, xlim = range(p3$x), ylim = range(p3$y),
           expand = FALSE) +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank())


# 4. Hemisphere-wide trend ----------------------------------------------------#
# Area-weighted mean of the local trend over land cells, computed in each draw.
# This is not the same as fixef(yrc). The SVC field is centered over the whole
# mesh, not over land, and the station network is dominated by the US.
hemi_trend <- colMeans(tr)
c(est = mean(hemi_trend), quantile(hemi_trend, c(0.025, 0.975)))
# Units: index units (days) per decade

# 5. Hemisphere-wide annual values -------------------------------------------#
# Predict every land cell in every year, then take the area-weighted mean.
# area = 1/n_land turns get_index()'s area-weighted sum into a mean. The SEs
# come from TMB's generalized delta method. get_index() builds its CIs from
# log_est and se on the log scale. That is fine for positive day-of-year
# values. For symmetric intervals, use est +/- 1.96 * se * est (delta method).
# If get_index() is slow, coarsen the grid (e.g., 100 km). Equal-area cells mean
# coarsening does not bias the mean.

# grid_yrs <- tidyr::expand_grid(lg[, c("x", "y", "station")], year = years) %>%
#   mutate(yrc = (year - yr_mean) / 10)
# p_idx <- predict(fit_use, newdata = grid_yrs, return_tmb_object = TRUE,
#                  re_form_iid = NA)
# hemi <- get_index(p_idx, area = rep(1 / n_land, nrow(grid_yrs)),
#                   bias_correct = FALSE)
# hemi

# Coarser land grid for the index only
st_crs(domain) <- laea_km
g200 <- st_make_grid(domain, cellsize = 200, what = "centers")
g200 <- g200[lengths(st_intersects(g200, domain)) > 0]
lg200 <- data.frame(st_coordinates(g200)) %>% rename(x = X, y = Y)
lg200_sf <- st_as_sf(lg200, coords = c("x", "y"), crs = laea_km, remove = FALSE)
lg200 <- lg200[lengths(st_intersects(lg200_sf, land)) > 0, ]
lg200$station <- factor(levels(dat$station)[1], levels = levels(dat$station))
n200 <- nrow(lg200)

# Time a 2-year run before committing to all 34 years
test <- tidyr::expand_grid(lg200, year = years[1:2]) %>%
  mutate(yrc = (year - yr_mean) / 10)
system.time({
  pt <- predict(fit3, newdata = test, return_tmb_object = TRUE, re_form_iid = NA)
  it <- get_index(pt, area = rep(1 / n200, nrow(test)), bias_correct = FALSE)
}) # Took 79 minutes for 2 years








# Trend line implied by the model, from the same draws as hemi_trend.
# This needs the hemispheric intercept at yrc = 0. That part does not cancel,
# so use fit1-style draws, or take it from the fitted index:
#   fitted line = mean(hemi$est) + (year - mean(years))/10 * mean(hemi_trend)
# (approximate when stations have unequal numbers of years)

ggplot(hemi, aes(year, est)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey80") +
  geom_line() + geom_point() +
  labs(x = NULL, y = "Hemisphere mean leaf index (DOY)") +
  theme_bw()
