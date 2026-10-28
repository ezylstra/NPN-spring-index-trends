library(dplyr)
library(stringr)
library(tidyr)
library(sf)
library(sdmTMB)
library(fmesher)
library(rnaturalearth)
library(ggplot2)
library(RColorBrewer)

# Read in results -------------------------------------------------------------#

# Set folder
folder <- "leaf-1991-2024-northamerica/"

# Read in sdmTMB model - necessary?
# fit1 <- readRDS(paste0(folder, "si-leaf-1991-2024-northamerica-fit.rds"))

# Read in saved objects (had been saved as a list)
results <- readRDS("leaf-1991-2024-northamerica/si-leaf-1991-2024-northamerica-results.rds")
list2env(results, envir = .GlobalEnv)
rm(results)

# Look at non-output objects --------------------------------------------------#

# Session info from model run
session

# Information about data
years # years included in data
yr_mean # mean year in dataset (used for centering year variable in model)
head(stations) # station IDs, with lat/lon and x/y in projected coords
n_by_year # number of stations with data each year

# Modeling choices
laea_na # projection used
n_vertices # Number of vertices in mesh

# Information about prediction grid
cellsize_km # Size of cell in prediction grid, in km
n_land # Number of cells in prediction grid (land only)

# Choices made when simulating values (local, hemisphere-wide estimates)
nsim  # Number of simulated values
batch # Number of simulated values per batch (done to reduce memory reqs)
seed # Seed used when simulating values

# Look at sanity check from model
fit_sanity
all(unlist(fit_sanity)) # TRUE if all individual checks were TRUE

# Model output ----------------------------------------------------------------#

# Estimates from model
fit_fixed
fit_ranpar

# Local trends ----------------------------------------------------------------#

head(lg) # rows match tr
# prediction cell coords + estimate of local trend (change per decade; trend_pt)
# trend_pt = par_hat[ib[2]] + as.vector(A %*% par_hat[iz])
  # or from a predict() call: fixef(fit1)["yrc"] + preds$zeta_s_yrc
  # This is the model's point estimate with no simulation noise

dim(tr) #rows match lg
# For each simulated draw, compute b + zeta(s). Results the same as predicting
# each cell at yrc = 0 and yrc = 1 in the same year and taking the difference. 
# The yearly field and the intercept field cancel, so each draw is exactly 
# b + zeta(s).

# Summarize simulated trend values
lg$trend_sim <- rowMeans(tr)
lg$sd    <- apply(tr, 1, sd)
lg$lwr   <- apply(tr, 1, quantile, 0.025)
lg$upr   <- apply(tr, 1, quantile, 0.975)
lg$p_pos <- rowMeans(tr > 0)

# Check that trend point estimate is similar to mean of simulated values:
cor(lg$trend_pt, lg$trend_sim)

# Create a map with local trends
lim1 <- max(abs(quantile(lg$trend_pt, c(0.02, 0.98))))

world <- ne_countries(scale = "medium", returnclass = "sf")
sf_use_s2(FALSE)
north <- st_crop(world, xmin = -180, xmax = -50, ymin = 20, ymax = 75)
north <- st_transform(north, laea_na)
land  <- st_union(st_make_valid(north))
sf_use_s2(TRUE)

trendmap <- ggplot() +
  geom_sf(data = land, fill = "gray", color = "grey30", linewidth = 0.25) +
  geom_tile(data = lg, aes(x, y, fill = trend_pt), width = 50, height = 50) +
  geom_sf(data = north, fill = NA, color = "grey30", linewidth = 0.25) +
  scale_fill_distiller(palette = "RdBu", limits = c(-lim1, lim1),
                       oob = scales::squish, direction = 1,
                       name = "Trend\n(per decade)") +
  coord_sf(crs = laea_na, xlim = range(lg$x), ylim = range(lg$y),
           expand = FALSE) +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white"))
trendmap

# Assess "significance": Bayesian FDR on the trend sign
  # p_pos is the proportion of draws > 0.
  # lfsr (local false sign rate) is the proportion of draws where the the sign 
  # of the trend is wrong and is min(p_pos, 1 − p_pos). [A cell with 
  # p_pos = 0.97 has lfsr = 0.03]
  # Rank the cells by lfsr from most to least certain, and keep adding cells 
  # while the average lfsr of the flagged set stays ≤ 0.05.
  # All cells have equal area (equal-area projection), so this is an area-based
  # rate.

q <- 0.05
lg$lfsr <- pmin(lg$p_pos, 1 - lg$p_pos)
o  <- order(lg$lfsr)
cm <- cumsum(lg$lfsr[o]) / seq_along(o)
k  <- max(c(0, which(cm <= q)))
lg$sig_fdr <- FALSE
lg$sig_fdr[o[seq_len(k)]] <- TRUE
mean(lg$sig_fdr)                 # proportion of land area flagged
mean(lg$lfsr[lg$sig_fdr])        # realized average lfsr of the flagged set

# Map: FDR-flagged cells
trend_fdr <- ggplot() +
  geom_sf(data = land, fill = "gray", color = "grey30", linewidth = 0.25) +
  geom_tile(data = lg, aes(x, y, fill = trend_pt, alpha = sig_fdr),
            width = 50, height = 50) +
  geom_sf(data = north, fill = NA, color = "grey30", linewidth = 0.25) +
  scale_alpha_manual(values = c("FALSE" = 0.15, "TRUE" = 1), guide = "none") +
  scale_fill_distiller(palette = "RdBu", limits = c(-lim1, lim1),
                       oob = scales::squish, direction = 1,
                       name = "Trend\n(per decade)") +
  coord_sf(crs = laea_na, xlim = range(lg$x), ylim = range(lg$y),
           expand = FALSE) +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),        
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white"))
trend_fdr

# Assess "significance": Simultaneous 95% band (stricter; all cells jointly) --#
# Max-standardized-deviation method: the band covers the whole trend surface
# with ~95% probability.
zmax <- apply(abs((tr - lg$trend_sim) / lg$sd), 2, max)
c95  <- quantile(zmax, 0.95)
c95 # compare with 1.96 for pointwise intervals
lg$sig_simul <- (lg$trend_sim - c95 * lg$sd > 0) | (lg$trend_sim + c95 * lg$sd < 0)
mean(lg$sig_simul)
# No cells flagged! (Maybe not too surprising since trends are small)

# Comparing methods:
  # They control different error rates, much like Benjamini–Hochberg versus 
  # Bonferroni:
  # Simultaneous band controls the chance of making any wrong call anywhere 
  # on the map at about 5%. It’s very strict and flags fewer cells.
  # FDR allows some wrong calls but limits them to about 5% of the flagged 
  # area. It’s less strict and flags more cells.

# Hemisphere-wide annual values and trend -------------------------------------#

# Notes on trends:
# The realized trend describes the land-averaged leaf index over study period, 
# including that period’s anomalies. The underlying trend is the 
# model’s underlying trend after the anomalies are treated as noise. For a 
# summary statement like “did leaf-out get earlier across North America from 
# 1991 to 2024?”, the realized version is probably a better choice.

# Model-estimated underlying trend
summary(underlying_trend) # per draw, days/decade
  # This is the model's underlying trend. It's calculated by summarizing 
  # the simulated values across cells in each year. Simulated values are 
  # the fixed trend effect + spatially varying yrc field (from tr object).
  # They don't take the yearly anomaly fields into account

# Annual values (from ann_sum, which comes from ann_draws)
ann_sum
  # Land-averaged annual values, one per year. Each value in ann_draws is the 
  # average over all land cells within a draw. This does take yearly anomaly
  # fields into account.

# Realized trend with proper uncertainty
  # OLS slope through the annual means, within each draw
# Mean, 95% confidence interval (computed from 1000 draws from the approximate
# joint distribution of the fixed and random effects)
c(est = mean(realized_trend), quantile(realized_trend, c(0.025, 0.975)))
# Probability the realized trend is positive
mean(realized_trend > 0)   

# Realized trend per year
slope_yr <- mean(realized_trend) / 10 # days per year

# Plot without a confidence interval around the estimated trend
hemiplot <- ggplot(ann_sum, aes(x = year, y = est)) +
  geom_linerange(aes(ymin = lwr, ymax = upr), color = "gray") +
  geom_point() +
  geom_abline(intercept = mean(ann_sum$est) - slope_yr * mean(ann_sum$year),
              slope = slope_yr, 
              linetype = 2) +
  labs(y = "Mean leaf index (DOY)", x = NULL)
hemiplot

# Calculate uncertainty around the trend
yr <- ann_sum$year
xcc <- (yr - mean(yr)) / 10   # Decades from the midpoint year

# years x nsim: each column is one draw's fitted line
line_draws <- sweep(outer(xcc, realized_trend), 2, colMeans(ann_draws), "+")

ann_sum$fit     <- rowMeans(line_draws)
ann_sum$fit_lwr <- apply(line_draws, 1, quantile, 0.025)
ann_sum$fit_upr <- apply(line_draws, 1, quantile, 0.975)

hemiplot_ci <- ggplot(ann_sum, aes(x = year)) +
  geom_ribbon(aes(ymin = fit_lwr, ymax = fit_upr), fill = "steelblue", alpha = 0.2) +
  geom_line(aes(y = fit), linetype = 2) +
  geom_linerange(aes(ymin = lwr, ymax = upr), color = "gray") +
  geom_point(aes(y = est)) +
  labs(y = "Mean leaf index (DOY)", x = NULL)
hemiplot_ci

# For reference, how the sweep function works ---------------------------------#
# m <- matrix(c(-1, 0, 1,     # draw 1: slope 1
#               -2, 0, 2),    # draw 2: slope 2
#             nrow = 3)
# sweep(m, 2, c(100, 120), "+")
