# We'll only use a small number of knots so this example runs quickly
# but you will likely want to use many more in applied situations.

library(sdmTMB)
library(ggplot2)
d <- pcod
pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30) # a coarse mesh for example speed
m <- sdmTMB(
  data = d, formula = density ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
  silent = F, 
  time = "year", spde = pcod_spde, family = tweedie(link = "log")
)

# Predictions at original data locations -------------------------------

predictions <- predict(m)
head(predictions)

predictions$resids <- residuals(m) # randomized quantile residuals

ggplot(predictions, aes(X, Y, col = resids)) + scale_colour_gradient2() +
  geom_point() + facet_wrap(~year)

hist(predictions$resids)

qqnorm(predictions$resids);abline(a = 0, b = 1)

# Predictions onto new data --------------------------------------------

predictions <- predict(m, newdata = qcs_grid)

# A short function for plotting our predictions:
plot_map <- function(dat, column = "est") {
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_raster() +
    facet_wrap(~year) +
    coord_fixed()
}

plot_map(predictions, "exp(est)") +
  scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Prediction (fixed effects + all random effects)")

plot_map(predictions, "exp(est_non_rf)") +
  ggtitle("Prediction (fixed effects and any time-varying effects)") +
  scale_fill_viridis_c(trans = "sqrt")

plot_map(predictions, "est_rf") +
  ggtitle("All random field estimates") +
  scale_fill_gradient2()

plot_map(predictions, "omega_s") +
  ggtitle("Spatial random effects only") +
  scale_fill_gradient2()

plot_map(predictions, "epsilon_st") +
  ggtitle("Spatiotemporal random effects only") +
  scale_fill_gradient2()

# Visualizing a marginal effect ----------------------------------------

nd <- data.frame(depth_scaled =
                   seq(min(d$depth_scaled), max(d$depth_scaled), length.out = 100))
nd$depth_scaled2 <- nd$depth_scaled^2

# You'll need at least one time element. If time isn't also a fixed effect
# then it doesn't matter what you pick:
nd$year <- 2003L
p <- predict(m, newdata = nd, se_fit = TRUE, re_form = NA)
ggplot(p, aes(depth_scaled, exp(est),
              ymin = exp(est - 1.96 * est_se), ymax = exp(est + 1.96 * est_se))) +
  geom_line() + geom_ribbon(alpha = 0.4)

# Plotting marginal effect of a spline ---------------------------------

m_gam <- sdmTMB(
  data = d, formula = density ~ 0 + as.factor(year) + s(depth_scaled, k = 3),
  time = "year", spde = pcod_spde, family = tweedie(link = "log")
)
nd <- data.frame(depth_scaled =
                   seq(min(d$depth_scaled), max(d$depth_scaled), length.out = 100))
nd$year <- 2003L
p <- predict(m_gam, newdata = nd, se_fit = TRUE, re_form = NA)
ggplot(p, aes(depth_scaled, exp(est),
              ymin = exp(est - 1.96 * est_se), ymax = exp(est + 1.96 * est_se))) +
  geom_line() + geom_ribbon(alpha = 0.4)

# Forecasting ----------------------------------------------------------
pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 15)

unique(d$year)
#> [1] 2003 2004 2005 2007 2009 2011 2013 2015 2017
m <- sdmTMB(
  data = d, formula = density ~ 1,
  ar1_fields = TRUE, # using an AR1 to have something to forecast with
  extra_time = 2019L,
  include_spatial = FALSE,
  time = "year", spde = pcod_spde, family = tweedie(link = "log")
)

# Add a year to our grid:
grid2019 <- qcs_grid[qcs_grid$year == max(qcs_grid$year), ]
grid2019$year <- 2019L # `L` because `year` is an integer in the data
qcsgrid_forecast <- rbind(qcs_grid, grid2019)

predictions <- predict(m, newdata = qcsgrid_forecast)
plot_map(predictions, "exp(est)") +
  scale_fill_viridis_c(trans = "log10")
plot_map(predictions, "epsilon_st") +
  scale_fill_gradient2()

# Estimating local trends ----------------------------------------------

pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 25)
m <- sdmTMB(data = pcod, formula = density ~ depth_scaled + depth_scaled2,
            spde = pcod_spde, family = tweedie(link = "log"),
            spatial_trend = TRUE, time = "year", spatial_only = TRUE)
p <- predict(m, newdata = qcs_grid)

plot_map(p, "zeta_s") +
  ggtitle("Spatial slopes") +
  scale_fill_gradient2()

plot_map(p, "est_rf") +
  ggtitle("Random field estimates") +
  scale_fill_gradient2()

plot_map(p, "exp(est_non_rf)") +
  ggtitle("Prediction (fixed effects only)") +
  scale_fill_viridis_c(trans = "sqrt")

plot_map(p, "exp(est)") +
  ggtitle("Prediction (fixed effects + all random effects)") +
  scale_fill_viridis_c(trans = "sqrt")
# }
