library(sdmTMB)
library(dplyr)
library(ggplot2)

rm(list = ls())

load("data/ALL_REA_FISH_RAW.rdata")

df = df %>% 
  subset(REGION == "MHI") %>%
  # subset(OBS_YEAR >= 2015)
  subset(ISLAND == "Oahu")

sp = df %>% 
  group_by(TAXONNAME) %>% 
  summarise(n = sum(BIOMASS_G_M2, na.rm = T)) %>% 
  mutate(freq = n/sum(n)) %>% 
  mutate(cumsum = cumsum(freq)) %>% 
  arrange(desc(freq))

sp

df$BIOMASS_G_M2 = ifelse(df$TAXONNAME == "Aprion virescens", df$BIOMASS_G_M2, 0)

df = df %>% 
  group_by(LONGITUDE, LATITUDE, OBS_YEAR, DEPTH) %>% 
  summarise(density = mean(BIOMASS_G_M2, na.rm = T))

plot(df$density, df$DEPTH)
hist(df$density)
summary(df$density)

zone <- (floor((df$LONGITUDE[1] + 180)/6) %% 60) + 1

xy_utm = as.data.frame(cbind(utm = project(as.matrix(df[, c("LONGITUDE", "LATITUDE")]),
                                           paste0("+proj=utm +zone=", zone))))

colnames(xy_utm) = c("X", "Y"); plot(xy_utm, pch = ".", bty = 'n')

df = cbind(df, xy_utm)

rea_spde <- make_mesh(df, c("X", "Y"), cutoff = 100) # a coarse mesh for speed

plot(rea_spde, pch = "."); axis(1); axis(2)

df$year = df$OBS_YEAR
df$depth = df$DEPTH
df$depth_scaled = as.numeric(scale(log(df$depth+1)))

# alt model 1, includes a single intercept spatial random field and another random field for spatially varying slopes the represent trends over time in space. Just estimates and intercept and accounts for all other variation through the random effects
m1 <- sdmTMB(
  data = df, 
  formula = density ~ 1,
  silent = F, 
  spatial_trend = T, 
  spatial_only = T, 
  time = "year", 
  spde = rea_spde, 
  family = tweedie(link = "log")
)

# alt model 2, add independent spatiotemporal random fields for each year
m2 <- sdmTMB(
  data = df, 
  formula = density ~ 1,
  silent = F, 
  spatial_trend = T, 
  spatial_only = F, 
  time = "year", 
  spde = rea_spde, 
  family = tweedie(link = "log")
)

# alt model 3, make the spatiotemporal random fields follow an AR1 process
m3 <- sdmTMB(
  data = df,
  formula = density ~ 1,
  silent = F, 
  spatial_trend = T, 
  spatial_only = F, 
  ar1_fields = T,
  time = "year", 
  spde = rea_spde, 
  family = tweedie(link = "log")
)

# look at gradients
max(m1$gradients)
max(m2$gradients)
max(m3$gradients)

df$residuals1 <- residuals(m1)
df$residuals2 <- residuals(m2)
df$residuals3 <- residuals(m3)

qqnorm(df$residuals1, ylim = c(-5, 5), xlim = c(-5, 5));abline(a = 0, b = 1)
qqnorm(df$residuals2, ylim = c(-5, 5), xlim = c(-5, 5));abline(a = 0, b = 1)
qqnorm(df$residuals3, ylim = c(-5, 5), xlim = c(-5, 5));abline(a = 0, b = 1)

plot_map_point <- function(dat, column = "est") {
  ggplot(dat, aes_string("X", "Y", colour = column)) +
    geom_point(alpha = 0.5, size = 5) +
    facet_wrap(~year) +
    coord_fixed() + 
    ggdark::dark_theme_void()
}

plot_map_point(df, "residuals1") + scale_color_gradient2()
plot_map_point(df, "residuals2") + scale_color_gradient2()
plot_map_point(df, "residuals3") + scale_color_gradient2()

#  extract some parameter estimates
sd1 <- as.data.frame(summary(TMB::sdreport(m1$tmb_obj)))
sd2 <- as.data.frame(summary(TMB::sdreport(m2$tmb_obj)))
sd3 <- as.data.frame(summary(TMB::sdreport(m3$tmb_obj)))

r1 <- m1$tmb_obj$report()
r2 <- m2$tmb_obj$report()
r3 <- m3$tmb_obj$report()

#the estimate and 95% confidence interval on the AR1 correlation parameter:
sd3$Estimate[row.names(sd3) == "ar1_phi"]
sd3$Estimate[row.names(sd3) == "ar1_phi"] + c(-2, 2) * sd3$`Std. Error`[row.names(sd3) == "ar1_phi"]

# prediction onto new data grid
load("data/Topography_NOAA_CRM_vol10.RData")

grid = topo

# grid$longitude = round(grid$x, digits = 2)
# grid$latitude = round(grid$y, digits = 2)

grid$longitude = round(grid$x, digits = 4)
grid$latitude = round(grid$y, digits = 4)

grid = grid %>% 
  group_by(longitude, latitude) %>% 
  subset(longitude > range(df$LONGITUDE)[1]) %>% 
  subset(longitude < range(df$LONGITUDE)[2]) %>%  
  subset(latitude > range(df$LATITUDE)[1]) %>%
  subset(latitude < range(df$LATITUDE)[2]) %>% 
  summarise(depth = mean(Topography, na.rm = T)*-1) 

zone <- (floor((grid$longitude[1] + 180)/6) %% 60) + 1
xy_utm = as.data.frame(cbind(utm = project(as.matrix(grid[, c("longitude", "latitude")]),
                                           paste0("+proj=utm +zone=", zone))))

colnames(xy_utm) = c("X", "Y"); plot(xy_utm, pch = ".")

grid = cbind(grid, xy_utm)

grid_year = NULL

year = as.vector(unique(df$year))

for (y in 1:length(year)) {
  
  # y = 1
  
  grid_y = grid  
  grid_y$year = year[[y]]  
  
  grid_year = rbind(grid_year, grid_y)
  
}

grid_year$depth_scaled = as.numeric(scale(log(grid_year$depth+1)))

p1 <- predict(m1, newdata = grid_year)
p2 <- predict(m2, newdata = grid_year)
p3 <- predict(m3, newdata = grid_year)

plot_map_raster <- function(dat, column = "est") {
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_tile(aes(height = 500, width = 500)) +
    facet_wrap(~year) +
    coord_fixed() +
    scale_fill_viridis_c() + 
    ggdark::dark_theme_minimal()
}

# pick out a single year to plot since they should all be the same for the slopes. Note that these are in log space.
plot_map_raster(filter(p1, year == 2015), "zeta_s")
plot_map_raster(filter(p2, year == 2015), "zeta_s")
plot_map_raster(filter(p3, year == 2015), "zeta_s")

#predictions including all fixed and random effects plotted in log space.
plot_map_raster(p1, "est")
plot_map_raster(p2, "est")
plot_map_raster(p3, "est")

# look at just the spatiotemporal random effects for models 2 and 3:
plot_map_raster(p2, "est_rf") + scale_fill_gradient2()
plot_map_raster(p3, "est_rf") + scale_fill_gradient2()

# single spatial random effects for all three models
plot_map_raster(filter(p1, year == 2015), "omega_s")
plot_map_raster(filter(p2, year == 2015), "omega_s")
plot_map_raster(filter(p3, year == 2015), "omega_s")

p1 <- predict(m1, 
              newdata = grid_year, 
              return_tmb_object = TRUE)
p2 <- predict(m2, 
              newdata = grid_year, 
              return_tmb_object = TRUE)
p3 <- predict(m3, 
              newdata = grid_year, 
              return_tmb_object = TRUE)

index1 <- get_index(p1, bias_correct = F)
index2 <- get_index(p2, bias_correct = F)
index3 <- get_index(p3, bias_correct = F)

index1$model = "m1"
index2$model = "m2"
index3$model = "m3"

index = rbind(index1, index2, index3)

ggplot(index, aes(year, est, color = model, fill = model)) + 
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  xlab('Year') + 
  ylab('Biomass estimate (metric tonnes)') + 
  facet_wrap(~model, scales = "free_y") + 
  ggdark::dark_theme_minimal()
  

# default model 
m0 <- sdmTMB(
  data = df, 
  formula = density ~ 0 + as.factor(year) + depth + depth_scaled,
  silent = F, 
  time = "year", 
  spde = rea_spde,
  family = tweedie(link = "log")
)

predictions <- predict(m, newdata = grid_year, return_tmb_object = TRUE, area = 4)

# predictions <- predict(m)
# head(predictions)
# 
# predictions$resids <- residuals(m) # randomized quantile residuals
# 
# ggplot(predictions, aes(X, Y, col = resids)) + 
#   scale_colour_gradient2() +
#   geom_point() + 
#   facet_wrap(~year)
# 
# hist(predictions$resids)
# 
# qqnorm(predictions$resids);abline(a = 0, b = 1)

# m1 <- run_extra_optimization(m, nlminb_loops = 0, newton_steps = 1)
# max(m1$gradients)

plot_map <- function(dat, column) {
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_tile(aes(height = 500, width = 500)) +
    facet_wrap(~year) +
    coord_fixed() + 
    ggdark::dark_theme_void()
}


plot_map(predictions$data, "exp(est)") +
  scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Prediction (fixed effects + all random effects)")

plot_map(predictions$data, "exp(est_non_rf)") +
  ggtitle("Prediction (fixed effects only)") +
  scale_fill_viridis_c(trans = "sqrt")

plot_map(predictions$data, "omega_s") +
  ggtitle("Spatial random effects only") +
  scale_fill_gradient2()

plot_map(predictions$data, "epsilon_st") +
  ggtitle("Spatiotemporal random effects only") +
  scale_fill_gradient2()

# not bias correcting for vignette-building speed:
index <- get_index(predictions, bias_correct = FALSE)

ggplot(index, aes(year, est)) + geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  xlab('Year') + ylab('Biomass estimate (metric tonnes)')

mutate(index, cv = sqrt(exp(se^2) - 1)) %>% 
  select(-log_est, -max_gradient, -bad_eig, -se) %>%
  knitr::kable(format = "pandoc", digits = c(0, 0, 0, 0, 2))
