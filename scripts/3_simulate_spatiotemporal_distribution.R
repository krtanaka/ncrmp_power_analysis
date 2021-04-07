##################################################
### Simulate spatial and temporal distribution ###
### change age-year-space covariance           ###
##################################################

library(SimSurvey)
library(patchwork)

rm(list = ls())

load("data/survey_grid_kt.RData")

set.seed(438)

### simulate a dynamic age structured populations ###
### a long-lived, spatially diffused              ###
long_diffused <- sim_abundance(ages = 1:500, 
                               years = 2000:2020,
                               R = sim_R(log_mean = log(3000),
                                         log_sd = 0.1),
                               Z = sim_Z(log_mean = log(0.05),
                                         log_sd = 0.1,
                                         phi_age = 0.1,
                                         phi_year = 0.1)) %>% 
  sim_distribution(grid = survey_grid_kt,
                   ays_covar = sim_ays_covar(range = 300, 
                                             phi_age = 0.9, 
                                             phi_year = 0.9),
                   depth_par = sim_parabola(mu = 10, 
                                            sigma = 10))

### simulate a dynamic age structured populations ###
### a short-lived, spatially clustered            ###
short_clustered <- sim_abundance(ages = 1:5, 
                                 years = 2000:2020,
                                 R = sim_R(log_mean = log(10000),
                                           log_sd = 0.1),
                                 Z = sim_Z(log_mean = log(0.8))) %>%
  sim_distribution(grid = survey_grid_kt,
                   ays_covar = sim_ays_covar(range = 50, 
                                             phi_year = 0.2,
                                             phi_age = 0.2),
                   depth_par = sim_parabola(mu = 10, sigma = 10))

sim = long_diffused
sim = short_clustered

head(sim$sp_N)
head(sim$grid_xy)

df = merge(sim$sp_N, sim$grid_xy)
df = df %>% group_by(x, y, year) %>% summarise(n = mean(N))

space =
  df %>% 
  ggplot(aes(x, y, fill = n, color = n)) + 
  # geom_tile() +
  geom_point(size = 1.5, alpha = 0.8) +
  facet_wrap(~year) + 
  scale_fill_viridis_c("g/sq.m", limits = c(0,  quantile(df$n, prob = 0.9))) +
  scale_color_viridis_c("g/sq.m", limits = c(0,  quantile(df$n, prob = 0.9))) +
  coord_fixed() + 
  ggdark::dark_theme_void()

time = sim$sp_N %>% 
  group_by(year) %>% 
  summarise(n = sum(N)) %>% 
  ggplot(aes(year, n, color = n)) + 
  geom_line() + 
  geom_point(size = 5) + 
  scale_color_viridis_c("")+ 
  ggdark::dark_theme_classic()

parabola_fun = sim_parabola(mu = 10, sigma = 10)
depth = parabola_fun(x = 0:30) %>% as.data.frame()
depth = cbind(depth, c(0:30))
colnames(depth) = c("response", "depth")
response = depth %>% 
  ggplot(aes(depth, response, color = response)) + 
  geom_point(size = 5) + 
  scale_color_viridis_c("")+ 
  ggdark::dark_theme_classic()

space + (time/response)
