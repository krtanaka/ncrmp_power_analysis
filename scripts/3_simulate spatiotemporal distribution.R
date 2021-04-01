##################################################
### Simulate spatial and temporal distribution ###
### change age-year-space covariance           ###
##################################################

library(SimSurvey)
library(patchwork)

rm(list = ls())

load("data/survey_grid_kt.RData")

sim1 <- sim_abundance(ages = 1:20, 
                     years = 1980:2020,
                     R = sim_R(log_mean = log(100)),
                     Z = sim_Z(log_mean = log(0.2))) %>%
  sim_distribution(grid = survey_grid_kt,
                   ays_covar = sim_ays_covar(range = 50, 
                                             phi_year = 0.9),
                   depth_par = sim_parabola(mu = 10, sigma = 10))

sim2 <- sim_abundance(ages = 1:5, 
                     years = 2000:2020,
                     R = sim_R(log_mean = log(500)),
                     Z = sim_Z(log_mean = log(0.8))) %>%
  sim_distribution(grid = survey_grid_kt,
                   ays_covar = sim_ays_covar(range = 2, 
                                             phi_year = 0.1,
                                             phi_age = 0.8,),
                   depth_par = sim_parabola(mu = 10, sigma = 10))

sim = sim1
sim = sim2

plot_distribution(sim)

head(sim$sp_N)
head(sim$grid_xy)

df = merge(sim$sp_N, sim$grid_xy)
df = df %>% group_by(x, y, year) %>% summarise(n = mean(N))

space = df %>% 
  ggplot(aes(x, y, fill = n)) + 
  geom_tile(aes(width = 0.05, height = 0.05)) + 
  facet_wrap(~year) + 
  scale_fill_viridis_c("g/sq.m") + 
  coord_fixed() + 
  ggdark::dark_theme_void()

time = sim$sp_N %>% 
  group_by(year) %>% 
  summarise(n = sum(N)) %>% 
  ggplot(aes(year, n, color = n)) + 
  geom_line() + 
  geom_point(size = 2) + 
  scale_color_viridis_c("")+ 
  ggdark::dark_theme_classic()

parabola_fun = sim_parabola(mu = 10, sigma = 10)
depth = parabola_fun(x = 0:30) %>% as.data.frame()
depth = cbind(depth, c(0:30))
colnames(depth) = c("response", "depth")
response = depth %>% 
  ggplot(aes(depth, response, color = response)) + 
  geom_point() + 
  scale_color_viridis_c("")+ 
  ggdark::dark_theme_classic()
  
space + (time/response)
  
