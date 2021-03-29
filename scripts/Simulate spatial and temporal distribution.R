library(patchwork)
# Simulate spatial and temporal distribution
# Simulate age-year-space covariance

sim <- sim_abundance(ages = 1:20, 
                     years = 1980:2020,
                     R = sim_R(log_mean = log(3e+07)),
                     Z = sim_Z(log_mean = log(0.2))) %>%
  sim_distribution(grid = survey_grid_kt,
                   ays_covar = sim_ays_covar(range = 8000, 
                                             phi_year = 0.99),
                   depth_par = sim_parabola(mu = 10, sigma = 10))

sim <- sim_abundance(ages = 1:5, 
                     years = 2000:2020,
                     R = sim_R(log_mean = log(1e+10)),
                     Z = sim_Z(log_mean = log(0.8))) %>%
  sim_distribution(grid = survey_grid_kt,
                   ays_covar = sim_ays_covar(range = 1000, 
                                             phi_year = 0.1),
                   depth_par = sim_parabola(mu = 10, sigma = 10))

# plot_distribution(sim)

head(sim$sp_N)
head(sim$grid_xy)

df = merge(sim$sp_N, sim$grid_xy)
df = df %>% group_by(x, y, year) %>% summarise(n = mean(N))

space = df %>% ggplot(aes(x, y, fill = n, color = n)) + 
  geom_raster() + 
  geom_point(size = 0.5) +
  facet_wrap(~year) + 
  scale_fill_viridis_c("g/sq.m") + 
  scale_color_viridis_c("g/sq.m") + 
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
  

#Define relationships with covariates
parabola_fun <- sim_parabola(mu = 15, sigma = 5, plot = TRUE)
parabola_fun(x = 0:30)

parabola_fun <- sim_parabola(mu = 50, sigma = 0.5, plot = TRUE)
parabola_fun(x = 0:100)

parabola_fun <- sim_parabola(mu = log(40), sigma = 0.5, log_space = TRUE, plot = TRUE)
parabola_fun(x = 1:1000)

parabola_fun <- sim_parabola(mu = c(50, 120), sigma = c(5, 3), plot = TRUE)
parabola_fun(x = rep(1:200, 2), age = rep(c(1, 2), each = 200))


#Make a depth stratified survey grid
r <- make_grid(res = c(10, 10),
               shelf_depth = 200,
               shelf_width = 500,
               depth_range = c(0, 1000),
               strat_breaks = seq(0, 1000, by = 40),
               strat_splits = 3,
               n_div = 5)
raster::plot(r)

p <- raster::rasterToPolygons(r$strat, dissolve = TRUE)
sp::plot(p)


