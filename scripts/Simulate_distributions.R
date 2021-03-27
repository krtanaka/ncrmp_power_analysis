
parabola_fun <- sim_parabola(mu = 200, sigma = 50, plot = TRUE)
parabola_fun(x = 0:500)

#long lived, sedentary species 
sim <- sim_abundance(ages = 1:20, 
                     R = sim_R(log_mean = log(3e+07)),
                     Z = sim_Z(log_mean = log(0.2))) %>% 
  sim_distribution(
    grid = make_grid(res = c(20, 20)),
    # grid = survey_grid, 
    ays_covar = sim_ays_covar(range = 8000, 
                              phi_year = 0.9),
    depth_par = sim_parabola(mu = 500, sigma = 70))

df = merge(sim$sp_N, sim$grid_xy)
df = df %>% group_by(x, y, year) %>% summarise(n = mean(N))

df %>% ggplot(aes(x, y, fill = n, color = n)) + 
  geom_raster() + 
  scale_fill_viridis_c() + 
  facet_wrap(~year) 

#short lived, dynamic species 
sim <- sim_abundance(ages = 1:5, 
                     R = sim_R(log_mean = log(1e+10)),
                     Z = sim_Z(log_mean = log(0.8))) %>%
  sim_distribution(grid = make_grid(res = c(20, 20)),
                   ays_covar = sim_ays_covar(range = 1000, 
                                             phi_year = 0.1))

df = merge(sim$sp_N, sim$grid_xy)
df = df %>% group_by(x, y, year) %>% summarise(n = mean(N))

df %>% ggplot(aes(x, y, fill = n, color = n)) + 
  geom_raster() + 
  scale_fill_viridis_c() + 
  facet_wrap(~year) 
