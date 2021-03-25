# Simulate spatial and temporal distribution
#Simulate age-year-space covariance

sim <- sim_abundance(ages = 1:2,
                     years = 1980:2020,
                     Z = sim_Z(log_mean = log(0.5), 
                               log_sd = 0.1, 
                               phi_age = 0.9, 
                               phi_year = 0.9),
                     R = sim_R(log_mean = log(150), 
                               log_sd = 0.1,
                               random_walk = TRUE),
                     N0 = sim_N0(),
                     growth = sim_vonB())%>%
  sim_distribution(grid = survey_grid_kt,
                   ays_covar = sim_ays_covar(phi_age = 0.8,
                                             phi_year = 0.1),
                   depth_par = sim_parabola(mu = 10,
                                            sigma = 50))

plot_distribution(sim)

head(sim$sp_N)
head(sim$grid_xy)

df = merge(sim$sp_N, sim$grid_xy)
df = df %>% group_by(x, y, year) %>% summarise(n = mean(N))

df %>% ggplot(aes(x, y, fill = n, color = n)) + 
  geom_raster() + 
  geom_point(size = 0.5) +
  facet_wrap(~year) + 
  scale_fill_viridis_c("g/sq.m") + 
  scale_color_viridis_c("g/sq.m") + 
  coord_fixed() + 
  ggdark::dark_theme_void()

#Define relationships with covariates
parabola_fun <- sim_parabola(alpha = 25, mu = 15, sigma = 5, plot = TRUE)
parabola_fun(x = 0:30)


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


