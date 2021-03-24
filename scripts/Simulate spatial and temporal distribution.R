# Simulate spatial and temporal distribution

sim <- sim_abundance(ages = 1:2,
                     years = 1:10) %>%
  sim_distribution(grid = survey_grid,
                   ays_covar = sim_ays_covar(phi_age = 0.8,
                                             phi_year = 0.1),
                   depth_par = sim_parabola(mu = 200,
                                            sigma = 50))
head(sim$sp_N)
head(sim$grid_xy)

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


#Simulate age-year-space covariance

#Define relationships with covariates
parabola_fun <- sim_parabola(alpha = 25, mu = 50, sigma = 5, plot = TRUE)
parabola_fun(x = 0:100)

df = merge(sim$sp_N, sim$grid_xy)
df = df %>% group_by(x, y, year) %>% summarise(n = mean(N))
df %>% ggplot(aes(x, y, fill = n)) + geom_raster() + facet_wrap(~year) + scale_fill_viridis_c()
