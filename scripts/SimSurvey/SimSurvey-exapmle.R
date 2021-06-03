library(SimSurvey)

make_grid(
  x_range = c(-140, 140),
  y_range = c(-140, 140),
  res = c(1, 1),
  shelf_depth = 200,
  shelf_width = 100,
  depth_range = c(0, 30),
  n_div = 1,
  strat_breaks = seq(0, 30, by = 3),
  strat_splits = 1,
  method = "spline"
)

r <- make_grid()
raster::plot(r)

p <- raster::rasterToPolygons(r$strat, dissolve = TRUE)
sp::plot(p)

sim <- sim_abundance(ages = 1:10, years = 1:5) %>%
  sim_distribution(grid = make_grid(res = c(10, 10))) %>%
  sim_survey(n_sims = 5, q = sim_logistic(k = 2, x0 = 3))

plot_survey(sim, which_year = 2, which_sim = 1)

