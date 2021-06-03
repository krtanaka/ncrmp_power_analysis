# Simulate stratified-random survey
sim <- sim_abundance(ages = 1:5, 
                     years = 2000:2020,
                     R = sim_R(log_mean = log(500)),
                     Z = sim_Z(log_mean = log(0.8))) %>%
  sim_distribution(grid = survey_grid_kt,
                   ays_covar = sim_ays_covar(),
                   depth_par = sim_parabola(mu = 15, sigma = 10)) %>%
  sim_survey(trawl_dim = c(0.005, 0.005))%>%
  run_strat() %>%
  strat_error()

plot_error_surface(sim)

