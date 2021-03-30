# Simulate stratified-random survey
sim <- sim_abundance(ages = 1:5, 
                     years = 2000:2020,
                     R = sim_R(log_mean = log(1e+10)),
                     Z = sim_Z(log_mean = log(0.8))) %>%
  sim_distribution(grid = survey_grid_kt,
                   ays_covar = sim_ays_covar(range = 8000, 
                                             phi_year = 0.99),
                   depth_par = sim_parabola(mu = 100, sigma = 10)) %>%
  sim_survey(trawl_dim = c(0.005, 0.005))

plot_survey(sim, which_year = 2000)

