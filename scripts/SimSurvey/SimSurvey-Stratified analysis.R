# Run stratified analysis on simulated data

sim <- sim_abundance(ages = 1:5, years = 1:5,
                     R = sim_R(log_mean = log(1e+7)),
                     growth = sim_vonB(length_group = 1)) %>%
  sim_distribution(grid = survey_grid,
                   ays_covar = sim_ays_covar(sd = 1)) %>%
  sim_survey(n_sims = 1, q = sim_logistic(k = 2, x0 = 3)) %>%
  run_strat()

