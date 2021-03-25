# Simulate stratified-random survey
sim <- sim_abundance(ages = 1:10, years = 1:5) %>%
sim_distribution(grid = survey_grid_kt) %>%
  sim_survey(n_sims = 5, 
             q = sim_logistic(k = 2, x0 = 3))

plot_survey(sim, which_year = 2, which_sim = 1)


