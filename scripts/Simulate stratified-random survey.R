# Simulate stratified-random survey
sim <- sim_abundance(ages = 1:10, years = 1:5) %>%
  sim_distribution(grid = survey_grid_kt) %>%
  sim_survey(n_sims = 5, q = sim_logistic(k = 2, x0 = 3))

plot_survey(sim, which_year = 2, which_sim = 1)

sim <- sim_abundance(ages = 1:5,
                     years = 1:10,
                     Z = sim_Z(log_mean = log(0.5), 
                               log_sd = 0.1, 
                               phi_age = 0.9, 
                               phi_year = 0.9),
                     R = sim_R(log_mean = log(150), 
                               log_sd = 0.1,
                               random_walk = TRUE),
                     N0 = sim_N0(N0 = 1),
                     growth = sim_vonB())%>%
  sim_distribution(grid = survey_grid,
                   ays_covar = sim_ays_covar(phi_age = 0.8,
                                             phi_year = 0.1),
                   depth_par = sim_parabola(alpha = 25, 
                                            mu = 15,
                                            sigma = 5)) %>%
  sim_survey(n_sims = 1, 
             q = sim_logistic(k = 2, x0 = 3))

plot_survey(sim, which_year = 2, which_sim = 1)


# Simulate stratified-random survey
sim <- sim_abundance(ages = 1:10, years = 1:5) %>%
  sim_distribution(grid = survey_grid)

# Simulate stratified-random survey
sim_kt <- sim_abundance(ages = 1:10, years = 1:5) %>%
  sim_distribution(grid = survey_grid_kt)

simsurvey(sim)