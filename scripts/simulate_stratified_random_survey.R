##########################################
### Simulate stratified-random survey  ###
##########################################

library(SimSurvey)

rm(list = ls())

load("data/survey_grid_kt.RData")

source('scripts/key_functions.R')

low_effort = sim_abundance() %>% 
  sim_distribution() %>% 
  sim_survey(n_sims = 1, 
             set_den = 1/1000,
             lengths_cap = 100, 
             ages_cap = 5)

high_effort = sim_abundance() %>% 
  sim_distribution() %>% 
  sim_survey(n_sims = 1, 
             set_den = 5/1000,
             lengths_cap = 500, 
             ages_cap = 25)

plot_survey(low_effort)
plot_survey(high_effort)

set.seed(438)

### simulate a dynamic age structured populations ###
### a long-lived, spatially diffused              ###
long_diffused <- sim_abundance(ages = 1:500, 
                               years = 2000:2020,
                               R = sim_R(log_mean = log(3000),
                                         log_sd = 0.1),
                               Z = sim_Z(log_mean = log(0.05),
                                         log_sd = 0.1,
                                         phi_age = 0.1,
                                         phi_year = 0.1)) %>% 
  sim_distribution(grid = survey_grid_kt,
                   ays_covar = sim_ays_covar(range = 300, 
                                             phi_age = 0.9, 
                                             phi_year = 0.9),
                   depth_par = sim_parabola(mu = 10, 
                                            sigma = 10)) %>% 
  sim_survey_rea(n_sims = 1, 
                 min_sets = 20, 
                 ages_cap = 10,
                 set_den = 2/1000, 
                 trawl_dim = c(0.001, 0.001))

plot_survey(long_diffused, which_year = 2000)

### simulate a dynamic age structured populations ###
### a short-lived, spatially clustered            ###
short_clustered <- sim_abundance(ages = 1:5, 
                                 years = 2000:2020,
                                 R = sim_R(log_mean = log(100000),
                                           log_sd = 0.1),
                                 Z = sim_Z(log_mean = log(0.8)),
                                 growth = sim_vonB(Linf = 30)) %>%
  sim_distribution(grid = survey_grid_kt,
                   ays_covar = sim_ays_covar(range = 300, 
                                             phi_year = 0.2,
                                             phi_age = 0.2),
                   depth_par = sim_parabola(mu = 10, sigma = 10)) %>% 
  sim_survey_rea(  min_sets = 90, 
                   ages_cap = 50,
                   lengths_cap = 10,
                   set_den = 2/1000)

plot_survey(short_clustered, which_year = 2019)
