##########################################
### Simulate stratified-random survey  ###
##########################################

library(SimSurvey)
library(raster)
library(data.table)

rm(list = ls())

load("data/survey_grid_kt.RData")

source('scripts/key_functions.R')

options(scipen=999)
options(digits = 2)


set.seed(438)

sim <- sim_abundance(ages = 1:5, 
                     years = 2018:2020,
                     R = sim_R(log_mean = log(3e+07),
                               log_sd = 0.1),
                     Z = sim_Z(log_mean = log(0.8)),
                     growth = sim_vonB(Linf = 30)) %>%
  sim_distribution(grid = survey_grid_kt,
                   ays_covar = sim_ays_covar(range = 300, 
                                             phi_year = 0.2,
                                             phi_age = 0.2),
                   depth_par = sim_parabola(mu = 10, sigma = 10)) %>% 
  sim_survey_rea(min_sets = 10, 
                 ages_cap = 500,
                 lengths_cap = 100,
                 set_den = 2/1000) %>% 
  run_strat()

sim %>% strat_error()

pop = sim_abundance(ages = 1:5, years = 1:5) %>% sim_distribution()
surveys = expand_surveys(set_den = c(0.5, 1)/1000,
                         lengths_cap = c(5, 10),
                         ages_cap = c(2, 5))
tests = test_surveys(pop, surveys = surveys, n_sims = 1, n_loops = 10)
plot_total_strat_fan(tests)
plot_length_strat_fan(tests, years = 1:5, lengths = 1:100)
plot_survey_rank(tests)
plot_error_surface(tests, plot_by = "rule")
