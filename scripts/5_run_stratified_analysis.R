##########################################
### Simulate stratified-random survey  ###
##########################################

library(SimSurvey)
library(raster)
library(data.table)

rm(list = ls())

load("data/survey_grid_kt.RData")

source('scripts/key_functions.R')

set.seed(438)

sim <- sim_abundance(ages = 1:5, 
                     years = 2018:2020,
                     R = sim_R(log_mean = log(100000),
                               log_sd = 0.1),
                     Z = sim_Z(log_mean = log(0.8)),
                     growth = sim_vonB(Linf = 30)) %>%
  sim_distribution(grid = survey_grid_kt,
                   ays_covar = sim_ays_covar(range = 300, 
                                             phi_year = 0.2,
                                             phi_age = 0.2),
                   depth_par = sim_parabola(mu = 10, sigma = 10)) %>% 
  sim_survey_rea(min_sets = 100, 
                 ages_cap = 100,
                 lengths_cap = 100,
                 set_den = 2/1000) %>% 
  run_strat()

