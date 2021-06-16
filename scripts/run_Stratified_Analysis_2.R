##########################################
### Simulate stratified-random survey  ###
##########################################

library(SimSurvey)
library(raster)
library(data.table)

rm(list = ls())

load("data/survey_grid_kt.RData")

options(scipen = 999, digits = 2)
set.seed(438)

sim <- sim_abundance() %>%
  sim_distribution(grid = make_grid(res = c(10, 10))) %>%
  # sim_distribution(grid = survey_grid_kt) %>% 
  sim_survey()

plot_survey(sim, which_year = 2, which_sim = 1)

sim_abundance(years = 1:5, ages = 1:5) %>% 
  # sim_distribution(grid = make_grid(res = c(10, 10))) %>% 
  sim_distribution(grid = survey_grid_kt) %>% 
  sim_survey(min_sets = 10, set_den = 2/1000) %>% 
  run_strat()
  


pop = sim_abundance(years = 1:5, ages = 1:5) %>% sim_distribution(grid = survey_grid_kt)

surveys = expand_surveys(set_den = c(0.5, 1, 2, 5, 10, 50)/100,
                         lengths_cap = c(5),
                         ages_cap = c(2))

tests = test_surveys(pop, surveys = surveys, n_sims = 2, n_loops = 5, cores = 6)

plot_total_strat_fan(tests, surveys = 6)
plot_length_strat_fan(tests)
plot_age_strat_fan(tests)

plot_survey_rank(tests, which_strat = "total")
