##########################################
### Simulate stratified-random survey  ###
##########################################

library(SimSurvey)
library(raster)
library(data.table)

rm(list = ls())

load("data/survey_grid_kt.RData")

options(scipen = 999, digits = 2)

sim <- sim_abundance(ages = 1:2, years = 1:5) %>%
  sim_distribution(grid = survey_grid_kt) %>%
  sim_survey(n_sims = 2,
             min_sets = 3)
plot_survey(sim, which_year = 1, which_sim = 1)

sim %>% strat_error()
sim3b %>% strat_error()

