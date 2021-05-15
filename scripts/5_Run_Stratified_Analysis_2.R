##########################################
### Simulate stratified-random survey  ###
##########################################

library(SimSurvey)
library(raster)
library(data.table)

rm(list = ls())

load("data/survey_grid_kt.RData")

sim <- sim_abundance(ages = 1:2, years = 1:3) %>%
  sim_distribution(grid = survey_grid_kt) %>%
  sim_survey(n_sims = 1)
plot_survey(sim, which_year = 2, which_sim = 1)

sim3a %>% strat_error()
sim3b %>% strat_error()

