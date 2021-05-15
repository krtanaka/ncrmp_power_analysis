#############################################################################
### Test sampling design of multiple surveys using a stratified analysis  ###
#############################################################################

library(SimSurvey)
library(dplyr)

rm(list = ls())

load("data/survey_grid_kt.RData")

source('scripts/key_functions.R')

options(scipen = 999)
options(digits = 2)

set.seed(438)

plot(survey_grid)
plot(survey_grid_kt)

# base populations with MHI survey grids
pop = sim_abundance() %>% 
  sim_distribution(grid = survey_grid_kt)

# Set-up a series of surveys from all combinations of settings supplied
surveys = expand_surveys(
  
  set_den = c(1, 2) / 1000,
  lengths_cap = c(100, 500),
  ages_cap = c(5, 20)
  
)

## This call runs 10 simulations of 8 different surveys over the same
## population, and then runs a stratified analysis and compares true vs
## estimated values. It may take a while to run.

tests = test_surveys(pop, surveys = surveys, n_sims = 1, n_loops = 10, cores = 8)
doParallel::stopImplicitCluster()

plot_total_strat_fan(tests)
plot_length_strat_fan(tests, years = 1:5, lengths = 1:100)
plot_survey_rank(tests)
plot_error_surface(tests, plot_by = "rule")

library(plotly)
tests$total_strat_error %>%
  filter(survey == 8, sim %in% 1:50) %>%
  group_by(sim) %>%
  plot_ly(x = ~year) %>%
  add_lines(y = ~I_hat, alpha = 0.5, name = "estimated") %>%
  add_lines(y = ~I, color = I("black"), name = "true") %>%
  layout(xaxis = list(title = "Year"),
         yaxis = list(title = "Abundance index"))

plot_total_strat_fan(tests, surveys = 1:8)
plot_length_strat_fan(tests, surveys = 1:8)
plot_age_strat_fan(tests, surveys = 1:8)
plot_age_strat_fan(tests, surveys = 1:8, select_by = "age")

plot_error_surface(tests, plot_by = "rule")
plot_error_surface(tests, plot_by = "samples")

plot_survey_rank(tests, which_strat = "length")
plot_survey_rank(tests, which_strat = "age")

