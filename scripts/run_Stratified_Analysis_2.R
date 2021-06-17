##########################################
### Simulate stratified-random survey  ###
##########################################

library(SimSurvey)
library(raster)
library(data.table)
library(ggplot2)

rm(list = ls())

load("data/survey_grid_kt.RData")

n_sims = 10
min_sets = 10
set_den = 5/1000

options(scipen = 999, digits = 2)
set.seed(100)

sim = sim_abundance(years = 1:10, ages = 1:5) %>% 
  # sim_distribution(grid = make_grid(res = c(100, 100))) %>%
  sim_distribution(grid = survey_grid_kt) %>%
  sim_survey(trawl_dim = c(0.01, 0.0353), n_sims = n_sims, min_sets = min_sets, set_den = set_den) %>% 
  run_strat() %>%
  strat_error()
sim$total_strat_error_stats
sim$total_strat_error
df = sim$total_strat_error

df %>% 
  ggplot() + 
  # geom_point(aes(year, I_hat, color = factor(sim), alpha = 0.5), show.legend = F) + 
  geom_line(aes(year, I_hat, color = factor(sim), alpha = 0.5), show.legend = F) +
  geom_point(aes(year, I), size = 5) + 
  geom_line(aes(year, I), size = 2) + 
  scale_fill_viridis_d() + 
  scale_color_viridis_d() +  
  ggdark::dark_theme_classic()
  
pop = sim_abundance(years = 1:5, ages = 1:5) %>% sim_distribution(grid = survey_grid_kt)

surveys = expand_surveys(set_den = c(0.5, 1, 2, 5, 10, 50)/100,
                         lengths_cap = c(5),
                         ages_cap = c(2))

tests = test_surveys(pop, surveys = surveys, n_sims = 2, n_loops = 5, cores = 6)

plot_total_strat_fan(tests, surveys = 6)
plot_length_strat_fan(tests)
plot_age_strat_fan(tests)

plot_survey_rank(tests, which_strat = "total")
