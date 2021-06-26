##########################################
### Simulate stratified-random survey  ###
##########################################

library(SimSurvey)
library(raster)
library(data.table)
library(ggplot2)
library(dplyr)

rm(list = ls())

set.seed(10)

islands = c("Hawaii", "Kahoolawe", "Kauai", "Lanai", "Maui", "Molokai", "Niihau", "Oahu" )[sample(1:8, 1)]
load(paste0("data/survey_grid_", islands, ".RData"))

n_sims = 1
min_sets = 2
set_den = 2/1000

# options(scipen = 999, digits = 2)

sim = sim_abundance(years = 2010:2020, ages = 1:5,
                    R = sim_R(log_mean = log(100),
                              log_sd = 0.9),
                    Z = sim_Z(log_mean = log(0.1))) %>% 
  sim_distribution(grid = survey_grid_kt) %>% 
  sim_sets(
    n_sims = n_sims,
    trawl_dim = c(0.01, 0.0353),
    min_sets = min_sets,
    set_den = set_den,
    resample_cells = FALSE
  )

sim %>% 
  # subset(strat == 2) %>% 
  ggplot(aes(x, y, color = factor(strat))) + 
  geom_point() + 
  facet_wrap(.~sim)

sim %>% 
  group_by(strat) %>% 
  summarise(n = n())

