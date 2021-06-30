##########################################
### Simulate stratified-random survey  ###
##########################################

library(SimSurvey)
library(raster)
library(data.table)
library(ggplot2)
library(dplyr)
library(scales)

rm(list = ls())

# set.seed(10)

islands = c("Hawaii", "Kahoolawe", "Kauai", "Lanai", "Maui", "Molokai", "Niihau", "Oahu" )[sample(1:8, 1)]
load(paste0("data/survey_grid_", islands, ".RData"))

n_sims = 10
min_sets = 1
set_den = 2/1000
trawl_dim = c(0.01, 0.0353)
resample_cells = FALSE

# options(scipen = 999, digits = 2)

sim = sim_abundance(years = 2010:2020, ages = 1:5,
                    R = sim_R(log_mean = log(100),
                              log_sd = 0.9),
                    Z = sim_Z(log_mean = log(0.1))) %>% 
  sim_distribution(grid = survey_grid_kt) 

strat_sets <- cell_sets <- NULL

cells <- data.table(rasterToPoints(sim$grid))

strat_det <- cells[, list(strat_cells = .N), by = "strat"]
strat_det$tow_area <- prod(trawl_dim)
strat_det$cell_area <- prod(res(sim$grid))
strat_det$strat_area <- strat_det$strat_cells * prod(res(sim$grid))
strat_det$strat_sets <- round(strat_det$strat_area * set_den)
strat_det$strat_sets[strat_det$strat_sets < min_sets] <- min_sets
strat_det$area_prop = rescale(strat_det$strat_area, to = c(0.1, 0.9))

cells <- merge(cells, strat_det, by = c("strat"))

i <- rep(seq(nrow(cells)), times = length(sim$years))
y <- rep(sim$years, each = nrow(cells))

cells <- cells[i, ]
cells$year <- y

i <- rep(seq(nrow(cells)), times = n_sims)
s <- rep(seq(n_sims), each = nrow(cells))

cells <- cells[i, ]
cells$sim <- s

sets = sample(cells)

# proportional sampling
subsets_prop <- cells[, .SD[sample(.N, strat_area, replace = resample_cells, prob = strat_area)], 
              by = c("sim", "year", "strat")]

subsets_prop[, `:=`(cell_sets, .N), by = c("sim", "year", "cell")]
subsets_prop$set <- seq(nrow(subsets_prop))
subsets_prop$design = "proportional"

#equal sampling
subsets_norm <- cells[, .SD[sample(.N, strat_sets, replace = resample_cells)], 
              by = c("sim", "year", "strat")]

subsets_norm[, `:=`(cell_sets, .N), by = c("sim", "year", "cell")]
subsets_norm$set <- seq(nrow(subsets_norm))
subsets_norm$design = "equal"

sim_sets = rbind(subsets_norm, subsets_prop)

sim_sets = sim_sets %>% 
  group_by(x, y, design) %>% 
  summarise(strat = mean(strat))

sim_sets %>% 
  # subset(strat == 1) %>%
  ggplot(aes(x, y, color = factor(strat))) + 
  geom_point() + 
  coord_fixed() + 
  scale_color_discrete("") + 
  facet_wrap(.~design)

sim_sets %>% 
  group_by(strat, design) %>% 
  summarise(n = n())

