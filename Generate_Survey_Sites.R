###############################################
### generate stratified-random survey sites ###
###############################################

library(SimSurvey)
library(raster)
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggdark)
library(colorRamps)

rm(list = ls())

load("data/modeled_survey_variability.RData") #modeled at grid scale

# set.seed(42)
# options(scipen = 999, digits = 2)

# pick an island
island = c("gua", "rot", "sai", "tin")[sample(1:4, 1)]; print(island)

# survey domain with sector & reef & hard_unknown & 3 depth bins
load(paste0("data/survey_grid_w_sector_reef/survey_grid_", island, ".RData")) 

# simulate stratified random surveys --------------------------------------

# n_sims = 100 # number of simulations
total_sample = 30 # total sample efforts you want to deploy
min_sets = 1 # minimum number of sets per strat
trawl_dim = c(0.01, 0.0353) # 0.000353 sq.km (353 sq.m) from two 15-m diameter survey cylinders
resample_cells = F

n <- id <- division <- strat <- N <- strat_sets <- cell_sets <- NULL

cells <- data.table(rasterToPoints(survey_grid_ncrmp))

# add modeled trophic biomass variability, summarize by strata
cells$sd = predict(g, cells); sd = cells[,c("strat", "sd")]; sd = sd %>% group_by(strat) %>% summarise(sd = mean(sd,na.rm = T))

strat_det <- cells[, list(strat_cells = .N), by = "strat"]; strat_det
strat_det$tow_area <- prod(trawl_dim); strat_det
strat_det$cell_area <- prod(res(survey_grid_ncrmp)); strat_det
strat_det$strat_area <- strat_det$strat_cells * prod(res(survey_grid_ncrmp)); strat_det
strat_det = right_join(strat_det, sd); strat_det

## allocate sampling units by area * sd
strat_det$weight = strat_det$strat_area * strat_det$sd; strat_det
strat_det$strat_sets = round((total_sample * strat_det$weight) / sum(strat_det$weight), 0); strat_det

# allocate sampling units by area
# strat_det$strat_sets = round((total_sample * strat_det$strat_area) / sum(strat_det$strat_area), 0); strat_det

# make sure minimum number of sets per strat is not 0 or 1
strat_det$strat_sets[strat_det$strat_sets < min_sets] <- min_sets; strat_det
strat_table = strat_det %>% dplyr::select(strat, strat_sets); strat_table

# add "strat" "strat_cells" "tow_area" ...
strat_det = strat_det[,c("strat", "strat_cells", "tow_area", "cell_area", "strat_area", "strat_sets" )]
cells <- merge(cells, strat_det, by = c("strat")) 

utmcoor <- SpatialPoints(cbind(cells$x, cells$y), proj4string = CRS("+proj=utm +units=km +zone=55"))
longlatcoor <- spTransform(utmcoor,CRS("+proj=longlat"))
cells$longitude <- coordinates(longlatcoor)[,1]
cells$latitude <- coordinates(longlatcoor)[,2]

# subset "cells" to create site locations
sets <- cells[, .SD[sample(.N, strat_sets, replace = resample_cells)], 
              by = c("strat")]

# count number of distinct sim*year*cell combinations
sets[, `:=`(cell_sets, .N), by = c("cell")]
sets$set <- seq(nrow(sets))
sets 

(bathymetry = cells %>% 
    ggplot(aes(x, y)) +
    geom_raster(aes(fill = depth)) + 
    scale_fill_viridis_c("Depth (m)") + 
    ylab("Northings (km)") + xlab("Eastings (km)") + 
    coord_fixed() +
    theme_minimal() + 
    theme(legend.position = "right"))

(strata = cells %>% 
    ggplot(aes(x, y)) +
    geom_raster(aes(fill = factor(strat))) + 
    scale_fill_viridis_d("Strata") + 
    ylab("Northings (km)") + xlab("Eastings (km)") + 
    coord_fixed() +
    theme_minimal() + 
    theme(legend.position = "right"))

(site_location = sets %>% 
    ggplot(aes(x, y, color = factor(strat))) + 
    geom_point() +   
    coord_fixed() +
    ylab("Northings (km)") + xlab("Eastings (km)") +
    scale_colour_viridis_d("Strata") + 
    theme_minimal() + 
    theme(legend.position = "right")+ 
    labs(
      title = "",
      subtitle = paste0(paste0("Number of strata = ", length(unique(cells$strat)), "\n", 
                               "Survey effort = ", total_sample, " sites"))))

(site_density = cells %>% 
    ggplot(aes(x, y)) +
    coord_fixed() +
    geom_raster(aes(fill = factor(strat_sets))) + 
    ylab("Northings (km)") + xlab("Eastings (km)") +
    # scale_fill_viridis_c("Site_allocation"))
    theme_minimal() + 
    theme(legend.position = "right"))
