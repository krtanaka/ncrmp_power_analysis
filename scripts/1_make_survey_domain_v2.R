########################################
### prep survey_grid for simulations ###
########################################

rm(list = ls())

library(dplyr)
library(ggplot2)
library(raster)
library(sp)
library(rgdal)
library(tidyr)
library(patchwork)
library(SimSurvey)
library(sf)
library(readr)

utm = read_csv('data/misc/ncrmp_utm_zones.csv')

islands = c("gua", "rot", "sai", "tin")[1]

for (isl in 1:length(islands)) {
  
  isl = 1
  
  load(paste0("data/gis_bathymetry/", islands[isl], ".RData"))
  bathymetry = rasterFromXYZ(topo)
  
  df = topo; rm(topo)
  df$cell = as.numeric(1:dim(df)[1])
  
  ###########################################################
  ### import sector/reefzones shapefile                   ###
  ### adjust resolutions and merge with bathymetry data   ###
  ### these are outputs from "convert_shp_to_data.frame.R ###
  ###########################################################
  
  if (isl == 1) {
    load("data/gis_sector/gua.RData"); sector = raster_and_table[[1]]; sector_name = raster_and_table[[2]]
    load("data/gis_reef/gua.RData"); reef = raster_and_table[[1]]; reef_name = raster_and_table[[2]]
    load("data/gis_hardsoft/gua.RData"); hardsoft = raster_and_table[[1]]; hardsoft_name = raster_and_table[[2]]
  }
  # if (isl == 2) { # using reef raster as a placeholder bc there are no sectors for this island
  #   load("data/gis_reef/rot.RData"); sector = raster_and_table[[1]] ; sector_name = raster_and_table[[2]] 
  #   load("data/gis_reef/rot.RData"); reef = raster_and_table[[1]]; reef_name = raster_and_table[[2]]
  #   load("data/gis_hardsoft/rot.RData"); hardsoft = raster_and_table[[1]]; hardsoft_name = raster_and_table[[2]]
  # }
  # if (isl == 3) { # using reef raster as a placeholder bc there are no sectors for this island
  #   load("data/gis_reef/sai.RData"); sector = raster_and_table[[1]] ; sector_name = raster_and_table[[2]] 
  #   load("data/gis_reef/sai.RData"); reef = raster_and_table[[1]]; reef_name = raster_and_table[[2]]
  #   load("data/gis_hardsoft/sai.RData"); hardsoft = raster_and_table[[1]]; hardsoft_name = raster_and_table[[2]]
  # }
  # if (isl == 4) { # using reef raster as a placeholder bc there are no sectors for this island
  #   load("data/gis_reef/tin.RData"); sector = raster_and_table[[1]] ; sector_name = raster_and_table[[2]] 
  #   load("data/gis_reef/tin.RData"); reef = raster_and_table[[1]]; reef_name = raster_and_table[[2]]
  #   load("data/gis_hardsoft/tin.RData"); hardsoft = raster_and_table[[1]]; hardsoft_name = raster_and_table[[2]]
  # }
  
  rm(raster_and_table)
  
  # if (isl %in% c(2:4)) {
  #   
  #   sector$z = 1
  #   sector_name = data.frame(SEC_NAME = "single_sector",
  #                            i = 1)
  # }
  
  crm_res = rasterFromXYZ(df[,c("x", "y", "cell")])
  dim(crm_res); crm_res  
  
  hardsoft = resample(hardsoft, crm_res, method = "bilinear") 
  sector = resample(sector, crm_res, method = "bilinear") 
  reef = resample(reef, crm_res, method = "bilinear") 
  bathymetry = resample(bathymetry, crm_res, method = "bilinear") 
  
  df = stack(hardsoft, sector, reef, bathymetry)
  df = as.data.frame(rasterToPoints(df))
  
  colnames(df) = c("longitude", "latitude", "hardsoft", "sector", "reef", "depth")
  
  df = na.omit(df)
  
  df$cell = 1:dim(df)[1]; df$cell = as.numeric(df$cell)
  df$division = as.numeric(1)
  
  df$depth_bin = ""
  df$depth_bin = ifelse(df$depth <= 0  & df$depth >= -6, 1L, df$depth_bin) 
  df$depth_bin = ifelse(df$depth < -6  & df$depth >= -18, 2L, df$depth_bin) 
  df$depth_bin = ifelse(df$depth < -18, 3L, df$depth_bin) 
  
  df$hardsoft = round(df$hardsoft, 0)
  df$reef = round(df$reef, 0)
  df$sector = round(df$sector, 0)
  
  df$sector = as.character(df$sector)
  df$reef = as.character(df$reef)
  
  df$depth = as.numeric(df$depth*-1)
  
  df$longitude = df$longitude * 0.001
  df$latitude = df$latitude * 0.001
  
  df = as.data.frame(df)
  
  colnames(sector_name) = c("sector", "sector_id")
  colnames(reef_name) = c("reef", "reef_id")
  colnames(hardsoft_name) = c("hardsoft", "hardsoft_id")
  
  df = merge(df, sector_name)
  df = merge(df, reef_name)
  df = merge(df, hardsoft_name)
  
  (depth = df %>% 
      ggplot( aes(longitude, latitude, fill = depth_bin)) + 
      geom_raster() + 
      scale_fill_discrete( "depth_bins") +
      ggdark::dark_theme_minimal())
  
  (sector = df %>% 
      ggplot( aes(longitude, latitude, fill = sector_id)) + 
      geom_raster() +
      scale_fill_discrete("sector") +
      ggdark::dark_theme_minimal())
  
  (reef = df %>% 
      ggplot( aes(longitude, latitude, fill = reef_id)) + 
      geom_raster() +
      scale_fill_discrete("reef") +
      ggdark::dark_theme_minimal())
  
  (hardsoft = df %>% 
      ggplot( aes(longitude, latitude, fill = hardsoft_id)) + 
      geom_raster() +
      scale_fill_discrete("hardsoft") +
      ggdark::dark_theme_minimal())
  
  df = df %>%
    subset(sector_id != "GUA_LAND") %>% # filter sector
    subset(reef_id %in% c("Backreef", "Forereef", "Lagoon")) %>% # filter land and Reef Crest/Reef Flat
    subset(hardsoft_id %in% c("Hard", "Unknown")) # filter hardsoft
  
  df$strat = paste(df$depth_bin, 
                   df$sector,
                   df$reef,
                   sep = "_")
  
  df$strat = as.numeric(as.factor(df$strat))
  
  (strata = df %>% 
      ggplot( aes(longitude, latitude, fill = factor(strat))) + 
      geom_raster() +
      scale_fill_discrete("Strata") +
      ggdark::dark_theme_minimal())
  
  cell = rasterFromXYZ(df[,c("longitude", "latitude", "cell")]); plot(cell)
  division = rasterFromXYZ(df[,c("longitude", "latitude", "division")]); plot(division)
  strat = rasterFromXYZ(df[,c("longitude", "latitude", "strat")]); plot(strat)
  depth = rasterFromXYZ(df[,c("longitude", "latitude", "depth")]); plot(depth)
  
  values = raster::values
  
  survey_grid_ncrmp = stack(cell, division, strat, depth)
  survey_grid_ncrmp$strat = round(survey_grid_ncrmp$strat, digits = 0)
  # values(survey_grid_ncrmp$strat) = ifelse(values(survey_grid_ncrmp$strat) > 3, 3, values(survey_grid_ncrmp$strat))
  values(survey_grid_ncrmp$division) = ifelse(is.na(values(survey_grid_ncrmp$division)), NA, 1)
  
  sp::spplot(survey_grid$cell) #SimSurvey example
  sp::spplot(survey_grid_ncrmp$cell)
  
  sp::spplot(survey_grid$division) #SimSurvey example
  sp::spplot(survey_grid_ncrmp$division)
  
  sp::spplot(survey_grid$strat) #SimSurvey example
  sp::spplot(survey_grid_ncrmp$strat)
  
  sp::spplot(survey_grid$depth) #SimSurvey example
  sp::spplot(survey_grid_ncrmp$depth)
  
  p <- raster::rasterToPolygons(survey_grid$strat, dissolve = TRUE)
  sp::plot(p); axis(1); axis(2)
  
  p <- raster::rasterToPolygons(survey_grid_ncrmp$strat, dissolve = TRUE)
  sp::plot(p); axis(1); axis(2)
  
  survey_grid_ncrmp = readAll(survey_grid_ncrmp)
  
  save(survey_grid_ncrmp, file = paste0("data/survey_grid_w_sector_reef/survey_grid_", islands[isl], ".RData"))

}
