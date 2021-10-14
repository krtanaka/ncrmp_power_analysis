###################################################
### prep survey_grid for subsequent simulations ###
###################################################

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

islands = c("gua", "rot", "sai", "tin")

for (isl in 1:length(islands)) {
  
  isl = 2
  
  load(paste0("data/gis_bathymetry/raster/", islands[isl], ".RData"))
  
  df = topo; rm(topo)
  
  # change to 1000 m res
  df$longitude = df$x
  df$latitude = df$y
  
  df = df %>% 
    group_by(longitude, latitude) %>% 
    summarise(depth = mean(depth))
  
  df$cell = 1:dim(df)[1]; df$cell = as.numeric(df$cell)
  df$division = as.numeric(1)
  
  ###########################################################
  ### import sector/reefzones shapefile                   ###
  ### adjust resolutions and merge with bathymetry data   ###
  ### these are outputs from "convert_shp_to_data.frame.R ###
  ###########################################################
  
  if (isl == 1) {
    load("data/gis_sector/raster/gua.RData"); sector = raster_and_table[[1]]; sector_name = raster_and_table[[2]]
    load("data/gis_reef/raster/gua.RData"); reef = raster_and_table[[1]]; reef_name = raster_and_table[[2]]
  }
  if (isl == 2) {
    # load("data/gis_sector/raster/kau.RData"); sector = raster_and_table[[1]]
    load("data/gis_reef/raster/rot.RData"); reef = raster_and_table[[1]]
  }
  if (isl == 3) {
    # load("data/gis_sector/raster/lan.RData"); sector = raster_and_table[[1]]
    load("data/gis_reef/raster/sai.RData"); reef = raster_and_table[[1]]
  }
  if (isl == 4) {
    # load("data/gis_sector/raster/timai.RData"); sector = raster_and_table[[1]]
    load("data/gis_reef/raster/tin.RData"); reef = raster_and_table[[1]]
  }
 
  rm(raster_and_table)
  
  sector = rasterToPoints(sector) %>% as.data.frame(); colnames(sector) = c("x", "y", "z")
  reef = rasterToPoints(reef) %>% as.data.frame(); colnames(reef) = c("x", "y", "z")
  
  # merge sectors -----------------------------------------------------------
 
  sector$lon = sector$x
  sector$lat = sector$y
  
  sector$sector_name = as.numeric(as.factor(sector$z))
  sector = as.matrix(sector[,c("lon", "lat", "sector_name")])
  e = extent(sector[,1:2])
  
  crm_res = rasterFromXYZ(df[,c("longitude", "latitude", "cell")])
  dim(crm_res); crm_res
  
  res = 10  # rasterize it, but be careful with resolutions, lower = better but more missing points
  r <- raster(e, ncol = round((dim(crm_res)[2]/res), digits = 0), nrow = round(dim(crm_res)[1]/res, digits = 0))
  sector <- rasterize(sector[, 1:2], r, sector[,3], fun = mean)
  dim(sector)
  sector
  default_proj = "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
  crs(sector) = default_proj
  plot(sector); summary(sector)
  
  sector = resample(sector, crm_res, method = "bilinear") 
  
  crm_res = as.data.frame(rasterToPoints(crm_res))
  sector = as.data.frame(rasterToPoints(sector))
  
  sector = left_join(crm_res, sector)
  colnames(sector) = c("longitude", "latitude", "cell", "sector")
  summary(sector)
  
  sector %>% 
    ggplot(aes(longitude, latitude, fill = factor(round(sector, 0)))) + 
    geom_tile() + 
    coord_fixed() + 
    ggdark::dark_theme_minimal()
  
  df = merge(df, sector, by = c("cell"))
  
  df$sector = round(df$sector, digits = 0)
  df = df[!is.na(df$sector), ]
  
  colnames(df)[2:3] = c("longitude", "latitude")
  
  # merge reefzones ---------------------------------------------------------
  
  reef$lon = reef$x
  reef$lat = reef$y
  
  reef$reef_zone = as.numeric(as.factor(reef$z))
  reef = as.matrix(reef[,c("lon", "lat", "reef_zone")])
  e = extent(reef[,1:2])
  
  crm_res = rasterFromXYZ(df[,c("longitude", "latitude", "cell")])
  dim(crm_res); crm_res
  
  res = 10  # rasterize it, but be careful with resolutions, lower = better but more missing points
  r <- raster(e, ncol = round((dim(crm_res)[2]/res), digits = 0), nrow = round(dim(crm_res)[1]/res, digits = 0))
  reef <- rasterize(reef[, 1:2], r, reef[,3], fun = mean)
  dim(reef)
  reef
  default_proj = "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
  crs(reef) = default_proj
  plot(reef)
  
  reef = resample(reef, crm_res, method = "bilinear") 
  
  crm_res = as.data.frame(rasterToPoints(crm_res))
  reef = as.data.frame(rasterToPoints(reef))
  
  reef = left_join(crm_res, reef)
  colnames(reef) = c("longitude", "latitude", "cell", "reef")
  summary(reef)
  
  reef %>% 
    ggplot(aes(longitude, latitude, fill = factor(round(reef, 0)))) + 
    geom_tile() + 
    coord_fixed() + 
    ggdark::dark_theme_minimal()
  
  df = merge(df, reef, by = "cell")
  
  df$reef = round(df$reef, digits = 0)
  df = df[!is.na(df$reef), ]
  
  df = df[ , -which(names(df) %in% c("longitude.y", "latitude.y", "longitude.y", "latitude.y"))]
  
  # make strata by depth * sector * reef type -------------------------------
  
  df$depth_bin = ""
  df$depth_bin = ifelse(df$depth <= 0  & df$depth >= -6, 1L, df$depth_bin) 
  df$depth_bin = ifelse(df$depth < -6  & df$depth >= -18, 2L, df$depth_bin) 
  df$depth_bin = ifelse(df$depth < -18, 3L, df$depth_bin) 
  
  df$sector = as.character(df$sector)
  df$reef = as.character(df$reef)
  
  df$strat = paste(df$depth_bin, 
                   df$sector,
                   df$reef,
                   sep = "_")
  
  df$strat = as.numeric(as.factor(df$strat))
  
  df$depth = as.numeric(df$depth*-1)
  
  colnames(df)[2:3] = c("longitude", "latitude")
  
  # utmcoor <- SpatialPoints(cbind(df$longitude, df$latitude), proj4string = CRS("+proj=utm +units=m +zone=55"))
  # longlatcoor <- spTransform(utmcoor,CRS("+proj=longlat"))
  # df$longitude <- coordinates(longlatcoor)[,1]
  # df$latitude <- coordinates(longlatcoor)[,2]
  
  colnames(sector_name) = c("sector_name", "sector")
  colnames(reef_name) = c("reef_name", "reef")
  
  df = merge(df, sector_name, all = T)
  df = merge(df, reef_name, all = T)

  (depth = df %>% 
      ggplot( aes(longitude, latitude, fill = depth)) + 
      # geom_tile(height = 0.001, width = 0.001) +
      geom_raster() + 
      scale_fill_gradientn(colours = colorRamps::matlab.like(100), "Bathymetry (m)") +
      coord_fixed() +
      ggdark::dark_theme_minimal() +
      theme(axis.title = element_blank(),
            legend.position = "right"))
  
  (sector = df %>% 
      ggplot( aes(longitude, latitude, fill = factor(sector))) + 
      geom_raster() +
      scale_fill_discrete("sector") +
      coord_fixed() +
      ggdark::dark_theme_minimal() +
      theme(axis.title = element_blank(),
            legend.position = "right"))
  
  (reef = df %>% 
      ggplot( aes(longitude, latitude, fill = as.factor(reef))) + 
      geom_raster() +
      scale_fill_discrete("reef") +
      coord_fixed() +
      theme_minimal() + 
      ggdark::dark_theme_minimal() +
      theme(axis.title = element_blank(),
            legend.position = "right"))
  
  (strata = df %>% 
      ggplot( aes(longitude, latitude, fill = factor(strat))) + 
      geom_raster() +
      scale_fill_discrete("Strata") +
      coord_fixed() +
      ggdark::dark_theme_minimal() +
      theme(axis.title = element_blank(),
            legend.position = "right"))
  
  df = as.data.frame(df)
  
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
  sp::plot(p)
  
  p <- raster::rasterToPolygons(survey_grid_ncrmp$strat, dissolve = TRUE)
  sp::plot(p)
  
  survey_grid_ncrmp = readAll(survey_grid_ncrmp)
  
  save(survey_grid_ncrmp, file = paste0("data/survey_grid_w_sector_reef/survey_grid_", islands[isl], ".RData"))
  
}

