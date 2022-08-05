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

load("data/crm/Topography_NOAA_CRM_vol10.RData")

df_base = topo; rm(topo)

df_base$longitude = df_base$x
df_base$latitude = df_base$y

# res = 2
# df_base$longitude = round(df_base$x, digits = res)
# df_base$latitude = round(df_base$y, digits = res)

MHI_extent = read.csv("data/misc/MHI_Extents.csv")

islands = MHI_extent$ISLAND
islands = islands[! islands %in% c("Kaula", "Lehua", "Molokini")] #remove islands that are too small
islands = islands[! islands %in% c("Kahoolawe")] # remove this island because its missing reef layer

for (il in 1:length(islands)) {
  
  # il = 7
  island = islands[il]
  extent = subset(MHI_extent, ISLAND == island)
  
  df <- df_base %>% 
    subset(
      longitude < extent$LEFT_XMIN &
        longitude > extent$RIGHT_XMAX &
        latitude > extent$TOP_YMAX & 
        latitude < extent$BOTTOM_YMIN) 
  
  df = df %>%
    group_by(longitude, latitude) %>%
    summarise(Topography = mean(Topography, na.rm = T))
  
  df$cell = 1:dim(df)[1]; df$cell = as.numeric(df$cell)
  df$division = as.numeric(1)
  
  # plot(df$longitude, df$latitude, pch = 20, bty = "l", ann = F, col = 2)
  
  #############################################################
  ### import sector/reefzones shapefile                     ###
  ### adjust resolutions and merge with crm bathymetry data ###
  ### these are outputs from "convert_shp_to_data.frame.R   ###
  #############################################################
  
  if (island == "Hawaii") {
    load("data/gis_sector/raster/HAW.RData"); sector = raster_and_table[[1]]
    load("data/gis_reef/raster/haw.RData"); reef = raster_and_table[[1]]
  }
  if (island == "Kauai") {
    load("data/gis_sector/raster/kau.RData"); sector = raster_and_table[[1]]
    load("data/gis_reef/raster/kau.RData"); reef = raster_and_table[[1]]
  }
  if (island == "Lanai") {
    load("data/gis_sector/raster/lan.RData"); sector = raster_and_table[[1]]
    load("data/gis_reef/raster/lan.RData"); reef = raster_and_table[[1]]
  }
  if (island == "Maui") {
    load("data/gis_sector/raster/mai.RData"); sector = raster_and_table[[1]]
    load("data/gis_reef/raster/mai.RData"); reef = raster_and_table[[1]]
  }
  if (island == "Molokai") {
    load("data/gis_sector/raster/mol.RData"); sector = raster_and_table[[1]]
    load("data/gis_reef/raster/mol.RData"); reef = raster_and_table[[1]]
  }
  if (island == "Niihau") {
    load("data/gis_sector/raster/nii.RData"); sector = raster_and_table[[1]]
    load("data/gis_reef/raster/nii.RData"); reef = raster_and_table[[1]]
  }
  if (island == "Oahu") {
    load("data/gis_sector/raster/oah.RData"); sector = raster_and_table[[1]]
    load("data/gis_reef/raster/oah.RData"); reef = raster_and_table[[1]]
  }
  # if (island == "Kahoolawe") {load("data/gis_sector/data.frame/kah.RData")}
  
  rm(raster_and_table)
  
  sector = rasterToPoints(sector) %>% as.data.frame(); colnames(sector) = c("x", "y", "z")
  reef = rasterToPoints(reef) %>% as.data.frame(); colnames(reef) = c("x", "y", "z")
  
  # plot(sector$lon, sector$lat, pch = ".", bty = "l", ann = F, col = 4)
  # points(reef$lon, reef$lat, pch = ".", bty = "l", ann = F, col = 2)
  
  # merge sectors -----------------------------------------------------------
  utmcoor <- SpatialPoints(cbind(sector$x, sector$y), proj4string = CRS("+proj=utm +units=m +zone=4"))
  longlatcoor <- spTransform(utmcoor,CRS("+proj=longlat"))
  sector$lon <- coordinates(longlatcoor)[,1]
  sector$lat <- coordinates(longlatcoor)[,2]
  rm(longlatcoor, utmcoor)
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
  plot(sector)
  
  sector = resample(sector, crm_res, method = "bilinear") 
  
  crm_res = as.data.frame(rasterToPoints(crm_res))
  sector = as.data.frame(rasterToPoints(sector))
  
  sector = left_join(crm_res, sector)
  colnames(sector) = c("longitude", "latitude", "cell", "sector")
  summary(sector)
  
  sector %>% 
    ggplot(aes(longitude, latitude, color = factor(round(sector, 0)))) + 
    geom_point()
  
  df = merge(df, sector, by = c("cell"))
  
  df$sector = round(df$sector, digits = 0)
  df = df[!is.na(df$sector), ]
  
  colnames(df)[2:3] = c("longitude", "latitude")
  
  # merge reefzones ---------------------------------------------------------
  utmcoor <- SpatialPoints(cbind(reef$x, reef$y), proj4string = CRS("+proj=utm +units=m +zone=4"))
  longlatcoor <- spTransform(utmcoor,CRS("+proj=longlat"))
  reef$lon <- coordinates(longlatcoor)[,1]
  reef$lat <- coordinates(longlatcoor)[,2]
  rm(longlatcoor, utmcoor)
  reef = reef %>% filter(!REEF_ZONE %in% c("Unknown", "Land", "Other", "Reef Crest/Reef Flat"))
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
    ggplot(aes(longitude, latitude, color = factor(round(reef, 0)))) + 
    geom_point()
  
  df = merge(df, reef, by = "cell")
  
  df$reef = round(df$reef, digits = 0)
  df = df[!is.na(df$reef), ]
  
  df = df[ , -which(names(df) %in% c("longitude.y", "latitude.y", "longitude.y", "latitude.y"))]
  
  # make strata by depth * sector * reef type -------------------------------
  
  df$depth_bin = ""
  df$depth_bin = ifelse(df$Topography <= 0  & df$Topography >= -6, 1L, df$depth_bin) 
  df$depth_bin = ifelse(df$Topography < -6  & df$Topography >= -18, 2L, df$depth_bin) 
  df$depth_bin = ifelse(df$Topography < -18, 3L, df$depth_bin) 
  
  df$sector = as.character(df$sector)
  df$reef = as.character(df$reef)
  
  df$strat = paste(df$depth_bin, df$sector, df$reef, sep = "_")
  df$strat = as.numeric(as.factor(df$strat))
  
  df$depth = as.numeric(df$Topography*-1)
  
  colnames(df)[2:3] = c("longitude", "latitude")
  
  (depth = df %>% 
      mutate(longitude = round(longitude, 2),
             latitude  = round(latitude, 2)) %>%
      group_by(longitude, latitude) %>%
      summarise(depth = mean(depth, na.rm = T)) %>%
      ggplot(aes(longitude, latitude, fill = depth)) +
      geom_raster() +
      scale_fill_gradientn(colours = colorRamps::matlab.like(100), "") +
      coord_fixed() +
      theme_pubr() +
      ggtitle("Depth (m)") + 
      theme(legend.position = c(0, 1),
            legend.justification = c(-0.1, 0.9),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank()))
  
  (depth = df %>% 
      mutate(longitude = round(longitude, 2),
             latitude  = round(latitude, 2)) %>% 
      group_by(longitude, latitude) %>% 
      summarise(depth_bin = mean(as.numeric(depth_bin), na.rm = T)) %>% 
      mutate(depth_bin = as.character(round(depth_bin, 0))) %>% 
      ggplot( aes(longitude, latitude, fill = factor(depth_bin))) + 
      geom_raster() +
      scale_fill_brewer(palette = "Blues", "") +
      coord_fixed() +
      theme_pubr() +
      ggtitle("Depth Bins") + 
      theme(legend.position = c(0, 1),
            legend.justification = c(-0.1, 0.9),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank()))
  
  (sector = df %>% 
      mutate(longitude = round(longitude, 2),
             latitude  = round(latitude, 2)) %>% 
      group_by(longitude, latitude) %>% 
      summarise(sector = mean(as.numeric(sector), na.rm = T)) %>% 
      mutate(sector = as.character(round(sector, 0))) %>% 
      ggplot( aes(longitude, latitude, fill = factor(sector))) + 
      geom_raster() +
      scale_fill_brewer(palette = "Oranges", "") +
      coord_fixed() +
      theme_pubr() +
      ggtitle("Island Sector") + 
      theme(legend.position = c(0, 1),
            legend.justification = c(-0.1, 0.9),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank()))
  
  (reef = df %>% 
      mutate(longitude = round(longitude, 2),
             latitude  = round(latitude, 2)) %>% 
      group_by(longitude, latitude) %>% 
      summarise(reef = mean(as.numeric(reef), na.rm = T)) %>% 
      mutate(reef = as.character(round(reef, 0))) %>% 
      ggplot( aes(longitude, latitude, fill = as.factor(reef))) + 
      geom_raster() +
      scale_fill_brewer(palette = "Greens", "", direction = 1) +
      coord_fixed() +
      theme_pubr() +
      ggtitle("Reef Type") + 
      theme(legend.position = c(0, 1),
            legend.justification = c(-0.1, 0.9),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank()))
  
  (strat = df %>% 
      mutate(longitude = round(longitude, 2),
             latitude  = round(latitude, 2)) %>% 
      group_by(longitude, latitude) %>% 
      summarise(strat = mean(as.numeric(strat), na.rm = T)) %>% 
      mutate(strat = as.character(round(strat, 0))) %>%
      ggplot( aes(longitude, latitude, fill = strat)) + 
      geom_raster() +
      # scale_fill_viridis_d("") +
      # scale_fill_viridis_c("") +
      scale_fill_brewer(palette = "Spectral", "") +
      coord_fixed() +
      theme_pubr() +
      ggtitle("Strata") + 
      theme(legend.position = c(0, 1),
            legend.justification = c(-0.1, 0.9),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank()))
  
  
  png(paste0("outputs/survey_grid_", islands[il], ".png"), height = 12, width = 10, res = 500, units = "in")
  ((depth + sector) / (reef + strat))
  dev.off()
  
  df = as.data.frame(df)
  
  cell = rasterFromXYZ(df[,c("longitude", "latitude", "cell")]); plot(cell)
  division = rasterFromXYZ(df[,c("longitude", "latitude", "division")]); plot(division)
  strat = rasterFromXYZ(df[,c("longitude", "latitude", "strat")]); plot(strat)
  depth = rasterFromXYZ(df[,c("longitude", "latitude", "depth")]); plot(depth)
  
  default_proj = "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
  
  ## UTM projection for Hawaii = 5, for Kauai-Maui = 4
  zone <- (floor((df$longitude[1] + 180)/6) %% 60) + 1
  
  ## http://www.gdsihawaii.com/hawpacgis/docs/HawaiiCooSys.pdf
  # utm_proj <- "+proj=utm +ellps=WGS84 +datum=WGS84 +units=km +no_defs"
  utm_proj <- "+proj=utm +zone=4 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"
  # utm_proj <- "+proj=utm +zone=5 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"
  
  crs(cell) = default_proj; cell = projectRaster(cell, crs = utm_proj); plot(cell)
  crs(division) = default_proj; division = projectRaster(division, crs = utm_proj); plot(division)
  crs(strat) = default_proj; strat = projectRaster(strat, crs = utm_proj); plot(strat)
  crs(depth) = default_proj; depth = projectRaster(depth, crs = utm_proj); plot(depth)
  
  survey_grid_kt = stack(cell, division, strat, depth)
  survey_grid_kt$strat = round(survey_grid_kt$strat, digits = 0)
  # values(survey_grid_kt$strat) = ifelse(values(survey_grid_kt$strat) > 3, 3, values(survey_grid_kt$strat))
  values(survey_grid_kt$division) = ifelse(is.na(values(survey_grid_kt$division)), NA, 1)
  
  sp::spplot(survey_grid$cell) #SimSurvey example
  sp::spplot(survey_grid_kt$cell)
  
  sp::spplot(survey_grid$division) #SimSurvey example
  sp::spplot(survey_grid_kt$division)
  
  sp::spplot(survey_grid$strat) #SimSurvey example
  sp::spplot(survey_grid_kt$strat)
  
  sp::spplot(survey_grid$depth) #SimSurvey example
  sp::spplot(survey_grid_kt$depth)
  
  p <- raster::rasterToPolygons(survey_grid$strat, dissolve = TRUE)
  sp::plot(p)
  
  p <- raster::rasterToPolygons(survey_grid_kt$strat, dissolve = TRUE)
  sp::plot(p)
  
  survey_grid_kt = readAll(survey_grid_kt)
  save(survey_grid_kt, file = paste0("data/survey_grid_w_sector_reef/survey_grid_", islands[il], ".RData"))
  
}

