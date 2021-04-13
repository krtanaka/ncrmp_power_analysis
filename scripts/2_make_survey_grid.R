#######################################
### prep survey_grid for simulation ###
#######################################

rm(list = ls())

library(dplyr)
library(ggplot2)
library(raster)
library(sp)
library(rgdal)
library(tidyr)
library(patchwork)
library(SimSurvey)

########################
### see raw GIS data ###
########################

load("data/HAW_Grid.RData")

Hawaii_Survey_Grid %>% 
  group_by(X, Y) %>% 
  summarise(X = mean(X),
            Y = mean(Y),
            depth = mean(DEPTH, na.rm = T)) %>% 
  na_if(-9999) %>% 
  ggplot( aes(X, Y, fill = depth)) + 
  geom_tile(aes(width = 0.005, height = 0.005)) +
  scale_fill_viridis_c("") +
  coord_fixed() +
  ggdark::dark_theme_minimal() + 
  theme(axis.title = element_blank())

#############################################################################
### Topography, NOAA Coastal Relief Model, 3 arc second, Vol. 10 (Hawaii) ###
### https://coastwatch.pfeg.noaa.gov/erddap/griddap/usgsCeCrm10.html      ###
### NOAA NGDC   (Dataset ID: usgsCeCrm10)                                 ###
#############################################################################
# topo = raster("G:/GIS/usgsCeCrm10.nc")
# topo = as.data.frame(rasterToPoints(topo))
# topo$Topography = ifelse(topo$Topography %in% c(-30:0), topo$Topography, NA)
# topo = topo %>% drop_na()
# 
# topo %>%
#   group_by(x, y) %>% 
#   summarise(x = mean(x),
#             y = mean(y),
#             Topography = mean(Topography, na.rm = T)) %>% 
#   ggplot(aes(x, y, fill = Topography)) +
#   geom_tile(aes(width = 0.005, height = 0.005)) +
#   scale_fill_viridis_c() +
#   coord_fixed() +
#   ggdark::dark_theme_minimal() + 
#   theme(axis.title = element_blank())
# 
# save(topo, file = 'data/Topography_NOAA_CRM_vol10.RData')

load("data/Topography_NOAA_CRM_vol10.RData")

df = topo

df$longitude = round(df$x, digits = 2)
df$latitude = round(df$y, digits = 2)

# df$longitude = round(df$x, digits = 3) 
# df$latitude = round(df$y, digits = 3) 

df = df %>% 
  group_by(longitude, latitude) %>% 
  summarise(Topography = mean(Topography, na.rm = T))

df$cell = 1:dim(df)[1]; df$cell = as.numeric(df$cell)
df$division = as.numeric(1)
df$strat = ""
df$strat = ifelse(df$Topography <= 0 & df$Topography >= -6, 1, df$strat)
df$strat = ifelse(df$Topography < -6 & df$Topography >= -18, 2, df$strat)
df$strat = ifelse(df$Topography < -18 & df$Topography >= -30, 3, df$strat)
df$strat = as.numeric(df$strat)
df$depth = as.numeric(df$Topography*-1)

# df <- df %>% subset(longitude < -154.8 & longitude > -156.2 & latitude > 18.8 & latitude < 20.4)
# df <- df %>% subset(longitude < -157.5 & longitude > -158.5 & latitude > 21 & latitude < 22)
 
depth = df %>% 
  ggplot( aes(longitude, latitude, fill = depth)) + 
  geom_tile(aes(width = 0.01, height = 0.01)) +
  scale_fill_viridis_c("") +
  coord_fixed() +
  ggdark::dark_theme_minimal() + 
  theme(axis.title = element_blank())

strat = df %>% 
  ggplot( aes(longitude, latitude, fill = as.factor(strat))) + 
  geom_tile(aes(width = 0.01, height = 0.01)) +
  scale_fill_viridis_d("") +
  coord_fixed() +
  ggdark::dark_theme_minimal() + 
  theme(axis.title = element_blank())

depth + strat

cell = rasterFromXYZ(df[,c("longitude", "latitude", "cell")]); plot(cell)
division = rasterFromXYZ(df[,c("longitude", "latitude", "division")]); plot(division)
strat = rasterFromXYZ(df[,c("longitude", "latitude", "strat")]); plot(strat)
depth = rasterFromXYZ(df[,c("longitude", "latitude", "depth")]); plot(depth)

default_proj = "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

## UTM projection for Hawaii = 5, for Kauai-Maui = 4
## http://www.gdsihawaii.com/hawpacgis/docs/HawaiiCooSys.pdf
# utm_proj <- "+proj=utm +ellps=WGS84 +datum=WGS84 +units=km +no_defs"
utm_proj <- "+proj=utm +zone=4 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"
# utm_proj <- "+proj=utm +zone=5 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"

crs(cell) = default_proj; cell = projectRaster(cell, crs = utm_proj); plot(cell)
crs(division) = default_proj; division = projectRaster(division, crs = utm_proj); plot(division)
crs(strat) = default_proj; strat = projectRaster(strat, crs = utm_proj); plot(strat)
crs(depth) = default_proj; depth = projectRaster(depth, crs = utm_proj); plot(depth)

survey_grid_kt = stack(cell, division, strat, depth)

sp::spplot(survey_grid) #SimSurvey example
sp::spplot(survey_grid_kt)

p <- raster::rasterToPolygons(survey_grid$strat, dissolve = TRUE)
sp::plot(p)

p <- raster::rasterToPolygons(survey_grid_kt$strat, dissolve = TRUE)
sp::plot(p)

survey_grid_kt = readAll(survey_grid_kt)
save(survey_grid_kt, file = "data/survey_grid_kt.RData")
