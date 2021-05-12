rm(list = ls())

library(readr)
library(dplyr)
library(ggplot2)
library(raster)

# df <- read_csv("data/jplmursst_strm15.csv")
# df$bathymetry_m = ifelse(df$bathymetry_m %in% c(-30:0), df$bathymetry_m, NA)
# # df$bathymetry_m = ifelse(df$bathymetry_m > 0, NA, df$bathymetry_m)
# 
# df = df %>% drop_na()
# 
topo = raster("/Users/Kisei.Tanaka/Desktop/usgsCeCrm10.nc")
topo = as.data.frame(rasterToPoints(topo))
topo$Topography = ifelse(topo$Topography %in% c(-30:0), topo$Topography, NA)
topo = topo %>% drop_na()
save(topo, file = 'data/Topography_NOAA_CRM_vol10.RData')

topo %>%
  ggplot(aes(x, y, fill = Topography, color = Topography)) +
  geom_tile() +
  scale_fill_viridis_c() +
  scale_color_viridis_c() +
  coord_fixed() +
  ggdark::dark_theme_void()

load("data/Topography_NOAA_CRM_vol10.RData"); df = topo

df$cell = 1:dim(df)[1]
df$division = 1
df$strat = 2
df$strat = ifelse(df$Topography %in% c(-10:0), 1, df$strat)
df$strat = ifelse(df$Topography %in% c(-30:-20), 3, df$strat)
df$depth = df$Topography*-1
df$longitude = df$x
df$latitude = df$y

df %>% 
  ggplot( aes(longitude, latitude, fill = depth)) + 
  geom_tile() +
  scale_fill_viridis_c() + 
  coord_fixed() + 
  ggdark::dark_theme_void()

cell = rasterFromXYZ(df[,c("longitude", "latitude", "cell")]); plot(cell)
division = rasterFromXYZ(df[,c("longitude", "latitude", "division")]); plot(division)
strat = rasterFromXYZ(df[,c("longitude", "latitude", "strat")]); plot(strat)
depth = rasterFromXYZ(df[,c("longitude", "latitude", "depth")]); plot(depth)

survey_grid_kt = stack(cell, division, strat, depth)

sp::spplot(survey_grid)
sp::spplot(survey_grid_kt)

p <- raster::rasterToPolygons(survey_grid$strat, dissolve = TRUE)
sp::plot(p)

p <- raster::rasterToPolygons(survey_grid_kt$strat, dissolve = TRUE)
sp::plot(p)
