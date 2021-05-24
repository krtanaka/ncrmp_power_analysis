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
library(sf)

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

rm(Hawaii_Survey_Grid)

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

df = topo; rm(topo)

df$longitude = df$x
df$latitude = df$y

# res = 2
# df$longitude = round(df$x, digits = res)
# df$latitude = round(df$y, digits = res)

# df <- df %>% subset(longitude < -154.8 & longitude > -156.2 & latitude > 18.8 & latitude < 20.4)
df <- df %>% subset(longitude < -157.5 & longitude > -158.5 & latitude > 21 & latitude < 22)

df = df %>%
  group_by(longitude, latitude) %>%
  summarise(Topography = mean(Topography, na.rm = T))

df$cell = 1:dim(df)[1]; df$cell = as.numeric(df$cell)
df$division = as.numeric(1)

###################################################
### import hard/soft bottom substrate shapefile ###
### adjust resolutions and merge with crm data  ###
###################################################
load("data/oah_hs_biogeo/oah_hs_biogeo_shp.RData")
utmcoor <- SpatialPoints(cbind(bottom_type$X, bottom_type$Y), proj4string = CRS("+proj=utm +zone=4"))
longlatcoor <- spTransform(utmcoor,CRS("+proj=longlat"))
bottom_type$lon <- coordinates(longlatcoor)[,1]
bottom_type$lat <- coordinates(longlatcoor)[,2]
rm(longlatcoor, utmcoor)
# bottom_type = bottom_type %>% filter(!HardSoft %in% c("Unknown", "Land", "Other"))
bottom_type = bottom_type %>% filter(!HardSoft %in% c("Land"))
bottom_type$HS = ifelse(bottom_type$HardSoft == "Hard", 1, 2)
bottom_type = as.matrix(bottom_type[,c("lon", "lat", "HS")])
e = extent(bottom_type[,1:2])

crm_res = rasterFromXYZ(df[,c("longitude", "latitude", "cell")])
plot(crm_res)
dim(crm_res)
crm_res

r <- raster(e, ncol = round((dim(crm_res)[2]/10), digits = 0), nrow = round(dim(crm_res)[1]/10, digits = 0))
bottom_type <- rasterize(bottom_type[, 1:2], r, bottom_type[,3], fun = mean)
plot(bottom_type)
dim(bottom_type)
bottom_type
default_proj = "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
crs(bottom_type) = default_proj
plot(bottom_type)

bottom_type = resample(bottom_type, crm_res, method = "bilinear") 

crm_res = as.data.frame(rasterToPoints(crm_res))
bottom_type = as.data.frame(rasterToPoints(bottom_type))
bottom_type = left_join(crm_res, bottom_type)
colnames(bottom_type) = c("longitude", "latitude", "cell", "substrate")
summary(bottom_type)

qplot(bottom_type$longitude, bottom_type$latitude, color = bottom_type$substrate)

df = merge(df, bottom_type, by = "cell")

df$substrate = round(df$substrate, digits = 0)
df = df[!is.na(df$substrate), ]

# make strata by depth * bottom type
df$strat = ""
df$strat = ifelse(df$Topography <= 0  & df$Topography >= -6  & df$substrate == 1, 1L, df$strat) # shallow & hard
df$strat = ifelse(df$Topography < -6  & df$Topography >= -18 & df$substrate == 1, 2L, df$strat) # mid & hard
df$strat = ifelse(df$Topography < -18 &                        df$substrate == 1, 3L, df$strat) # deep & hard
df$strat = ifelse(df$Topography <= 0  & df$Topography >= -6  & df$substrate == 2, 4L, df$strat) # shallow & soft
df$strat = ifelse(df$Topography < -6  & df$Topography >= -18 & df$substrate == 2, 5L, df$strat) # mid & soft
df$strat = ifelse(df$Topography < -18 &                        df$substrate == 2, 6L, df$strat) # deep & soft
df$strat = as.numeric(df$strat)
df$depth = as.numeric(df$Topography*-1)

colnames(df)[2:3] = c("longitude", "latitude")

depth = df %>% 
  ggplot( aes(longitude, latitude, fill = depth)) + 
  geom_tile(aes(width = 0.005, height = 0.005)) +
  # scale_fill_viridis_c("") +
  scale_fill_gradientn(colours = colorRamps::matlab.like(100), "Bathymetry(m)") +
  coord_fixed() +
  ggdark::dark_theme_minimal() + 
  theme(axis.title = element_blank(),
        legend.position = "bottom")

substrate = df %>% 
  ggplot( aes(longitude, latitude, fill = factor(substrate))) + 
  geom_tile(aes(width = 0.005, height = 0.005)) +
  scale_fill_discrete("Bottom type") +
  # scale_fill_gradientn(colours = colorRamps::matlab.like(100), "Bathymetry(m)") + 
  coord_fixed() +
  ggdark::dark_theme_minimal() + 
  theme(axis.title = element_blank(),
        legend.position = "bottom")

strat = df %>% 
  ggplot( aes(longitude, latitude, fill = as.factor(strat))) + 
  geom_tile(aes(width = 0.005, height = 0.005)) +
  # scale_fill_viridis_d("Strata") +
  # scale_fill_gradientn(colours = colorRamps::matlab.like(6), "Strata") + 
  scale_fill_discrete("Strata") +
  coord_fixed() +
  ggdark::dark_theme_minimal() + 
  theme(axis.title = element_blank(),
        legend.position = "bottom")

png(paste0("/Users/Kisei.Tanaka/Desktop/strata.png"), res = 100, height = 5, width = 10, units = "in")
depth + substrate + strat
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
save(survey_grid_kt, file = "data/survey_grid_kt.RData")
