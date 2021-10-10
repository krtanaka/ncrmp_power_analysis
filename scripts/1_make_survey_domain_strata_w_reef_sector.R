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

load("data/gis_bathymetry/raster/gua_nthmp_dem_10m_mosaic.RData")

df = topo; rm(topo)

# change to 100 m res
df$longitude = round(df$x*0.1, 0)
df$latitude = round(df$y*0.1, 0)

df = df %>% 
  group_by(longitude, latitude) %>% 
  summarise(gua_nthmp_dem_10m_mosaic = mean(gua_nthmp_dem_10m_mosaic))

df$cell = 1:dim(df)[1]; df$cell = as.numeric(df$cell)
df$division = as.numeric(1)

#############################################################
### import sector/reefzones shapefile                     ###
### adjust resolutions and merge with crm bathymetry data ###
### these are outputs from "convert_shp_to_data.frame.R   ###
#############################################################

load("data/gis_sector/gua_base_land_openwater_mpa__100.RData"); sector = raster_and_table[[1]]
# load("data/gis_reef/raster/haw.RData"); reef = raster_and_table[[1]]

rm(raster_and_table)

sector = rasterToPoints(sector) %>% as.data.frame(); colnames(sector) = c("x", "y", "z")
reef = rasterToPoints(reef) %>% as.data.frame(); colnames(reef) = c("x", "y", "z")

# merge sectors -----------------------------------------------------------
sector$lon = sector$x
sector$lat = sector$y

sector$lon = round(sector$x*0.1, 0)
sector$lat = round(sector$y*0.1, 0)
sector = sector %>% 
  group_by(lon, lat) %>% 
  summarise(z = round(mean(z), 0))

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
  mutate(longitude = round(longitude*0.01, 0),
         latitude = round(latitude*0.01, 0)) %>%
  group_by(longitude, latitude) %>%
  summarise(sector = mean(sector, na.rm = T)) %>%
  ggplot(aes(longitude, latitude, fill = factor(round(sector, 0)))) + 
  geom_tile() + 
  coord_fixed() + 
  ggdark::dark_theme_minimal()

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
# reef = reef %>% filter(!REEF_ZONE %in% c("Unknown", "Land", "Other", "Reef Crest/Reef Flat"))
reef$reef_zone = as.numeric(as.factor(reef$z))
reef = as.matrix(reef[,c("lon", "lat", "reef_zone")])
e = extent(reef[,1:2])

crm_res = rasterFromXYZ(df[,c("longitude", "latitude", "cell")])
dim(crm_res); crm_res

res = 10  # rasterize it, but be careful with resolutions, lower = better but more missing points
r <- raster(e, ncol = round((dim(crm_res)[2]/res), digits = 0), nrow = round(dim(crm_res)[1]/res, digits = 0))
reef <- rasterize(reef[, 1:2], r, reef[,3], fun = mean)
reef %>% 
  rasterToPoints(spatial = T) %>% 
  as.data.frame() %>%
  ggplot(aes(x, y, fill = layer)) + 
  geom_raster() + 
  coord_fixed() +
  # scale_fill_viridis_b("") + 
  ggdark::dark_theme_minimal()
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
  geom_point() + 
  coord_fixed() + 
  # scale_fill_viridis_b("") + 
  ggdark::dark_theme_minimal()

df = merge(df, reef, by = "cell")

df$reef = round(df$reef, digits = 0)
df = df[!is.na(df$reef), ]

df = df[ , -which(names(df) %in% c("longitude.y", "latitude.y", "longitude.y", "latitude.y"))]

# make strata by depth * sector * reef type -------------------------------

df$depth_bin = ""
df$depth_bin = ifelse(df$gua_nthmp_dem_10m_mosaic <= 0  & df$gua_nthmp_dem_10m_mosaic >= -6, 1L, df$depth_bin) 
df$depth_bin = ifelse(df$gua_nthmp_dem_10m_mosaic < -6  & df$gua_nthmp_dem_10m_mosaic >= -18, 2L, df$depth_bin) 
df$depth_bin = ifelse(df$gua_nthmp_dem_10m_mosaic < -18, 3L, df$depth_bin) 

df$sector = as.character(df$sector)
df$reef = as.character(df$reef)

df$strat = paste(df$depth_bin, df$sector
                 # , df$reef
                 , sep = "_")
df$strat = as.numeric(as.factor(df$strat))

df$depth = as.numeric(df$gua_nthmp_dem_10m_mosaic*-1)

colnames(df)[2:3] = c("longitude", "latitude")

# utmcoor <- SpatialPoints(cbind(df$longitude, df$latitude), proj4string = CRS("+proj=utm +units=m +zone=55"))
# longlatcoor <- spTransform(utmcoor,CRS("+proj=longlat"))
# df$longitude <- coordinates(longlatcoor)[,1]
# df$latitude <- coordinates(longlatcoor)[,2]

(depth = df %>% 
    ggplot( aes(longitude.y, latitude.y, fill = depth)) + 
    geom_raster(interpolate = T) +
    scale_fill_gradientn(colours = colorRamps::matlab.like(100), "Bathymetry (m)") +
    coord_fixed() +
    ggdark::dark_theme_minimal() +
    theme(axis.title = element_blank(),
          legend.position = "bottom"))

(sector = df %>% 
    ggplot( aes(longitude.y, latitude.y, fill = factor(sector))) + 
    geom_raster() +
    scale_fill_discrete("sector") +
    coord_fixed() +
    ggdark::dark_theme_minimal() +
    theme(axis.title = element_blank(),
          legend.position = "bottom"))

reef = df %>% 
  ggplot( aes(longitude, latitude, fill = as.factor(reef))) + 
  geom_tile(aes(width = 0.005, height = 0.005)) +
  scale_fill_discrete("reef") +
  coord_fixed() +
  theme_minimal() + 
  ggdark::dark_theme_minimal() +
  theme(axis.title = element_blank(),
        legend.position = "bottom")

# pdf(paste0("outputs/survey_grid_", islands[il], ".pdf"), height = 8, width = 10)
p =  depth + sector + reef
print(p)
# dev.off()

df = as.data.frame(df)

cell = rasterFromXYZ(df[,c("longitude.y", "latitude.y", "cell")]); plot(cell)
division = rasterFromXYZ(df[,c("longitude.y", "latitude.y", "division")]); plot(division)
strat = rasterFromXYZ(df[,c("longitude.y", "latitude.y", "strat")]); plot(strat)
depth = rasterFromXYZ(df[,c("longitude.y", "latitude.y", "depth")]); plot(depth)

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

