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
utmcoor <- SpatialPoints(cbind(bottom_type$X, bottom_type$Y), proj4string = CRS("+proj=utm +units=m +zone=4"))
longlatcoor <- spTransform(utmcoor,CRS("+proj=longlat"))
bottom_type$lon <- coordinates(longlatcoor)[,1]
bottom_type$lat <- coordinates(longlatcoor)[,2]
rm(longlatcoor, utmcoor)
bottom_type = bottom_type %>% filter(!HardSoft %in% c("Unknown", "Land", "Other"))
# bottom_type = bottom_type %>% filter(!HardSoft %in% c("Land"))
bottom_type$hs = ifelse(bottom_type$HardSoft == "Hard", 1, 2)
bottom_type = as.matrix(bottom_type[,c("lon", "lat", "hs")])
e = extent(bottom_type[,1:2])

crm_res = rasterFromXYZ(df[,c("longitude", "latitude", "cell")])
# plot(crm_res, col = rainbow(100))
crm_res %>% 
  rasterToPoints(spatial = T) %>% 
  as.data.frame() %>% 
  ggplot(aes(x, y, fill = cell)) + 
  geom_raster() + 
  coord_fixed() + 
  scale_fill_gradientn(colors = cm.colors(10)) + 
  ggdark::dark_mode()
dim(crm_res)
crm_res

# rasterize it, but be careful with resolutions
res = 10
r <- raster(e, ncol = round((dim(crm_res)[2]/res), digits = 0), nrow = round(dim(crm_res)[1]/res, digits = 0))
bottom_type <- rasterize(bottom_type[, 1:2], r, bottom_type[,3], fun = mean)
bottom_type %>% 
  rasterToPoints(spatial = T) %>% 
  as.data.frame() %>% 
  ggplot(aes(x, y, fill = layer)) + 
  geom_raster() + 
  coord_fixed() + 
  scale_fill_viridis_b("hard_soft") + 
  ggdark::dark_mode()
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

bottom_type %>% 
  ggplot(aes(longitude, latitude, fill = substrate)) + 
  geom_raster(interpolate = T) + 
  coord_fixed() + 
  scale_fill_viridis_b("hard_soft") + 
  ggdark::dark_mode()

df = merge(df, bottom_type, by = "cell")

df$substrate = round(df$substrate, digits = 0)
df = df[!is.na(df$substrate), ]
