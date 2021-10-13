#############################################################################
### Topography, NOAA Coastal Relief Model, 3 arc second, Vol. 10 (Hawaii) ###
### https://coastwatch.pfeg.noaa.gov/erddap/griddap/usgsCeCrm10.html      ###
### NOAA NGDC   (Dataset ID: usgsCeCrm10)                                 ###
#############################################################################

rm(list = ls())

library(dplyr)
library(ggplot2)
library(raster)
library(tidyr)

topo = raster("L:/ktanaka/GIS/bathymetry/gua_nthmp_dem_10m_mosaic.tif") # Guam
topo = raster("L:/ktanaka/GIS/bathymetry/sai_mb_5m.tif") #Saipan
topo = raster("L:/ktanaka/GIS/bathymetry/Rota_5m_bathymetry.asc") 
topo = raster("L:/ktanaka/GIS/bathymetry/tinian_5m.asc") 

topo[topo <= -30] <- NA
topo[topo >= 0] <- NA

topo = as.data.frame(rasterToPoints(topo)) %>% drop_na()
colnames(topo) = c("x", "y", "depth")

topo %>%
  mutate(x = round(x*0.1, 0),
         y = round(y*0.1, 0)) %>%
  group_by(x, y) %>%
  summarise(depth = mean(depth, na.rm = T)) %>%
  ggplot(aes(x, y, fill = depth)) +
  geom_raster() + 
  scale_fill_viridis_c() +
  ggdark::dark_theme_minimal() +
  coord_fixed() +
  theme(axis.title = element_blank())

save(topo, file = 'data/gis_bathymetry/raster/gua_nthmp_dem_10m_mosaic.tif.RData')
save(topo, file = 'data/gis_bathymetry/raster/sai_mb_5m.tif.RData')
save(topo, file = 'data/gis_bathymetry/raster/rota_5m_bathymetry.asc.RData')
save(topo, file = 'data/gis_bathymetry/raster/tinian_5m.asc.RData')


