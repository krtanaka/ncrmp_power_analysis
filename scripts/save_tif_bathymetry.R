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

# topo = raster("N:/GIS/Projects/CommonMaps/01_Preprocess/MARI/GUA/bathymetry/gua_mb_ld_5m_mosaic.tif")
topo = raster("N:/GIS/Projects/CommonMaps/01_Preprocess/MARI/GUA/bathymetry/gua_nthmp_dem_10m_mosaic.tif")

topo[topo <= -30] <- NA
topo[topo >= 0] <- NA

topo = as.data.frame(rasterToPoints(topo)) %>% drop_na()

topo %>%
  mutate(x = round(x*0.1, 0),
         y = round(y*0.1, 0)) %>%
  group_by(x, y) %>%
  summarise(bathymetry = mean(bathymetry, na.rm = T)) %>%
  ggplot(aes(x, y, fill = bathymetry)) +
  geom_raster() + 
  scale_fill_viridis_c() +
  ggdark::dark_theme_minimal() +
  coord_fixed() +
  theme(axis.title = element_blank())

save(topo, file = 'data/gis_bathymetry/raster/gua_nthmp_dem_10m_mosaic.RData')
