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

topo = raster("N:/GIS/Projects/CommonMaps/01_Preprocess/MARI/GUA/bathymetry/gua_mb_ld_5m_mosaic.tif")
topo = raster("N:/GIS/Projects/CommonMaps/01_Preprocess/MARI/GUA/bathymetry/gua_nthmp_dem_10m_mosaic.tif")

topo[topo <= -30] <- NA
topo[topo >= 0] <- NA

topo = as.data.frame(rasterToPoints(topo)) %>% drop_na()
# topo$Topography = ifelse(topo$gua_mb_ld_5m_mosaic %in% c(-30:0), topo$gua_mb_ld_5m_mosaic, NA)
# topo = topo %>% drop_na() %>% dplyr::select(x, y, Topography)

topo %>%
  group_by(x, y) %>%
  summarise(x = mean(x),
            y = mean(y),
            Topography = mean(gua_nthmp_dem_10m_mosaic, na.rm = T)) %>%
  ggplot(aes(x, y, fill = Topography, color = Topography)) +
  # geom_hex(bins = 50)
  # geom_tile(aes(width = 500, height = 500)) +
  geom_raster() + 
  # geom_point() + 
  scale_fill_viridis_c() +
  coord_fixed() +
  ggdark::dark_theme_minimal() +
  theme(axis.title = element_blank())

save(topo, file = 'data/gis_bathymetry/raster/gua_nthmp_dem_10m_mosaic.RData')
