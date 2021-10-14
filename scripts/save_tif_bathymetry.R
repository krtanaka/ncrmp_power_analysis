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
library(marmap)
library(lattice)

topo = raster("L:/ktanaka/GIS/bathymetry/gua_nthmp_dem_10m_mosaic.tif") # Guam

topo <- aggregate(topo, fact = 10)
res(topo)

topo = raster("L:/ktanaka/GIS/bathymetry/sai_mb_5m.tif") # Saipan
topo = raster("L:/ktanaka/GIS/bathymetry/Rota_5m_bathymetry.asc") # Rota
topo = raster("L:/ktanaka/GIS/bathymetry/tinian_5m.asc") # Tinian

#get some sample data
data(meuse.grid)
gridded(meuse.grid) <- ~x+y
meuse.raster <- raster(meuse.grid)
res(meuse.raster)
#[1] 40 40

#aggregate from 40x40 resolution to 120x120 (factor = 3)
meuse.raster.aggregate <- aggregate(meuse.raster, fact=3)
res(meuse.raster.aggregate)
#[1] 120 120

#disaggregate from 40x40 resolution to 10x10 (factor = 4)
meuse.raster.disaggregate <- disaggregate(meuse.raster, fact=4)
res(meuse.raster.disaggregate)
#[1] 10 10

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

wireframe(unclass(as.bathy(topo)), shade = TRUE, aspect = c(1/2, 0.1))


