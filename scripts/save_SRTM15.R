################################################################################
### Dataset Title: 	SRTM15_PLUS Estimated Topography, 15 seconds, Global, v1 ###
### https://coastwatch.pfeg.noaa.gov/erddap/griddap/srtm15plus.html          ###
################################################################################

library(dplyr)
library(ggplot2)
library(raster)
library(tidyr)
library(marmap)
library(lattice)
library(rayshader)

topo = raster("data/gis_bathymetry/raster/srtm15.nc")
topo = as.data.frame(rasterToPoints(topo))
topo$Altitude = ifelse(topo$Altitude %in% c(-30:0), topo$Altitude, NA)
topo = topo %>% drop_na() %>% dplyr::select(x, y, Altitude)

topo %>%
  ggplot(aes(x, y, fill = Altitude)) +
  # geom_tile(aes(width = 0.005, height = 0.005)) +
  geom_raster() + 
  scale_fill_viridis_c() +
  coord_fixed() +
  ggdark::dark_theme_minimal() +
  theme(axis.title = element_blank())

save(topo, file = 'data/gis_bathymetry/raster/Topography_SRTM15.RData')

wireframe(unclass(as.bathy(topo)), 
          shade = T, 
          aspect = c(1/2, 0.05),
          par.box = c(col = "gray"),
          scales = list(arrows = FALSE, col = "transparent"), # col="black" is required 
          zlab = "", 
          xlab = "",
          ylab = "")
