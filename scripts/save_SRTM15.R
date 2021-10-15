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
topo = raster("data/gis_bathymetry/raster/srtm30.nc")
topo = raster("data/gis_bathymetry/raster/etopo.nc")
topo = raster("data/gis_bathymetry/raster/gebco.nc")

topo = as.data.frame(rasterToPoints(topo))
topo$Altitude = ifelse(topo$Altitude %in% c(-3000:0), topo$Altitude, NA)
topo$Altitude = ifelse(topo$Altitude < 0, topo$Altitude, NA)
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
          aspect = c(1/0.3, 0.1),
          par.box = c(col = "gray"),
          scales = list(arrows = FALSE, col = "transparent"), # col="black" is required 
          zlab = "", 
          xlab = "",
          ylab = "")

topo = ggplot(topo) +
  geom_tile(aes(x=x,y=y,fill=depth)) +
  scale_fill_viridis() + 
  coord_fixed() + 
  theme_void()

plot_gg(topo, 
        multicore = T,
        height = 5,
        width = 6,
        scale = 30, 
        raytrace = TRUE,
        # windowsize = c(1400, 866), 
        zoom = 0.5, 
        phi = 30, 
        theta = 30)

Sys.sleep(0.2)

render_snapshot(clear = TRUE)
