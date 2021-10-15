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

topo = raster("etopo180_e146_1835_4408.nc")
topo = raster("GEBCO_2020_634c_9bc8_a85c.nc")
topo = raster("usgsCeCrm10_Lon0360_84b5_0453_d127.nc")
topo = raster("srtm15plus_d126_a144_18da.nc")

topo = as.data.frame(rasterToPoints(topo))
names(topo)[3] = "Altitude"
# topo$Altitude = ifelse(topo$Altitude %in% c(-3000:0), topo$Altitude, NA)
topo$Altitude = ifelse(topo$Altitude >= 0, topo$Altitude, NA)
# topo$Altitude = ifelse(topo$Altitude <= 0, topo$Altitude, 0)
topo = topo %>% drop_na() %>% dplyr::select(x, y, Altitude)

topo %>%
  ggplot(aes(x, y, fill = Altitude)) +
  # geom_tile(aes(width = 0.005, height = 0.005)) +
  geom_raster() + 
  scale_fill_viridis_c() +
  coord_fixed() +
  ggdark::dark_theme_minimal() +
  theme(axis.title = element_blank())

# save(topo, file = 'data/gis_bathymetry/raster/Topography_SRTM15.RData')

wireframe(unclass(as.bathy(topo)), 
          shade = T,
          aspect = c(length(unique(topo$y))/length(unique(topo$x)), 0.1),
          par.box = c(col = "transparent"),
          # scales = list(arrows = FALSE, col = "transparent"), # col="black" is required
          # par.settings = list(axis.line = list(col = 'transparent')),
          light.source = c(10,0,10),
          zlab = "", 
          xlab = "",
          ylab = "",
          perspective = T,
          screen = list(z = -15, x = -55),
          zoom = 1.3)

topo = topo %>% 
  ggplot(aes(x=x,y=y,fill=Altitude)) + 
  geom_tile() +
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
