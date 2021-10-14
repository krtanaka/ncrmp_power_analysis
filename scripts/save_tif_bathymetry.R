rm(list = ls())

library(dplyr)
library(ggplot2)
library(raster)
library(tidyr)
library(marmap)
library(lattice)

islands = c("gua", "rot", "sai", "tin")

for (isl in 1:length(islands)) {
  
  isl = 4
  
  if (islands[isl] == "gua")  topo = raster("L:/ktanaka/GIS/bathymetry/gua_nthmp_dem_10m_mosaic.tif") # Guam
  if (islands[isl] == "rot") topo = raster("L:/ktanaka/GIS/bathymetry/Rota_5m_bathymetry.asc") # Rota
  if (islands[isl] == "sai") topo = raster("L:/ktanaka/GIS/bathymetry/sai_mb_5m.tif") # Saipan
  if (islands[isl] == "tin") topo = raster("L:/ktanaka/GIS/bathymetry/tinian_5m.asc") # Tinian
  
  topo[topo <= -30] <- NA
  topo[topo >= 0] <- NA
  
  topo <- aggregate(topo, fact = 100/res(topo))
  res(topo)
  
  topo = as.data.frame(rasterToPoints(topo)) %>% drop_na()
  colnames(topo) = c("x", "y", "depth")
  
  topo %>%
    ggplot(aes(x, y, fill = depth)) +
    geom_raster() +
    scale_fill_viridis_c() +
    ggdark::dark_theme_minimal() +
    coord_fixed() +
    theme(axis.title = element_blank())
  
  save(topo, file = paste0('data/gis_bathymetry/raster/', islands[isl], '.RData'))
  
}

wireframe(unclass(as.bathy(topo)), shade = T, aspect = c(1/2, 0.1))
