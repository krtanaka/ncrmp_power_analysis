library(raster)
library(rgdal)
library(rasterVis)
library(stringr)

rm(list = ls())

type = c("fish", "benthic")[2]

# Load a SpatialPolygonsDataFrame example shapefile
if (type == "fish") dat <- readOGR('data/eco_zones/Fish_AUTO_1X/Fish_AUTO_1X_sectorshapefile.shp', stringsAsFactors = F)
if (type == "benthic") dat <- readOGR('data/eco_zones/BenthicCover_AUTO_1X/BenthicCover_AUTO_1X_sectorshapefile.shp', stringsAsFactors = F)

dat = dat[1]
class(dat)

zones = dat$ClustID

isl = c('Hawaii', 'Kauai', 'Lanai', 'Maui', 'Molokai', 'Niihau', 'Oahu')

for (i in 1:length(isl)) {
  
  # i = 3
  
  zones_i = str_subset(zones, isl[i], negate = F)
  
  dat_i <- dat[dat$ClustID %in% zones_i,]
  
  # plot(dat_i, pch = ".")
  
  # get names
  nam <- unique(dat_i$ClustID)
  
  # create a data.frame
  nam_df <- data.frame(ID = 1:length(nam), nam = nam)
  
  # Place IDs
  dat_i$ID <- nam_df$ID[match(dat_i$ClustID, nam_df$nam)]
  
  # Define RasterLayer object
  r.raster <- raster()
  
  # Define raster extent
  extent(r.raster) <- extent(dat_i)
  
  # Define pixel size, too coarse then you can't get all IDs
  res(r.raster) <- 0.005
  
  # rasterize
  ras <- rasterize(x = dat_i, y = r.raster, field = "ID")
  
  # ratify raster
  r <- ratify(ras)
  
  # Create levels
  rat <- levels(r)[[1]]
  rat$names <- nam_df$nam
  rat$IDs <- nam_df$ID
  levels(r) <- rat
  
  # levelplot(r)
  plot(r)
  
  tom = r
  
  save(tom, file = paste0("data/eco_zones/", isl[i], "_", type, "_ecozones.RData"))
  
}
