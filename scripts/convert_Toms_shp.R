library(raster)
library(rgdal)
library(rasterVis)

rm(list = ls())

type = c("fish", "benthic")[2]

# Load a SpatialPolygonsDataFrame example shapefile
if (type == "fish") dat <- readOGR('data/eco_zones/Fish_AUTO_1X/Fish_AUTO_1X_sectorshapefile.shp', stringsAsFactors = F)
if (type == "benthic") dat <- readOGR('data/eco_zones/BenthicCover_AUTO_1X/BenthicCover_AUTO_1X_sectorshapefile.shp', stringsAsFactors = F)

dat = dat[1]
class(dat)

# get names
nam <- unique(dat$ClustID)

# create a data.frame
nam_df <- data.frame(ID = 1:length(nam), nam = nam)

# Place IDs
dat$ID <- nam_df$ID[match(dat$ClustID,nam_df$nam)]

# Define RasterLayer object
r.raster <- raster()

# Define raster extent
extent(r.raster) <- extent(dat)

# Define pixel size, too coarse then you can't get all IDs
res(r.raster) <- 0.01

# rasterize
ras <- rasterize(x = dat, y = r.raster, field = "ID")

# ratify raster
r <- ratify(ras)

# Create levels
rat <- levels(r)[[1]]
rat$names <- nam_df$nam
rat$IDs <- nam_df$ID
levels(r) <- rat

levelplot(r)
plot(r)

tom = r

save(tom, file = paste0("data/eco_zones/", type, "_ecozones.RData"))
