library(raster)
library(rgdal)
library(rgeos)
library(pbapply)

# p <- shapefile('N:/GIS/Projects/CommonMaps/01_Preprocess/MHI/HAW/hardsoft/biogeo/haw_hs_biogeo_shp.shp')
# p
# pgeo <- spTransform(p, CRS('+proj=longlat +datum=WGS84'))
# 
# ext <- floor(extent(pgeo))
# rr <- raster(ext, res = 0.1)
# rr <- rasterize(pgeo, rr, field = 1)
# plot(rr)


# Import shapefile
sp_df <- readOGR('N:/GIS/Projects/CommonMaps/01_Preprocess/MHI/HAW/hardsoft/biogeo/haw_hs_biogeo_shp.shp') 
sp_df <- readOGR('N:/GIS/Projects/CommonMaps/01_Preprocess/MHI/OAH/hardsoft/biogeo/oah_hs_biogeo_shp.shp') 
sp_df <- readOGR('N:/GIS/Projects/CommonMaps/01_Preprocess/MHI/KAH/hardsoft/biogeo/kah_hs_biogeo_shp.shp') 

# Raster template 
r <- raster(extent(sp_df))
projection(r) <- proj4string(sp_df)
res(r) <- 1000 # xxx m spatial resolution

# Per pixel, identify ID covering largest area
r_val <-  pbsapply(1:ncell(r), function(i) {
  
  r_dupl <- r
  r_dupl[i] <- 1
  p <- rasterToPolygons(r_dupl) # Current cell -> polygon
  sp_df_crp <- crop(sp_df, p)   # Crop initial polygons by current cell extent
  
  if (is.null(sp_df_crp)) {                   # Case 1: no polygon intersecting current cell
    
    return(NA)
    
  } else if (nrow(sp_df_crp@data) < 2) {     # Case 2: one polygon intersecting current cell  
    
    return(rownames(sp_df_crp@data)) 
    
  } else {                                   # Case 3: multiple polygons intersecting current cell
    
    areas <- gArea(sp_df_crp, byid = TRUE)
    index <- which.max(areas)
    
    return(rownames(sp_df_crp@data)[index])
  }
}
)

# Write ID values covering the largest area per pixel into raster template
r[] <- as.numeric(r_val)
plot(r)
plot(sp_df, border = "grey45", add = TRUE)
