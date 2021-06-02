library(raster)
library(rgdal)
library(rgeos)
library(pbapply)

rm(list = ls())

shp_list = list.files(path = "G:/GIS/hardsoft/MHI/", pattern = "shp.shp")

for (shp_i in 1:length(shp_list)) {
  
  shp_i = 1
  
  # Import shapefile
  sp_df <- readOGR(paste0("G:/GIS/hardsoft/MHI/", shp_list[shp_i]))[4]

  # combine multiple hard layers...
  sp_df_hard <- sp_df[sp_df$HardSoft %in% c("Hard"),]
  how_many_hard_layers = dim(sp_df_hard)[1]
  sp_df_hard_merged = sp_df_hard[1,]
  
  for (h in 2:how_many_hard_layers) {
    
    # h = 5
    
    # sp_df_h = sp_df_hard[h,]
    
    sp_df_hard_merged = bind(sp_df_hard_merged, sp_df_hard[h,])
    
    print(h)
    
  }
  
  sp_df_hard_merged <- aggregate(sp_df_hard_merged, dissolve = T)
  sp_df_hard_merged@polygons[[1]]@ID = "1"
  pid <- sapply(slot(sp_df_hard_merged, "polygons"), function(x) slot(x, "ID"))
  p.df <- data.frame( ID=1:length(sp_df_hard_merged), row.names = pid)
  sp_df_hard_merged <- SpatialPolygonsDataFrame(sp_df_hard_merged, p.df)
  class(sp_df_hard_merged) 
  plot(sp_df_hard_merged)
  sp_df_hard_merged@data$ID = 1
  
  # combine multiple soft layers...
  sp_df_soft <- sp_df[sp_df$HardSoft %in% c("Soft"),]
  how_many_soft_layers = dim(sp_df_soft)[1]
  sp_df_soft_merged = sp_df_soft[1,]
  
  for (h in 2:how_many_soft_layers) {
    
    # h = 5
    
    sp_df_h = sp_df_soft[h,]
    
    sp_df_soft_merged = rbind(sp_df_soft_merged, sp_df_h)
    print(h)
    
  }
  
  sp_df_soft_merged <- aggregate(sp_df_soft_merged, dissolve = T)
  sp_df_soft_merged@polygons[[1]]@ID = "2"
  pid <- sapply(slot(sp_df_soft_merged, "polygons"), function(x) slot(x, "ID"))
  p.df <- data.frame( ID=1:length(sp_df_soft_merged), row.names = pid)
  sp_df_soft_merged <- SpatialPolygonsDataFrame(sp_df_soft_merged, p.df)
  class(sp_df_soft_merged) 
  plot(sp_df_soft_merged)
  sp_df_soft_merged@data$ID = 2
  
  sp_df = rbind(sp_df_hard_merged, sp_df_soft_merged)

  # hard_or_soft = sp_df$HardSoft
  # hard_or_soft = data.frame(int = c(1:length(hard_or_soft)), hs = hard_or_soft)
  # hard_or_soft$binary = ifelse(hard_or_soft$hs == "Hard", 1, 2)
  # hard_id = hard_or_soft$int[hard_or_soft$binary==1]
  # soft_id = hard_or_soft$int[hard_or_soft$binary==2]
  
  # Raster template 
  r <- raster(extent(sp_df))
  projection(r) <- proj4string(sp_df)
  res(r) <- 10000 # spatial resolution in m
  
  # Per pixel, identify ID covering largest area
  r_val <-  pbsapply(1:ncell(r), function(i) {
    
    r_dupl <- r
    r_dupl[i] <- 1
    p <- rasterToPolygons(r_dupl) # Current cell -> polygon
    sp_df_crp <- crop(sp_df, p)   # Crop initial polygons by current cell extent
    
    # Case 1: no polygon intersecting current cell
    if (is.null(sp_df_crp)) {                   
      
      return(NA)
      
      
      # Case 2: one polygon intersecting current cell  
    } else if (nrow(sp_df_crp@data) < 2) {     
      
      return(rownames(sp_df_crp@data)) 
      
      
      # Case 3: multiple polygons intersecting current cell
    } else {                                 
      
      areas <- gArea(sp_df_crp, byid = TRUE)
      index <- which.max(areas)
      
      return(rownames(sp_df_crp@data)[index])
    }
  }
  )
  
  r_val
  r_val = ifelse(r_val %in% hard_id, 1, r_val)
  r_val = ifelse(r_val %in% soft_id, 2, r_val)
  
  # Write ID values covering the largest area per pixel into raster template
  r[] <- as.numeric(r_val)
  plot(r)
  plot(sp_df, border = "grey45", add = TRUE)
  
  island_name = substr(shp_list[shp_i],1,nchar(shp_list[shp_i])-18)
  
  r = readAll(r)
  
  save(r, file = paste0("data/hardsoft_", island_name, ".Rdata"))
  
}
