library(raster)
library(rgdal)
library(rgeos)
library(pbapply)

rm(list = ls())

shp_list = list.files(path = "G:/GIS/hardsoft/MHI/", pattern = "shp.shp")

for (shp_i in 1:length(shp_list)) {
  
  # shp_i = 1

  # Import shapefile
  sp_df <- readOGR(paste0("G:/GIS/hardsoft/MHI/", shp_list[shp_i]))[4]

  # first separate hard and soft substrate...
  sp_df_h <- sp_df[sp_df$HardSoft %in% c("Hard"),]
  sp_df_s <- sp_df[sp_df$HardSoft %in% c("Soft"),]
  
  if (length(sp_df_s) == 0) sp_df_s = sp_df_h
  
  # Raster template 
  r <- raster(extent(sp_df))
  projection(r) <- proj4string(sp_df)
  res(r) <- 100 # spatial resolution in m
  
  # Per pixel, identify ID covering largest area
  r_val_h <-  pbsapply(1:ncell(r), function(i) {
    
    r_dupl <- r
    r_dupl[i] <- 1
    p <- rasterToPolygons(r_dupl) # Current cell -> polygon
    sp_df_crp <- crop(sp_df_h, p)   # Crop initial polygons by current cell extent
    
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
  })
  r_val_s <-  pbsapply(1:ncell(r), function(i) {
    
    r_dupl <- r
    r_dupl[i] <- 1
    p <- rasterToPolygons(r_dupl) # Current cell -> polygon
    sp_df_crp <- crop(sp_df_s, p)   # Crop initial polygons by current cell extent
    
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
  })
  
  r_val_h = ifelse(is.na(r_val_h), NA, "1")
  r_val_s = ifelse(is.na(r_val_s), NA, "2")
  
  r_val = data.frame(h = r_val_h, s = r_val_s)
  r_val$hs <- NA
  r_val$hs[r_val$h == 1 & r_val$s == 2] <- "1"
  r_val$hs[r_val$h == 1 & is.na(r_val$s) == T] <- "1"
  r_val$hs[is.na(r_val$h) == T & r_val$s == 2] <- "2"
  
  r_val = r_val$hs
  
  # Write ID values covering the largest area per pixel into raster template
  r[] <- as.numeric(r_val)
  plot(r)
  plot(sp_df, border = "grey45", add = TRUE)
  
  island_name = substr(shp_list[shp_i],1,nchar(shp_list[shp_i])-18)
  
  r = readAll(r)
  
  save(r, file = paste0("data/hardsoft_", island_name, ".RData"))
  
}
