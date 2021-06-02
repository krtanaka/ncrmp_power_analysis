library(raster)
library(rgdal)
library(rgeos)
library(pbapply)

rm(list = ls())

shp_list = list.files(path = "G:/GIS/hardsoft/MHI/", pattern = "shp.shp"); shp_list

for (shp_i in 1:length(shp_list)) {
  
  # shp_i = 8
  
  # Import shapefile
  df <- readOGR(paste0("G:/GIS/hardsoft/MHI/", shp_list[shp_i]))[4]
  df@data
  table = data.frame(df@data, i = 1:length(df))
  hard_i = unique(table[table$HardSoft == "Hard",]$i)
  soft_i = unique(table[table$HardSoft == "Soft",]$i)
  land_i = unique(table[table$HardSoft == "Land",]$i)
  ukwn_i = unique(table[table$HardSoft == "Unknown",]$i)
  
  # Raster template 
  r <- raster(extent(df))
  projection(r) <- proj4string(df)
  res(r) <- 500 # spatial resolution in m
  
  # Per pixel, identify ID covering largest area
  r_val <-  pbsapply(1:ncell(r), function(i) {
    
    r_dupl <- r
    r_dupl[i] <- 1
    p <- rasterToPolygons(r_dupl) # Current cell -> polygon
    sp_df_crp <- crop(df, p)   # Crop initial polygons by current cell extent
    
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
  
  r_val = ifelse(r_val %in% hard_i, "1", r_val)
  r_val = ifelse(r_val %in% soft_i, "2", r_val)
  r_val = ifelse(r_val %in% land_i, "3", r_val)
  r_val = ifelse(r_val %in% ukwn_i, "4", r_val)
  
  # Write ID values covering the largest area per pixel into raster template
  r[] <- as.numeric(r_val)
  plot(r)
  plot(df, border = "grey45", add = TRUE)
  
  island_name = substr(shp_list[shp_i],1,nchar(shp_list[shp_i])-18)
  
  r = readAll(r)
  
  save(r, file = paste0("data/hardsoft_", island_name, ".RData"))
  
}
