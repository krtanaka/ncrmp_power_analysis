library(raster)
library(rgdal)
library(rgeos)
library(pbapply)
library(jubilee)
library(future.apply)

rm(list = ls())

# jubilee.mcsapply() = works in Linux and MAC, faster
# future_lapply)() = for Windows, needs plan(multisession) 
 
plan(multisession) 
  
####################################
### Hard / Soft Bottom Substrate ###
####################################

shp_list = list.files(path = "G:/GIS/hardsoft/MHI/", pattern = "shp.shp"); shp_list
shp_list = list.files(path = "/mnt/ldrive/ktanaka/GIS/hardsoft/MHI/", pattern = "shp.shp"); shp_list


for (shp_i in 1:length(shp_list)) {
  
  start = Sys.time()
  
  # shp_i = 1
  
  # Import shapefile
  df <- readOGR(paste0("G:/GIS/hardsoft/MHI/", shp_list[shp_i]))[4]
  df@data
  table = data.frame(df@data, i = 0:(length(df)-1)); table
  hard_i = unique(table[table$HardSoft == "Hard",]$i); hard_i
  soft_i = unique(table[table$HardSoft == "Soft",]$i); soft_i
  land_i = unique(table[table$HardSoft == "Land",]$i); land_i
  ukwn_i = unique(table[table$HardSoft == "Unknown",]$i); ukwn_i  
  othr_i = unique(table[table$HardSoft == "Other",]$i); othr_i  
  
  # Raster template 
  r <- raster(extent(df))
  projection(r) <- proj4string(df)
  res(r) <- 5000 # spatial resolution in m
  
  # Per pixel, identify ID covering largest area, try jubilee.mcsapply() or pbsapply(), or future_lapply()
  r_val <-  simplify2array(future_lapply(1:ncell(r), function(i) {
    # r_val <-  jubilee.mcsapply(1:ncell(r), mc.cores = 48, function(i) {
    
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
  }))
  
  r_val = ifelse(r_val %in% hard_i, "Hard", r_val)
  r_val = ifelse(r_val %in% soft_i, "Soft", r_val)
  r_val = ifelse(r_val %in% land_i, "Land", r_val)
  r_val = ifelse(r_val %in% ukwn_i, "Unknown", r_val)
  r_val = ifelse(r_val %in% othr_i, "Other", r_val)
  
  r_val = gsub("Land", NA, r_val)
  r_val = gsub("Unknown", NA, r_val)
  r_val = gsub("Other", NA, r_val)
  
  # Write ID values covering the largest area per pixel into raster template
  r[] <- as.numeric(r_val)
  # plot(r, col = topo.colors(2))
  # plot(df, border = "grey45", add = TRUE)
  
  island_name = substr(shp_list[shp_i],1,nchar(shp_list[shp_i])-18)
  
  r = readAll(r)
  
  save(r, file = paste0("data/hardsoft_", island_name, ".RData"))
  
  end = Sys.time()
  
  time = end - start
  
  print(paste0(island_name, "...done...took ", time, "..."))
  
}


##################################
### regional sub-island sector ###
##################################

shp_list = list.files(path = "G:/GIS/sector/MHI/", pattern = ".shp"); shp_list

plan(multisession) 

for (shp_i in 1:length(shp_list)) {
  
  start = Sys.time()
  
  shp_i = 6
  
  # Import shapefile
  df <- readOGR(paste0("G:/GIS/sector/MHI/", shp_list[shp_i]))[1]
  df@data
  table = data.frame(df@data, i = 0:(length(df)-1)); table
  
  # Raster template 
  r <- raster(extent(df))
  projection(r) <- proj4string(df)
  res(r) <- 5000 # spatial resolution in m
  
  # Per pixel, identify ID covering largest area, try jubilee.mcsapply() or pbsapply(), or future_lapply()
  r_val <-  simplify2array(future_lapply(1:ncell(r), function(i) {
    # r_val <-  jubilee.mcsapply(1:ncell(r), mc.cores = 48, function(i) {
    
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
  }))
  
  # for (t in 1:nrow(table)) {
  #   
  #   # t = 1
  #   r_val = ifelse(r_val %in% table[t,2], table[t,1], r_val)
  #   
  # }
  
  # Write ID values covering the largest area per pixel into raster template
  r[] <- as.numeric(r_val)
  # plot(r, col = topo.colors(length(unique(r))))
  # plot(df, border = "grey45", add = TRUE)
  
  island_name = substr(shp_list[shp_i], 1, nchar(shp_list[shp_i])-12)
  
  r = readAll(r)
  
  save(r, file = paste0("data/island_sector_", island_name, ".RData"))
  
  end = Sys.time()
  
  time = end - start
  
  print(paste0(island_name, "...done...took ", time, "..."))
  
}
