library(raster)
library(rgdal)
library(rgeos)
library(pbapply)
library(jubilee)
library(future.apply)
library(dplyr)

rm(list = ls())

### Instruction ###
# 1. Get on pifsc.onaga.gov
# 2. Use multicore and run this script
# 3. Get on Pifsc64, pull & push

# jubilee.mcsapply() = works in Linux and MAC, faster
# future_lapply)() = for Windows, needs plan(multisession) 

# plan(multisession) # Uncomment if you are running this on Windows OS

spatial_resolution = 100 # spatial resolution in m
cores = 64 # number of cores to use

shp_path = "L:/ktanaka/GIS" # pc
# shp_path = "/mnt/ldrive/ktanaka/GIS/" # Onaga
# shp_path = "N:/GIS/Projects/CommonMaps/01_Preprocess/MARI/GUA/"

##################################
### Hard/Soft Bottom Substrate ###
##################################
shp_list = list.files(path = paste0(shp_path, "/hardsoft/"), pattern = "\\.shp$", full.names = T); shp_list
shp_list = shp_list[c(1)]; shp_list

for (shp_i in 1:length(shp_list)) {
  
  start = Sys.time()

  # shp_i = 1
  
  # Import shapefile
  df <- readOGR(shp_list[shp_i])
  df <- df[df$HardSoft != "Land",]
  df <- df[df$HardSoft != "Other",]
  # df <- df[df$HardSoft != "Unknown",]
  
  df@data
  table = data.frame(df@data, i = 0:(length(df)-1)); table
  
  table = table %>% 
    group_by(HardSoft) %>% 
    summarise(i = paste0(i, collapse = ",")) 
  
  # Raster template 
  r <- raster(extent(df))
  projection(r) <- proj4string(df)
  res(r) <- spatial_resolution # spatial resolution in m
  
  # Per pixel, identify ID covering largest area, try jubilee.mcsapply() or pbsapply(), or future_lapply()
  r_val <-  simplify2array(future_lapply(1:ncell(r), function(i) {
  # r_val <-  jubilee.mcsapply(1:ncell(r), mc.cores = cores, function(i) {
    
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
  
  # Write ID values covering the largest area per pixel into raster template
  r[] <- as.numeric(r_val)
  plot(r, col = rainbow(length(unique(r_val))))
  plot(df, border = "grey45", add = TRUE)
  
  # island_name = tolower(substr(shp_list[shp_i], 25, nchar(shp_list[shp_i])-19))
  island_name = tolower(substr(shp_list[shp_i], 25, nchar(shp_list[shp_i])))
  
  r = readAll(r)
  
  raster_and_table = list(r, table)
  
  save(raster_and_table, file = paste0("data/gis_hardsoft/raster/", island_name, "_", spatial_resolution, ".RData"))
  
  end = Sys.time()
  
  time = end - start
  
  print(paste0(island_name, "...done...took ", time, "..."))
  
}

#################
### Reef Zone ###
#################
shp_list = list.files(path = paste0(shp_path, "/reefzone/"), pattern = "\\.shp$", full.names = T); shp_list
shp_list = shp_list[c(1)]; shp_list

for (shp_i in 1:length(shp_list)) {
  
  start = Sys.time()
  
  shp_i = 1 
  
  # Import shapefile
  df <- readOGR(shp_list[shp_i])
  df <- df[df$Zone_Code != "LND",]
  df <- df[df$Zone_Code != "Other",]
  df <- df[df$Zone_Code != "RCF",]
  df <- df[df$Zone_Code != "LAG",]
  df <- df[df$Zone_Code != "UNK",]
  
  df@data
  table = data.frame(df@data, i = 0:(length(df)-1)); table
  
  table = table %>% 
    group_by(Zone_Code) %>% 
    summarise(i = paste0(i, collapse = ",")) 
  
  # Raster template 
  r <- raster(extent(df))
  projection(r) <- proj4string(df)
  res(r) <- spatial_resolution # spatial resolution in m
  
  # Per pixel, identify ID covering largest area, try jubilee.mcsapply() or pbsapply(), or future_lapply()
  # r_val <-  simplify2array(future_lapply(1:ncell(r), function(i) {
  r_val <-  jubilee.mcsapply(1:ncell(r), mc.cores = cores, function(i) {
    
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
  
  # Write ID values covering the largest area per pixel into raster template
  r[] <- as.numeric(r_val)
  plot(r, col = rainbow(length(unique(r_val))))
  plot(df, border = "grey45", add = TRUE)
  
  # island_name = tolower(substr(shp_list[shp_i], 25, nchar(shp_list[shp_i])-19))
  island_name = tolower(substr(shp_list[shp_i], 25, nchar(shp_list[shp_i])))
  
  r = readAll(r)
  
  raster_and_table = list(r, table)
  
  save(raster_and_table, file = paste0("data/gis_hardsoft/raster/", island_name, "_", spatial_resolution, ".RData"))
  
  end = Sys.time()
  
  time = end - start
  
  print(paste0(island_name, "...done...took ", time, "..."))
  
}

##################################
### Regional Sub-Island Sector ###
##################################
shp_list = list.files(path = paste0(shp_path, "/sector/"), pattern = "\\.shp$", full.names = T); shp_list
shp_list = shp_list[c(1)]; shp_list

for (shp_i in 1:length(shp_list)) {
  
  start = Sys.time()
  
  shp_i = 1
  
  # Import shapefile
  df <- readOGR(shp_list[shp_i])
  
  df@data
  table = data.frame(df@data, i = 0:(length(df)-1)); table
  
  table = table %>% 
    # group_by(SEC_NAME) %>% 
    group_by(OID_) %>% 
    summarise(i = paste0(i, collapse = ",")) 
  
  # Raster template 
  r <- raster(extent(df))
  projection(r) <- proj4string(df)
  res(r) <- spatial_resolution # spatial resolution in m
  
  # Per pixel, identify ID covering largest area, try jubilee.mcsapply() or pbsapply(), or future_lapply()
  # r_val <-  simplify2array(future_lapply(1:ncell(r), function(i) {
  r_val <-  jubilee.mcsapply(1:ncell(r), mc.cores = cores, function(i) {
    
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
  
  # Write ID values covering the largest area per pixel into raster template
  r[] <- as.numeric(r_val)
  plot(r, col = topo.colors(length(unique(r))))
  plot(df, border = "grey45", add = TRUE)
  
  # island_name = tolower(substr(shp_list[shp_i], 25, nchar(shp_list[shp_i])-19))
  island_name = tolower(substr(shp_list[shp_i], 23, nchar(shp_list[shp_i])))
  
  r = readAll(r)
  
  raster_and_table = list(r, table)
  
  save(raster_and_table, file = paste0("data/gis_sector/raster/", island_name, "_", spatial_resolution, ".RData"))
  
  end = Sys.time()
  
  time = end - start
  
  print(paste0(island_name, "...done...took ", time, "..."))
  
}

######################
### Marine Reserve ###
######################
shp_list = list.files(path = paste0(shp_path, "/reserve"), pattern = "\\.shp$", full.names = T); shp_list
shp_list = shp_list[c(1)]; shp_list

for (shp_i in 1:length(shp_list)) {
  
  start = Sys.time()
  
  shp_i = 1
  
  # Import shapefile
  df <- readOGR(shp_list[shp_i])
  df <- df[df$Sector != "Land",]
  df <- df[df$Sector != "Other",]

  df@data
  table = data.frame(df@data, i = 0:(length(df)-1)); table
  
  table = table %>% 
    group_by(Sector) %>% 
    summarise(i = paste0(i, collapse = ",")) 
  
  # Raster template 
  r <- raster(extent(df))
  projection(r) <- proj4string(df)
  res(r) <- spatial_resolution # spatial resolution in m
  
  # Per pixel, identify ID covering largest area, try jubilee.mcsapply() or pbsapply(), or future_lapply()
  # r_val <-  simplify2array(future_lapply(1:ncell(r), function(i) {
  r_val <-  jubilee.mcsapply(1:ncell(r), mc.cores = cores, function(i) {
    
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
  
  # Write ID values covering the largest area per pixel into raster template
  r[] <- as.numeric(r_val)
  plot(r, col = topo.colors(length(unique(r))))
  plot(df, border = "grey45", add = TRUE)
  
  # island_name = tolower(substr(shp_list[shp_i], 25, nchar(shp_list[shp_i])-19))
  island_name = tolower(substr(shp_list[shp_i], 24, nchar(shp_list[shp_i])))
  
  r = readAll(r)
  
  raster_and_table = list(r, table)
  
  save(raster_and_table, file = paste0("data/gis_reserve/raster/", island_name, "_", spatial_resolution, ".RData"))
  
  end = Sys.time()
  
  time = end - start
  
  print(paste0(island_name, "...done...took ", time, "..."))
  
}
