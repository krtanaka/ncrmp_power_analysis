library(raster)
library(rgdal)
library(rgeos)
library(pbapply)
library(jubilee)
library(future.apply)
library(dplyr)

rm(list = ls())

spatial_resolution = 100 # spatial resolution in m

shp_path = "L:/ktanaka/GIS" # pc

##################################
### Hard/Soft Bottom Substrate ###
##################################
shp_list = list.files(path = paste0(shp_path, "/hardsoft/"), pattern = "\\.shp$", full.names = T); shp_list
shp_list = shp_list[c(2, 12:14)]; shp_list

for (shp_i in 1:length(shp_list)) {
  
  start = Sys.time()

  # shp_i = 3

  dat <- shapefile(shp_list[shp_i])
  
  # dat <- spTransform(dat, CRS('+proj=utm +zone=55 +datum=WGS84 +units=km +no_defs'))
  dat <- spTransform(dat, CRS('+proj=utm +zone=55 +datum=WGS84 +units=m +no_defs'))
  # dat <- spTransform(dat, CRS('+proj=longlat +datum=WGS84'))

  dat <- dat[dat$HardSoft %in% c("Hard", "hard", "Unknown", "unknown"),]
  
  # get names
  nam <- unique(dat$HardSoft)
  
  # create a data.frame
  nam_df <- data.frame(ID = 1:length(nam), nam = nam)
  
  # Place IDs
  dat$ID <- nam_df$ID[match(dat$HardSoft,nam_df$nam)]
  
  # Define RasterLayer object
  r.raster <- raster()
  
  # Define raster extent
  extent(r.raster) <- extent(dat)
  
  # Define pixel size
  res(r.raster) <- spatial_resolution
  
  # rasterize
  ras <- rasterize(x = dat, y = r.raster, field = "ID")
  
  # ratify raster
  r <- ratify(ras)
  
  # Create levels
  rat <- levels(r)[[1]]
  rat$names <- nam_df$nam
  rat$IDs <- nam_df$ID
  levels(r) <- rat
  
  # rasterVis::levelplot(r)
  
  island_name = tolower(substr(shp_list[shp_i], 25, nchar(shp_list[shp_i])-18))
  # island_name = tolower(substr(shp_list[shp_i], 36, nchar(shp_list[shp_i])))
  
  r = readAll(r)
  
  save(r, file = paste0("data/gis_hardsoft/raster_alt/", island_name, "_", spatial_resolution, ".RData"))
  
  end = Sys.time()
  
  time = end - start
  
  print(paste0(island_name, "...done...took ", time, "..."))
  
}

#################
### Reef Zone ###
#################
shp_list = list.files(path = paste0(shp_path, "/reefzone/"), pattern = "\\.shp$", full.names = T); shp_list
shp_list = shp_list[c(1, 9:11)]; shp_list
shp_list = shp_list[1]; shp_list

for (shp_i in 1:length(shp_list)) {
  
  start = Sys.time()
  
  # shp_i = 1
  
  # Import shapefile
  df <- readOGR(shp_list[shp_i])
  
  names(df)[2] = "ZONE_CODE"
  
  table(df$ZONE_CODE)
  
  # df <- df[df$ZONE_CODE %in% c("Backreef", "Forereef", "Lagoon", "BRF", "FRF", "LAG"),]

  plot(df, pch = ".")
  
  df@data
  table = data.frame(df@data, i = 0:(length(df)-1)); table
  
  # table = table %>% 
  #   group_by(ZONE_CODE) %>% 
  #   summarise(i = paste0(i, collapse = ",")) 
  
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
  island_name = tolower(substr(shp_list[shp_i], 36, nchar(shp_list[shp_i])))
  
  r = readAll(r)
  
  raster_and_table = list(r, table)
  
  save(raster_and_table, file = paste0("/mnt/ldrive/ktanaka/ncrmp_power_analysis/data/gis_reef/", island_name, "_", spatial_resolution, ".RData"))
  
  end = Sys.time()
  
  time = end - start
  
  print(paste0(island_name, "...done...took ", time, "..."))
  
}

##################################
### Regional Sub-Island Sector ###
##################################
shp_list = list.files(path = paste0(shp_path, "/sector/"), pattern = "\\.shp$", full.names = T); shp_list
shp_list = shp_list[1]; shp_list

for (shp_i in 1:length(shp_list)) {
  
  start = Sys.time()
  
  # shp_i = 1
  
  # Import shapefile
  df <- readOGR(shp_list[shp_i])
  # df <- df[df$Sector != c("Land"),]
  # df <- df[df$Sector != c("Harbor"),]
  
  df@data
  table = data.frame(df@data, i = 0:(length(df)-1)); table
  
  # table = table %>% 
  #   group_by(SEC_NAME) %>%
  #   summarise(i = paste0(i, collapse = ",")) 
  
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
  island_name = tolower(substr(shp_list[shp_i], 34, nchar(shp_list[shp_i])))
  
  r = readAll(r)
  
  raster_and_table = list(r, table)
  
  save(raster_and_table, file = paste0("/mnt/ldrive/ktanaka/ncrmp_power_analysis/data/gis_sector/raster/", island_name, "_", spatial_resolution, ".RData"))
  
  end = Sys.time()
  
  time = end - start
  
  print(paste0(island_name, "...done...took ", time, "..."))
  
}

######################
### Marine Reserve ###
######################
shp_list = list.files(path = paste0(shp_path, "/reserve"), pattern = "\\.shp$", full.names = T); shp_list

for (shp_i in 1:length(shp_list)) {
  
  start = Sys.time()
  
  shp_i = 1
  
  # Import shapefile
  df <- readOGR(shp_list[shp_i])
  # df <- df[df$Sector != "Land",]
  # df <- df[df$Sector != "Other",]

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
  island_name = tolower(substr(shp_list[shp_i], 34, nchar(shp_list[shp_i])))
  
  r = readAll(r)
  
  raster_and_table = list(r, table)
  
  save(raster_and_table, file = paste0("/mnt/ldrive/ktanaka/ncrmp_power_analysis/data/gis_reserve/raster/", island_name, "_", spatial_resolution, ".RData"))
  
  end = Sys.time()
  
  time = end - start
  
  print(paste0(island_name, "...done...took ", time, "..."))
  
}