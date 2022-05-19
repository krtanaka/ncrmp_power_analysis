###################################################
### prep survey_grid for subsequent simulations ###
###################################################

rm(list = ls())

library(raster)
library(rgdal)
library(rasterVis)
library(stringr)

type = c("fish", "benthic")[1]

# Load a SpatialPolygonsDataFrame example shapefile
if (type == "fish") dat <- readOGR('data/gis_eco_zones/fish/Fish_AUTO_1X/Fish_AUTO_1X_sectorshapefile.shp', stringsAsFactors = F)
if (type == "benthic") dat <- readOGR('data/gis_eco_zones/benthic/BenthicCover_AUTO_1X/BenthicCover_AUTO_1X_sectorshapefile.shp', stringsAsFactors = F)

dat = dat[1]
class(dat)

zones = dat$ClustID

load("data/crm/Topography_NOAA_CRM_vol10.RData")

df_base = topo; rm(topo)

df_base$longitude = df_base$x
df_base$latitude = df_base$y

MHI_extent = read.csv("data/misc/MHI_Extents.csv")

islands = MHI_extent$ISLAND
islands = islands[! islands %in% c("Kaula", "Lehua", "Molokini")] # remove islands that are too small
islands = islands[! islands %in% c("Kahoolawe")] %>% as.vector() # remove this island because its missing reef layer

for (i in 1:length(islands)) {
  
  # i = 7
  
  island = islands[i]
  extent = subset(MHI_extent, ISLAND == island)
  
  df <- df_base %>% 
    subset(
      longitude < extent$LEFT_XMIN &
        longitude > extent$RIGHT_XMAX &
        latitude > extent$TOP_YMAX & 
        latitude < extent$BOTTOM_YMIN) 
  
  df = df %>%
    group_by(longitude, latitude) %>%
    summarise(Topography = mean(Topography, na.rm = T))
  
  df$cell = 1:dim(df)[1]; df$cell = as.numeric(df$cell)
  df$division = as.numeric(1)
  crm_res = rasterFromXYZ(df[,c("longitude", "latitude", "cell")])
  
  if (islands[i] == "Hawaii") il = "haw"
  if (islands[i] == "Kauai") il = "kau"
  if (islands[i] == "Lanai") il = "lan"
  if (islands[i] == "Maui") il = "mai"
  if (islands[i] == "Molokai") il = "mol"
  if (islands[i] == "Niihau") il = "nii"
  if (islands[i] == "Oahu") il = "oah"
  
  zones_i = str_subset(zones, islands[i], negate = F)
  
  dat_i <- dat[dat$ClustID %in% zones_i,]
  
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
  
  save(tom, file = paste0("data/gis_eco_zones/", type, "/", il, ".RData"))
  
  default_proj = "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
  crs(tom) = default_proj
  
  tom = resample(tom, crm_res, method = "bilinear") 
  
  crm_res = as.data.frame(rasterToPoints(crm_res))
  tom = as.data.frame(rasterToPoints(tom))
  
  tom = left_join(crm_res, tom) %>% na.omit()
  colnames(tom) = c("longitude", "latitude", "cell", "zone")
  summary(tom)
  
  df = merge(df, tom, by = c("cell"))
  
  df$zone = round(df$zone, digits = 0)
  df = df[!is.na(df$zone), ]
  
  colnames(df)[2:3] = c("longitude", "latitude")
  df$zone = as.character(df$zone)
  
  df$depth_bin = ""
  df$depth_bin = ifelse(df$Topography <= 0  & df$Topography >= -6, 1L, df$depth_bin) 
  df$depth_bin = ifelse(df$Topography < -6  & df$Topography >= -18, 2L, df$depth_bin) 
  df$depth_bin = ifelse(df$Topography < -18, 3L, df$depth_bin) 
  
  df$strat = df$zone # Strata is zone * depth bins
  # df$strat = paste(df$depth_bin, df$zone, sep = "_") # Strata is zone * depth bins
  df$strat = as.numeric(as.factor(df$strat))
  
  df$depth = as.numeric(df$Topography*-1)
  
  df$strat = as.numeric(as.factor(df$strat))
  
  df = as.data.frame(df)
  
  df %>% 
    ggplot(aes(longitude, latitude, color = factor(strat)))+ 
    geom_point() + 
    coord_fixed() + 
    scale_color_viridis_d("") + 
    ggdark::dark_theme_minimal()
  
  cell = rasterFromXYZ(df[,c("longitude", "latitude", "cell")]); plot(cell)
  division = rasterFromXYZ(df[,c("longitude", "latitude", "division")]); plot(division)
  strat = rasterFromXYZ(df[,c("longitude", "latitude", "strat")]); plot(strat)
  depth = rasterFromXYZ(df[,c("longitude", "latitude", "depth")]); plot(depth)
  
  default_proj = "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
  
  ## UTM projection for Hawaii = 5, for Kauai-Maui = 4
  zone <- (floor((df$longitude[1] + 180)/6) %% 60) + 1
  
  ## http://www.gdsihawaii.com/hawpacgis/docs/HawaiiCooSys.pdf
  # utm_proj <- "+proj=utm +ellps=WGS84 +datum=WGS84 +units=km +no_defs"
  utm_proj <- "+proj=utm +zone=4 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"
  # utm_proj <- "+proj=utm +zone=5 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"
  
  crs(cell) = default_proj; cell = projectRaster(cell, crs = utm_proj); plot(cell)
  crs(division) = default_proj; division = projectRaster(division, crs = utm_proj); plot(division)
  crs(strat) = default_proj; strat = projectRaster(strat, crs = utm_proj); plot(strat)
  crs(depth) = default_proj; depth = projectRaster(depth, crs = utm_proj); plot(depth)
  
  survey_grid_kt = stack(cell, division, strat, depth)
  survey_grid_kt$strat = round(survey_grid_kt$strat, digits = 0)
  # values(survey_grid_kt$strat) = ifelse(values(survey_grid_kt$strat) > 3, 3, values(survey_grid_kt$strat))
  values(survey_grid_kt$division) = ifelse(is.na(values(survey_grid_kt$division)), NA, 1)
  
  sp::spplot(survey_grid$cell) #SimSurvey example
  sp::spplot(survey_grid_kt$cell)
  
  sp::spplot(survey_grid$division) #SimSurvey example
  sp::spplot(survey_grid_kt$division)
  
  sp::spplot(survey_grid$strat) #SimSurvey example
  sp::spplot(survey_grid_kt$strat)
  
  sp::spplot(survey_grid$depth) #SimSurvey example
  sp::spplot(survey_grid_kt$depth)
  
  p <- raster::rasterToPolygons(survey_grid$strat, dissolve = TRUE)
  sp::plot(p)
  
  p <- raster::rasterToPolygons(survey_grid_kt$strat, dissolve = TRUE)
  sp::plot(p)
  
  survey_grid_kt = readAll(survey_grid_kt)
  save(survey_grid_kt, file = paste0("data/survey_grid_w_zones/", type, "/survey_grid_", islands[i], ".RData"))
  
}


