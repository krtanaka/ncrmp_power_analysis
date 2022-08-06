library(ggplot2)
library(dplyr)
library(colorRamps)
library(ggpubr)
library(gridExtra)
library(sp)
library(sf)
library(rgdal)
library(rmapshaper)
library(raster)
require(gridExtra)
library(cowplot)
library(grid)

rm(list = ls())

islands = c("Kauai", #1
            # "Lehua", #2
            "Niihau", #3
            # "Kaula", #4
            "Oahu", #5
            "Molokai", #6
            "Maui", #7
            "Lanai", #8
            # "Molokini", #9
            # "Kahoolawe", #10
            "Hawaii")

load("data/MHI_islands_shp.RData")
crs(ISL_bounds) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
world = ISL_bounds[which(ISL_bounds$ISLAND %in% toupper(islands)),]
world = st_transform(st_as_sf(world))
world <- ms_simplify(world, keep = 0.01, keep_shapes = F)

list = list.files(path = "outputs/", pattern = "MHI_2010-2019"); list

trophic = NULL

for (i in 1:4) {
  
  # i = 1
  
  load(paste0("outputs/", list[i]))
  sp = strsplit(list[i], split = "_")[[1]][3]; sp
  response_scale = strsplit(list[i], split = "_")[[1]][4]; response_scale
  sdm = sdm_output[,c("X", "Y", "year", "est", "zeta_s")]; rm(sdm_output)
  colnames(sdm)[1:2] = c("x", "y")
  sdm$est = exp(sdm$est); hist(sdm$est);summary(sdm$est)
  sdm$sp = sp
  trophic = rbind(trophic, sdm)
  
}

# options(digits = 1)
# options(scipen=10000)

trophic = trophic %>% subset(year >= 2010)

utmcoor <- SpatialPoints(cbind(trophic$x, trophic$y), proj4string = CRS(paste0("+proj=utm +units=km +zone=",4)))
longlatcoor <- spTransform(utmcoor,CRS("+proj=longlat"))
trophic$lon <- coordinates(longlatcoor)[,1]
trophic$lat <- coordinates(longlatcoor)[,2]

MHI_extent = read.csv("data/misc/MHI_Extents.csv")

islands = MHI_extent$ISLAND
islands = islands[! islands %in% c("Kaula", "Lehua", "Molokini")] #remove islands that are too small
islands = islands[! islands %in% c("Kahoolawe")] # remove this island because its missing reef layer

extent = subset(MHI_extent, ISLAND == islands[6])

(PISCIVORE = trophic %>% 
    subset(sp == "PISCIVORE" & 
      lon < extent$LEFT_XMIN &
        lon > extent$RIGHT_XMAX &
        lat > extent$TOP_YMAX & 
        lat < extent$BOTTOM_YMIN) %>% 
    mutate(lon = round(lon, 2),
           lat = round(lat, 2)) %>%
    group_by(lon, lat) %>%
    summarise(est = mean( est)) %>%
    ggplot(aes(lon, lat, fill = est)) + 
    geom_tile() + 
    scale_fill_viridis_c("g_m2") +
    theme_pubr() +
    ggtitle("Piscivore") + 
    theme(legend.position = c(0, 1),
          legend.justification = c(-0.1, 0.9),
          panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray90"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray90"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(face = "bold")))

(PLANKTIVORE = trophic %>% 
    subset(sp == "PLANKTIVORE" & 
             lon < extent$LEFT_XMIN &
             lon > extent$RIGHT_XMAX &
             lat > extent$TOP_YMAX & 
             lat < extent$BOTTOM_YMIN) %>% 
    mutate(lon = round(lon, 2),
           lat = round(lat, 2)) %>%
    group_by(lon, lat) %>%
    summarise(est = mean( est)) %>%
    ggplot(aes(lon, lat, fill = est)) + 
    geom_tile() + 
    scale_fill_viridis_c("g_m2") +
    theme_pubr() +
    ggtitle("Planktivore") + 
    theme(legend.position = c(0, 1),
          legend.justification = c(-0.1, 0.9),
          panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray90"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray90"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(face = "bold")))

(PRIMARY = trophic %>% 
    subset(sp == "PRIMARY" & 
             lon < extent$LEFT_XMIN &
             lon > extent$RIGHT_XMAX &
             lat > extent$TOP_YMAX & 
             lat < extent$BOTTOM_YMIN) %>% 
    mutate(lon = round(lon, 2),
           lat = round(lat, 2)) %>%
    group_by(lon, lat) %>%
    summarise(est = mean( est)) %>%
    ggplot(aes(lon, lat, fill = est)) + 
    geom_tile() + 
    scale_fill_viridis_c("g_m2") +
    theme_pubr() +
    ggtitle("Primary") + 
    theme(legend.position = c(0, 1),
          legend.justification = c(-0.1, 0.9),
          panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray90"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray90"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(face = "bold")))

(SECONDARY = trophic %>% 
    subset(sp == "SECONDARY" & 
             lon < extent$LEFT_XMIN &
             lon > extent$RIGHT_XMAX &
             lat > extent$TOP_YMAX & 
             lat < extent$BOTTOM_YMIN) %>% 
    mutate(lon = round(lon, 2),
           lat = round(lat, 2)) %>%
    group_by(lon, lat) %>%
    summarise(est = mean( est)) %>%
    ggplot(aes(lon, lat, fill = est)) + 
    geom_tile() + 
    scale_fill_viridis_c("g_m2") +
    theme_pubr() +
    ggtitle("Secondary") + 
    theme(legend.position = c(0, 1),
          legend.justification = c(-0.1, 0.9),
          panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray90"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray90"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(face = "bold")))


png(paste0("outputs/sdm_output_", islands[6], ".png"), height = 12, width = 10, res = 500, units = "in")
((PISCIVORE + PLANKTIVORE) / (PRIMARY + SECONDARY))
dev.off()
