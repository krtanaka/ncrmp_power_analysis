library(sf)
library(ggplot2)
library(rnaturalearth)
library(marmap)
library(raster)
library(rmapshaper)
library(ggrepel)
library(dplyr)
library(ggOceanMaps)
library(metR)

rm(list = ls())

world <- ne_countries(scale = "large", returnclass = "sf")

b = marmap::getNOAA.bathy(lon1 = min(-161),
                          lon2 = max(-154),
                          lat1 = min(18),
                          lat2 = max(23),
                          resolution = 1)

b = marmap::fortify.bathy(b)

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
world <- ms_simplify(world, keep = 0.001, keep_shapes = F)

load("data/rea/ALL_REA_FISH_RAW_SST.RData")

label = df %>% 
  subset(REGION == "MHI" & ISLAND %in% islands) %>%
  group_by(ISLAND) %>% 
  summarise(lat = mean(LATITUDE),
            lon = mean(LONGITUDE))

# (f1a = basemap(limits = c(-161, -154, 18, 23), 
#                land.col = "gray20", 
#                bathy.style = "poly_greys",
#                bathy.size = 0.5,
#                land.border.col = NA, 
#                bathymetry = TRUE) + 
#     geom_label_repel(data = label,
#                      aes(x = lon, y = lat, label = ISLAND),
#                      fontface = "bold",
#                      nudge_x = c(0.5, 0.5, 0.5, 0.5, 0.5),
#                      nudge_y = c(0.5, 0.5, 0.5, 0.5, 0.5),
#                      label.size = NA,
#                      label.padding=.1,
#                      na.rm=TRUE,
#                      fill = alpha(c("white"),0.8)) +
#     scale_fill_viridis_d("Depth(m)", direction = -1) +
#     scale_x_longitude() +
#     scale_y_latitude() +
#     labs(x = "", y = "") + 
#     theme_pubr(I(20)) + 
#     theme(legend.position = "right"))

(f1a <- ggplot() +
    geom_sf(data = world) +
    coord_sf(crs = st_crs(4135), # old hawaii projection code
             xlim = c(-161, -154),
             ylim = c(18, 23), expand = F) +
    # geom_point(data = df, aes(LONGITUDE, LATITUDE, color = factor(OBS_YEAR)), alpha = 0.5, size = 0.1) +
    # scale_color_manual(values = matlab.like(9), "") + 
    geom_contour(data = b,
                 aes(x = x, y = y, z = z),
                 breaks = seq(-10000, 0, by = 100),
                 size = c(0.1),
                 alpha = 0.8,
                 colour = matlab.like(127637)) +
    geom_label_repel(data = label, 
                     aes(x = lon, y = lat, label = ISLAND), 
                     label.size = NA,
                     fontface = "bold",   
                     label.padding = 0.5, 
                     na.rm = T,
                     fill = alpha(c("white"), 0.8),
                     nudge_x = c(0.5, 0.5, 0.5, 0.5, 0.5),
                     nudge_y = c(0.5, 0.5, 0.5, 0.5, 0.5)) +
    theme_half_open() + 
    theme(panel.background = element_rect(fill = "white"), # bg of the panel
          plot.background = element_rect(fill = "white"), # bg of the plot
          axis.title = element_blank()))

df = df %>%
  subset(REGION == "MHI" & ISLAND %in% islands & OBS_YEAR > 2009) %>%
  group_by(LONGITUDE, LATITUDE, ISLAND, OBS_YEAR) %>%
  summarise(n = n()) %>% 
  group_by(ISLAND, OBS_YEAR) %>%
  summarise(n = n())

(f1b = df %>% 
    ggplot(aes(OBS_YEAR, n, fill = ISLAND, color = ISLAND, group = ISLAND)) + 
    geom_point(size = 5, shape = 21, alpha = 0.7) +
    geom_line(show.legend = F) + 
    scale_fill_manual(values = matlab.like(7), "") +
    scale_color_manual(values = matlab.like(7), "") +
    labs(y = "Sampling effort (n)", x = "Year") + 
    theme_half_open() + 
    scale_x_continuous(breaks = c(2010, 2012, 2013, 2015, 2016, 2019), 
                       labels = c(2010, 2012, 2013, 2015, 2016, 2019)) +
    theme(legend.position = c(0,1),
          legend.justification = c(-0.1, 0.9),
          panel.background = element_rect(fill = alpha('white', 0.5)), # bg of the panel
          plot.background = element_rect(fill = alpha('white', 0.5)), # bg of the plot
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)))

library(cowplot)

f1 <-
  ggdraw() +
  draw_plot(f1a) +
  draw_plot(f1b, x = 0.04, y = 0.14, width = 0.5, height = 0.4)

# Can save the plot with ggsave()
ggsave(filename = "/Users/kisei/Desktop/Fig1.png", 
       plot = f1,
       width = 12, 
       height = 12,
       units = "in",
       dpi = 300)
