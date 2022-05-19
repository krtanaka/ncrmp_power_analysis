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
library(colorRamps)
library(cowplot)

rm(list = ls())

world <- ne_countries(scale = "small", returnclass = "sf")

a <- st_as_sf(data.frame(plot_id = 1, lat = -156.8, long = 20), 
              coords = c("lat", "long"), crs = 4326)

(globe = ggplot(data = world) +
    geom_sf() +
    geom_sf(data = a, shape = 0, size = 10, color = "red", stroke = 2) +
    coord_sf(crs = "+proj=laea +lat_0=35 +lon_0=-156.8 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs") + 
    theme_bw() +
    theme(panel.background = element_rect(fill = 'white')))#transparent plot bg)

b = marmap::getNOAA.bathy(lon1 = -160.5,
                          lon2 = -153,
                          lat1 = 17,
                          lat2 = 23,
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
world <- ms_simplify(world, keep = 0.1, keep_shapes = F)

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
             xlim = c(-160.5, -154.7),
             ylim = c(18.8, 22.5), expand = F) +
    # geom_point(data = df, aes(LONGITUDE, LATITUDE, color = factor(OBS_YEAR)), alpha = 0.5, size = 0.1) +
    # scale_color_manual(values = matlab.like(9), "") + 
    # geom_contour(data = b,
    #              aes(x = x, y = y, z = z),
    #              breaks = seq(-10000, 0, by = 500),
    #              size = c(0.2),
    #              alpha = 0.8,
    #              colour = matlab.like(29771)) +
    geom_contour(data = b,
                 aes(x = x, y = y, z = z, colour = stat(level)),
                 breaks = seq(-1000, 0, by = 100),
                 size = c(0.5)) +
    # scale_colour_distiller(palette = "RdYlBu", direction = -1, "Depth (m)") +
    scale_colour_viridis_c("Depth (m)") +
    geom_label_repel(data = label, 
                     aes(x = lon, y = lat, label = ISLAND), 
                     label.size = NA,
                     fontface = "bold",   
                     label.padding = 0.5, 
                     na.rm = T,
                     fill = alpha(c("white"), 0.9),
                     nudge_x = c(0.5, 0.5, 0.5, 0.5, 0.5),
                     nudge_y = c(0.5, 0.5, 0.5, 0.5, 0.5)) +
    theme_linedraw() + 
    theme(legend.position = c(0.84,1),
          legend.justification = c(1.1,1.1),
          panel.background = element_rect(fill = "white"), # bg of the panel
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
    scale_fill_viridis_d("") + 
    scale_color_viridis_d("") + 
    labs(y = "Sampling effort (n)", x = "") + 
    theme_bw() +
    scale_x_continuous(breaks = c(2010, 2012, 2013, 2015, 2016, 2019), 
                       labels = c(2010, 2012, 2013, 2015, 2016, 2019)) +
    theme(legend.position = c(0,1),
          legend.justification = c(-0.1, 0.9),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = alpha('white', 0.8)), # bg of the panel
          plot.background = element_rect(fill = alpha('white', 0.8)), # bg of the plot
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)))

library(cowplot)

f1 <-
  ggdraw() +
  draw_plot(f1a) +
  draw_plot(f1b, x = 0.045, y = 0.124, width = 0.5, height = 0.4) + 
  draw_plot(globe, x = 0.816, y = 0.705, width = 0.2, height = 0.2)

ggsave(filename = "outputs/fig1.png", 
       plot = f1,
       width = 12, 
       height = 10,
       units = "in",
       dpi = 500)
