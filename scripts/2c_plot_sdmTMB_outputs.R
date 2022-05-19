library(ggplot2)
library(dplyr)
library(colorRamps)
library(ggpubr)
library(gridExtra)
library(sp)
library(sf)
library(rgdal)
library(rmapshaper)

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

options(digits = 1)

trophic = trophic %>% subset(year >= 2010)

utmcoor <- SpatialPoints(cbind(trophic$x, trophic$y), proj4string = CRS(paste0("+proj=utm +units=km +zone=",4)))
longlatcoor <- spTransform(utmcoor,CRS("+proj=longlat"))
trophic$lon <- coordinates(longlatcoor)[,1]
trophic$lat <- coordinates(longlatcoor)[,2]

df = trophic %>% 
  subset(sp == "PISCIVORE") %>%
  mutate(lon = round(lon, 1),
         lat = round(lat, 1)) %>%
  group_by(lon, lat) %>% 
  summarise(est = mean(zeta_s)) %>% 
  mutate(abs_est = abs(est))

(c1 = ggplot() +
    geom_sf(data = world) +
    geom_point(data = df,  aes_string("lon", "lat", color = "est"), size = scale(df$abs_est)) +
    coord_sf(crs = st_crs(4135)) +
    scale_color_distiller(palette ="RdBu",
                          direction = -1,
                          # limits = c(quantile(df$est, 0.999)*-1, quantile(df$est, 0.999)) +
                          "") + 
    ggtitle("(a) Linear trend 2010-2019") + 
    theme_classic() +
    annotate("text",  x = Inf, y = Inf, label = "Piscivore", vjust = 1.2, hjust = 1.2, size = 4) +  
    theme(legend.position = c(0, 0), 
          legend.justification = c(-0.1, -0.1),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.box.background = element_rect(fill = "transparent", colour = "transparent"),
          # legend.text = element_text(color = "white", size = 12),
          legend.key.size = unit(0.5, "cm"),
          # panel.background = element_rect(fill = "gray10", colour = "gray10"),
          panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray80"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray80"),
          # axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(face = "bold")))

df = trophic %>% 
  subset(sp == "PLANKTIVORE") %>%
  mutate(lon = round(lon, 1),
         lat = round(lat, 1)) %>%
  group_by(lon, lat) %>% 
  summarise(est = mean(zeta_s)) %>% 
  mutate(abs_est = abs(est))

(c2 = ggplot() +
    geom_sf(data = world) +
    geom_point(data = df,  aes_string("lon", "lat", color = "est"), size = scale(df$abs_est)) +
    coord_sf(crs = st_crs(4135)) +
    scale_color_distiller(palette ="RdBu",
                          direction = -1,
                          # limits = c(quantile(df$est, 0.999)*-1, quantile(df$est, 0.999)) +
                          "") + 
    ggtitle("") + 
    theme_classic() +
    annotate("text",  x = Inf, y = Inf, label = "Piscivore", vjust = 1.2, hjust = 1.2, size = 4) + 
    theme(legend.position = c(0, 0), 
          legend.justification = c(-0.1, -0.1),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.box.background = element_rect(fill = "transparent", colour = "transparent"),
          # legend.text = element_text(color = "white", size = 12),
          legend.key.size = unit(0.5, "cm"),
          # panel.background = element_rect(fill = "gray10", colour = "gray10"),
          panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray80"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray80"),
          # axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(face = "bold")))

df = trophic %>% 
  subset(sp == "PRIMARY") %>%
  mutate(lon = round(lon, 1),
         lat = round(lat, 1)) %>%
  group_by(lon, lat) %>% 
  summarise(est = mean(zeta_s)) %>% 
  mutate(abs_est = abs(est))

(c3 = ggplot() +
    geom_sf(data = world) +
    geom_point(data = df,  aes_string("lon", "lat", color = "est"), size = scale(df$abs_est)) +
    coord_sf(crs = st_crs(4135)) +
    scale_color_distiller(palette ="RdBu",
                          direction = -1,
                          # limits = c(quantile(df$est, 0.999)*-1, quantile(df$est, 0.999)) +
                          "") + 
    ggtitle(" ") + 
    theme_classic() +
    annotate("text",  x = Inf, y = Inf, label = "Primary", vjust = 1.2, hjust = 1.2, size = 4) +  
    theme(legend.position = c(0, 0), 
          legend.justification = c(-0.1, -0.1),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.box.background = element_rect(fill = "transparent", colour = "transparent"),
          # legend.text = element_text(color = "white", size = 12),
          legend.key.size = unit(0.5, "cm"),
          # panel.background = element_rect(fill = "gray10", colour = "gray10"),
          panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray80"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray80"),
          # axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(face = "bold")))

df = trophic %>%
  subset(sp == "SECONDARY") %>%
  mutate(lon = round(lon, 1),
         lat = round(lat, 1)) %>%
  group_by(lon, lat) %>% 
  summarise(est = mean(zeta_s)) %>% 
  mutate(abs_est = abs(est))

(c4 = ggplot() +
    geom_sf(data = world) +
    geom_point(data = df,  aes_string("lon", "lat", color = "est"), size = scale(df$abs_est)) +
    coord_sf(crs = st_crs(4135)) +
    scale_color_distiller(palette ="RdBu",
                          direction = -1,
                          # limits = c(quantile(df$est, 0.999)*-1, quantile(df$est, 0.999)) +
                          "") + 
    ggtitle("  ") + 
    theme_classic() +
    annotate("text",  x = Inf, y = Inf, label = "Secondary", vjust = 1.2, hjust = 1.2, size = 4) + 
    theme(legend.position = c(0, 0), 
          legend.justification = c(-0.1, -0.1),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.box.background = element_rect(fill = "transparent", colour = "transparent"),
          # legend.text = element_text(color = "white", size = 12),
          legend.key.size = unit(0.5, "cm"),
          # panel.background = element_rect(fill = "gray10", colour = "gray10"),
          panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray80"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray80"),
          # axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(face = "bold")))

png('outputs/fig2a.png', height = 3, width = 17, units = "in", res = 500)
grid.arrange(c1, c2, c3, c4, nrow = 1)
dev.off()

(c_supp = trophic %>% 
    group_by(x, y, sp) %>% 
    summarise(est = mean(zeta_s)) %>% 
    ggplot(aes(x, y, fill = est)) + 
    geom_tile(height = 0.8, width = 0.8) +
    coord_fixed() + 
    facet_wrap(~sp, nrow = 2) +
    scale_fill_gradient2("Linear trend") +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.background = element_rect(fill = "gray10", colour = "gray10"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "gray20"), 
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "gray20")))

# png('outputs/fig2_supp.png', height = 7, width = 10, units = "in", res = 500)
# (c_supp)
# dev.off()

df = trophic %>% 
  subset(sp == "PISCIVORE") %>%
  mutate(lon = round(lon, 1),
         lat = round(lat, 1)) %>%
  group_by(lon, lat) %>% 
  summarise(est = median(est)) %>%  
  mutate(abs_est = abs(est))

(m1 = ggplot() +
    geom_sf(data = world) +
    geom_point(data = df,  aes_string("lon", "lat", color = "est"), size = scale(df$abs_est)) +
    coord_sf(crs = st_crs(4135)) +
    scale_color_gradientn(colours = matlab.like(100), "") + 
    scale_color_viridis_c("") + 
    ggtitle("(b) Median biomass 2010-2019") + 
    theme_classic() +
    annotate("text",  x = Inf, y = Inf, label = "Piscivore", vjust = 1.2, hjust = 1.2, size = 4) + 
    theme(legend.position = c(0, 0), 
          legend.justification = c(-0.1, -0.1),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.box.background = element_rect(fill = "transparent", colour = "transparent"),
          # legend.text = element_text(color = "white", size = 12),
          legend.key.size = unit(0.5, "cm"),
          # panel.background = element_rect(fill = "gray10", colour = "gray10"),
          panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray80"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray80"),
          # axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(face = "bold")))

df = trophic %>% 
  subset(sp == "PLANKTIVORE") %>%
  mutate(lon = round(lon, 1),
         lat = round(lat, 1)) %>%
  group_by(lon, lat) %>% 
  summarise(est = median(est)) %>%  
  mutate(abs_est = abs(est))

(m2 = ggplot() +
    geom_sf(data = world) +
    geom_point(data = df,  aes_string("lon", "lat", color = "est"), size = scale(df$abs_est)) +
    coord_sf(crs = st_crs(4135)) +
    scale_color_gradientn(colours = matlab.like(100), "") + 
    scale_color_viridis_c("") + 
    ggtitle(" ") +
    theme_classic() +
    annotate("text",  x = Inf, y = Inf, label = "Planktivore", vjust = 1.2, hjust = 1.2, size = 4) + 
    theme(legend.position = c(0, 0), 
          legend.justification = c(-0.1, -0.1),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.box.background = element_rect(fill = "transparent", colour = "transparent"),
          # legend.text = element_text(color = "white", size = 12),
          legend.key.size = unit(0.5, "cm"),
          # panel.background = element_rect(fill = "gray10", colour = "gray10"),
          panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray80"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray80"),
          # axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(face = "bold")))

df = trophic %>% 
  subset(sp == "PRIMARY") %>%
  mutate(lon = round(lon, 1),
         lat = round(lat, 1)) %>%
  group_by(lon, lat) %>% 
  summarise(est = median(est)) %>%  
  mutate(abs_est = abs(est))

(m3 = ggplot() +
    geom_sf(data = world) +
    geom_point(data = df,  aes_string("lon", "lat", color = "est"), size = scale(df$abs_est)) +
    coord_sf(crs = st_crs(4135)) +
    scale_color_gradientn(colours = matlab.like(100), "") + 
    scale_color_viridis_c("") + 
    ggtitle(" ") +
    theme_classic() +
    annotate("text",  x = Inf, y = Inf, label = "Primary", vjust = 1.2, hjust = 1.2, size = 4) + 
    theme(legend.position = c(0, 0), 
          legend.justification = c(-0.1, -0.1),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.box.background = element_rect(fill = "transparent", colour = "transparent"),
          # legend.text = element_text(color = "white", size = 12),
          legend.key.size = unit(0.5, "cm"),
          # panel.background = element_rect(fill = "gray10", colour = "gray10"),
          panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray80"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray80"),
          # axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(face = "bold")))

df = trophic %>% 
  subset(sp == "SECONDARY") %>%
  mutate(lon = round(lon, 1),
         lat = round(lat, 1)) %>%
  group_by(lon, lat) %>% 
  summarise(est = median(est)) %>%  
  mutate(abs_est = abs(est))

(m4 = ggplot() +
    geom_sf(data = world) +
    geom_point(data = df,  aes_string("lon", "lat", color = "est"), size = scale(df$abs_est)) +
    coord_sf(crs = st_crs(4135)) +
    scale_color_gradientn(colours = matlab.like(100), "") + 
    scale_color_viridis_c("") + 
    ggtitle(" ") +
    theme_classic() +
    annotate("text",  x = Inf, y = Inf, label = "Secondary", vjust = 1.2, hjust = 1.2, size = 4) + 
    theme(legend.position = c(0, 0), 
          legend.justification = c(-0.1, -0.1),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.box.background = element_rect(fill = "transparent", colour = "transparent"),
          # legend.text = element_text(color = "white", size = 12),
          legend.key.size = unit(0.5, "cm"),
          # panel.background = element_rect(fill = "gray10", colour = "gray10"),
          panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray80"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray80"),
          # axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(face = "bold")))

png('outputs/fig2b.png', height = 3, width = 17, units = "in", res = 500)
grid.arrange(m1, m2, m3, m4, nrow = 1)
dev.off()

(fig2_supp = trophic %>% 
    group_by(lon, lat, sp) %>% 
    summarise(est = median(est)) %>%  
    ggplot(aes(lon, lat, fill = est)) + 
    geom_tile(height = 0.01, width = 0.01) +
    coord_fixed() + 
    facet_wrap(~sp, nrow = 2) +
    scale_fill_gradientn("", colors = matlab.like(100)) +
    theme(
      # axis.line = element_blank(),
      # axis.text = element_blank(),
      # axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.background = element_rect(fill = "gray10", colour = "gray10"),
      panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "gray20"), 
      panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "gray20")))

# png('outputs/fig2_supp.png', height = 7, width = 10, units = "in", res = 500)
# (fig2_supp)
# dev.off()

(fig3 = trophic %>%
    group_by(year, sp) %>% 
    summarise(mean_est = median(est),
              sd = sd(est, na.rm = T)) %>%
    ggplot(aes(year, mean_est)) +
    geom_pointrange(
      aes(ymin = ifelse(mean_est - sd < 0, 0, mean_est - sd), 
          ymax = mean_est+sd, color = sp),
      position = position_dodge(0.5)) + 
    # scale_color_discrete("", ) + 
    scale_color_viridis_d("") +
    # scale_color_manual(values = terrain.colors(4), "") +
    # scale_x_continuous(expand = c(0, 0)) +
    # scale_y_continuous(expand = c(0, 0)) +
    labs(y = expression("Biomass (g) per " ~ m^2~""), x = "") +
    scale_x_continuous(breaks = c(2010, 2012, 2013, 2015, 2016, 2019), 
                       labels = c(2010, 2012, 2013, 2015, 2016, 2019)) + 
    theme_classic() +
    theme(legend.position = c(1,1),
          legend.justification = c(1, 0.8),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          # panel.background = element_rect(fill = "gray90", colour = "gray90"),
          panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray100"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray100")))

png("outputs/fig3.png", height = 5, width = 6, units = "in", res = 500)
(fig3)
dev.off()

(fig3_supp = trophic %>% 
    # subset(sp == "SECONDARY") %>% 
    group_by(year, sp) %>% 
    summarise(mean_est = median(est),
              sd = sd(est, na.rm = T)) %>%
    ggplot(aes(year, mean_est)) + 
    geom_line(show.legend = F) +
    geom_point(size = 3, show.legend = F) + 
    geom_errorbar(aes(ymin = ifelse(mean_est - sd < 0, 0, mean_est - sd), 
                      ymax = mean_est+sd), 
                  width=.2,
                  show.legend = F,
                  position = position_dodge(0.05)) + 
    labs(y = expression("Biomass (g) per " ~ m^2~""), x = "") +
    facet_wrap(~sp) + 
    theme_half_open() +
    theme(legend.position = "right"))

# png('outputs/fig4_supp.png', height = 5, width = 8, units = "in", res = 500)
# (fig4_supp)
# dev.off()

