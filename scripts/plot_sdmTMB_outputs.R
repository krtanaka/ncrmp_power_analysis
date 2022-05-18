library(ggplot2)
library(dplyr)
library(colorRamps)
library(ggpubr)
library(gridExtra)
library(sp)

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
    coord_sf(crs = st_crs(4135), # old hawaii projection code
             xlim = c(-160.5, -154.7),
             ylim = c(18.8, 22.5), expand = F) +
    scale_color_distiller(palette ="RdBu",
                          direction = -1,
                          # limits = c(quantile(df$est, 0.999)*-1, quantile(df$est, 0.999)) +
                          "") + 
    ggtitle("Piscivore") + 
    theme_half_open() + 
    # guides(color=guide_legend(), size = guide_legend()) + 
    theme(legend.position = c(0, 0), 
          legend.justification = c(-0.1, -0.1),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.box.background = element_rect(fill = "transparent", colour = "transparent"),
          legend.text = element_text(color = "white", size = 12),
          legend.key.size = unit(0.5, "cm"),
          panel.background = element_rect(fill = "gray10", colour = "gray10"),
          panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray20"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray20"),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())
)

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
    coord_sf(crs = st_crs(4135), # old hawaii projection code
             xlim = c(-160.5, -154.7),
             ylim = c(18.8, 22.5), expand = F) +
    scale_color_distiller(palette ="RdBu",
                          direction = -1,
                          # limits = c(quantile(df$est, 0.999)*-1, quantile(df$est, 0.999)) +
                          "") + 
    ggtitle("Planktivore") + 
    theme_half_open() + 
    # guides(color=guide_legend(), size = guide_legend()) + 
    theme(legend.position = c(0, 0), 
          legend.justification = c(-0.1, -0.1),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.box.background = element_rect(fill = "transparent", colour = "transparent"),
          legend.text = element_text(color = "white", size = 12),
          legend.key.size = unit(0.5, "cm"),
          panel.background = element_rect(fill = "gray10", colour = "gray10"),
          panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray20"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray20"),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())
)

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
    coord_sf(crs = st_crs(4135), # old hawaii projection code
             xlim = c(-160.5, -154.7),
             ylim = c(18.8, 22.5), expand = F) +
    scale_color_distiller(palette ="RdBu",
                          direction = -1,
                          # limits = c(quantile(df$est, 0.999)*-1, quantile(df$est, 0.999)) +
                          "") + 
    ggtitle("Primary") + 
    theme_half_open() + 
    # guides(color=guide_legend(), size = guide_legend()) + 
    theme(legend.position = c(0, 0), 
          legend.justification = c(-0.1, -0.1),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.box.background = element_rect(fill = "transparent", colour = "transparent"),
          legend.text = element_text(color = "white", size = 12),
          legend.key.size = unit(0.5, "cm"),
          panel.background = element_rect(fill = "gray10", colour = "gray10"),
          panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray20"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray20"),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())
)

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
    coord_sf(crs = st_crs(4135), # old hawaii projection code
             xlim = c(-160.5, -154.7),
             ylim = c(18.8, 22.5), expand = F) +
    scale_color_distiller(palette ="RdBu",
                          direction = -1,
                          # limits = c(quantile(df$est, 0.999)*-1, quantile(df$est, 0.999)) +
                          "") + 
    ggtitle("Secondary") + 
    theme_half_open() + 
    # guides(color=guide_legend(), size = guide_legend()) + 
    theme(legend.position = c(0, 0), 
          legend.justification = c(-0.1, -0.1),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.box.background = element_rect(fill = "transparent", colour = "transparent"),
          legend.text = element_text(color = "white", size = 12),
          legend.key.size = unit(0.5, "cm"),
          panel.background = element_rect(fill = "gray10", colour = "gray10"),
          panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray20"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray20"),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())
)

png('outputs/fig2a.png', height = 7, width = 10, units = "in", res = 100)
grid.arrange(c1, c2, c3, c4, nrow = 2)
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

png('outputs/fig2_supp.png', height = 7, width = 10, units = "in", res = 100)
(c_supp)
dev.off()

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
    coord_sf(crs = st_crs(4135), # old hawaii projection code
             xlim = c(-160.5, -154.7),
             ylim = c(18.8, 22.5), expand = F) +
    scale_color_gradientn(colours = matlab.like(100), "") + 
    ggtitle("Piscivore") +
    theme_half_open() + 
    theme(legend.position = c(0, 0), 
          legend.justification = c(-0.1, -0.1),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.box.background = element_rect(fill = "transparent", colour = "transparent"),
          legend.text = element_text(color = "white", size = 12),
          legend.key.size = unit(0.5, "cm"),
          panel.background = element_rect(fill = "gray10", colour = "gray10"),
          panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray20"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray20"),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()))

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
    coord_sf(crs = st_crs(4135), # old hawaii projection code
             xlim = c(-160.5, -154.7),
             ylim = c(18.8, 22.5), expand = F) +
    scale_color_gradientn(colours = matlab.like(100), "") + 
    ggtitle("Planktivore") +
    theme_half_open() + 
    theme(legend.position = c(0, 0), 
          legend.justification = c(-0.1, -0.1),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.box.background = element_rect(fill = "transparent", colour = "transparent"),
          legend.text = element_text(color = "white", size = 12),
          legend.key.size = unit(0.5, "cm"),
          panel.background = element_rect(fill = "gray10", colour = "gray10"),
          panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray20"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray20"),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()))

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
    coord_sf(crs = st_crs(4135), # old hawaii projection code
             xlim = c(-160.5, -154.7),
             ylim = c(18.8, 22.5), expand = F) +
    scale_color_gradientn(colours = matlab.like(100), "") + 
    ggtitle("Primary") +
    theme_half_open() + 
    theme(legend.position = c(0, 0), 
          legend.justification = c(-0.1, -0.1),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.box.background = element_rect(fill = "transparent", colour = "transparent"),
          legend.text = element_text(color = "white", size = 12),
          legend.key.size = unit(0.5, "cm"),
          panel.background = element_rect(fill = "gray10", colour = "gray10"),
          panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray20"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray20"),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()))

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
    coord_sf(crs = st_crs(4135), # old hawaii projection code
             xlim = c(-160.5, -154.7),
             ylim = c(18.8, 22.5), expand = F) +
    scale_color_gradientn(colours = matlab.like(100), "") + 
    ggtitle("Secondary") +
    theme_half_open() + 
    theme(legend.position = c(0, 0), 
          legend.justification = c(-0.1, -0.1),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.box.background = element_rect(fill = "transparent", colour = "transparent"),
          legend.text = element_text(color = "white", size = 12),
          legend.key.size = unit(0.5, "cm"),
          panel.background = element_rect(fill = "gray10", colour = "gray10"),
          panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray20"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray20"),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()))

(m4 = trophic %>% 
    subset(sp == "PRIMARY") %>% 
    group_by(x, y, sp) %>% 
    summarise(est = median(est)) %>%  
    ggplot(aes(x, y, fill = est)) + 
    geom_tile(height = 0.8, width = 0.8) +
    coord_fixed() + 
    facet_grid(~sp) +
    ylab("Northings (km)") + 
    xlab("Eastings (km)") + 
    # scale_fill_gradient2("") +
    # scale_fill_viridis_c("") +
    scale_fill_gradientn("", colors = matlab.like(100)) +
    theme_half_open() +
    theme(legend.position = c(0, 0), 
          legend.justification = c(-0.1, -0.1),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.text = element_text(color = "white", size = 12),
          legend.key.size = unit(0.5, "cm"),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.background = element_rect(fill = "gray10", colour = "gray10"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "gray20"), 
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "gray20"))+ 
    ggtitle(""))

png(paste0('/Users/', Sys.info()[7], '/Desktop/fig4b.png'), height = 3, width = 16, units = "in", res = 100)
grid.arrange(m1, m2, m3, m4, nrow = 2)
dev.off()

(m_supp = trophic %>% 
    group_by(x, y, sp) %>% 
    summarise(est = median(est)) %>%  
    ggplot(aes(x, y, fill = est)) + 
    geom_tile(height = 0.8, width = 0.8) +
    coord_fixed() + 
    facet_wrap(~sp, nrow = 2) +
    ylab("Northings (km)") + 
    xlab("Eastings (km)") + 
    scale_fill_gradientn("", colors = matlab.like(100)) +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.background = element_rect(fill = "gray10", colour = "gray10"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "gray20"), 
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "gray20")))

png('outputs/fig3_supp.png', height = 7, width = 10, units = "in", res = 100)
(m_supp)
dev.off()


(t1 = trophic %>% 
    subset(sp == "PISCIVORE") %>%
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
    ylab("g/353 sq.m") + 
    xlab("Year") + 
    facet_grid(~sp) + 
    theme_half_open() + 
    theme(legend.position = "right"))

(t2 = trophic %>% 
    subset(sp == "PLANKTIVORE") %>% 
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
    ylab("g/353 sq.m") + 
    xlab("Year") + 
    facet_grid(~sp) + 
    theme_half_open() + 
    theme(legend.position = "right"))

(t3 = trophic %>% 
    subset(sp == "PRIMARY") %>% 
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
    ylab("g/353 sq.m") + 
    xlab("Year") + 
    facet_grid(~sp) + 
    theme_half_open() + 
    theme(legend.position = "right"))

(t4 = trophic %>% 
    subset(sp == "SECONDARY") %>% 
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
    ylab("g/353 sq.m") + 
    xlab("Year") + 
    facet_grid(~sp) + 
    theme_half_open() + 
    theme(legend.position = "right"))

png(paste0('/Users/', Sys.info()[7], '/Desktop/trend.png'), height = 3, width = 18, units = "in", res = 100)
grid.arrange(t1, t2, t3, t4, nrow = 1)
dev.off()

