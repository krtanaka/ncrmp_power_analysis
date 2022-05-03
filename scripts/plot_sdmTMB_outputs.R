library(ggplot2)
library(dplyr)
library(colorRamps)
library(ggpubr)
library(gridExtra)

rm(list = ls())

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

(c1 = trophic %>% 
    # subset(year == 2010) %>% 
    # subset(sp == "PISCIVORE") %>% 
    group_by(x, y, sp) %>% 
    summarise(est = mean(zeta_s)) %>%  
    ggplot(aes(x, y, fill = est, color = est)) + 
    geom_tile(height = 0.8, width = 0.8) +
    coord_fixed() + 
    facet_grid(~ sp) +
    ylab("Northings (km)") + 
    xlab("Eastings (km)") + 
    scale_fill_gradient2("Linear trend") +
    scale_color_gradient2("Linear trend") +
    # scale_fill_viridis_c("Linear trend") + 
    # scale_color_viridis_c("Linear trend") + 
    # ggdark::dark_theme_minimal() +
    theme_half_open() +
    ggtitle("(a)"))

(c1 = trophic %>% 
    subset(year == 2010) %>% 
    subset(sp == "PISCIVORE") %>% 
    group_by(x, y, sp) %>% 
    summarise(est = mean(zeta_s)) %>%  
    ggplot(aes(x, y, fill = est, color = est)) + 
    geom_tile(height = 0.8, width = 0.8) +
    # geom_point(alpha = 0.5, size = 0.5) +
    coord_fixed() + 
    facet_grid(~ sp) +
    ylab("Northings (km)") + 
    xlab("Eastings (km)") + 
    # scale_fill_gradient2("Linear trend") +
    # scale_color_gradient2("Linear trend") +
    scale_fill_viridis_c("Linear trend") + 
    scale_color_viridis_c("Linear trend") + 
    # ggdark::dark_theme_minimal() +
    theme_half_open() +
    theme(legend.position = c(0, 0), 
          legend.justification = c(-0.1, -0.1))+ 
    ggtitle("(a)"))

(c2 = trophic %>% 
    subset(year == 2010) %>% 
    subset(sp == "PLANKTIVORE") %>% 
    group_by(x, y, sp) %>% 
    summarise(est = mean(zeta_s)) %>%  
    ggplot(aes(x, y, fill = est, color = est)) + 
    geom_tile(height = 0.8, width = 0.8) +
    # geom_point(alpha = 0.5, size = 0.5) +
    coord_fixed() + 
    facet_grid(~ sp) + 
    ylab("Northings (km)") + 
    xlab("Eastings (km)") + 
    # scale_fill_gradient2("Linear trend") +
    # scale_color_gradient2("Linear trend") +
    scale_fill_viridis_c("Linear trend") + 
    scale_color_viridis_c("Linear trend") + 
    # ggdark::dark_theme_minimal() +
    theme_half_open() +
    theme(legend.position = c(0, 0), 
          legend.justification = c(-0.1, -0.1))+ 
    ggtitle(""))

(c3 = trophic %>% 
    subset(year == 2010) %>% 
    subset(sp == "PRIMARY") %>% 
    group_by(x, y, sp) %>% 
    summarise(est = mean(zeta_s)) %>%  
    ggplot(aes(x, y, fill = est, color = est)) + 
    geom_tile(height = 0.8, width = 0.8) +
    # geom_point(alpha = 0.5, size = 0.5) +
    coord_fixed() + 
    facet_grid(~ sp) + 
    ylab("Northings (km)") + 
    xlab("Eastings (km)") + 
    # scale_fill_gradient2("Linear trend") +
    # scale_color_gradient2("Linear trend") +
    scale_fill_viridis_c("Linear trend") + 
    scale_color_viridis_c("Linear trend") + 
    # ggdark::dark_theme_minimal() +
    theme_half_open() +
    theme(legend.position = c(0, 0), 
          legend.justification = c(-0.1, -0.1))+ 
    ggtitle(""))

(c4 = trophic %>% 
    subset(year == 2010) %>% 
    subset(sp == "SECONDARY") %>% 
    group_by(x, y, sp) %>% 
    summarise(est = mean(zeta_s)) %>%  
    ggplot(aes(x, y, fill = est, color = est)) + 
    geom_tile(height = 0.8, width = 0.8) +
    # geom_point(alpha = 0.5, size = 0.5) +
    coord_fixed() + 
    facet_grid(~ sp) + 
    ylab("Northings (km)") + 
    xlab("Eastings (km)") + 
    # scale_fill_gradient2("Linear trend") +
    # scale_color_gradient2("Linear trend") +
    scale_fill_viridis_c("Linear trend") + 
    scale_color_viridis_c("Linear trend") + 
    # ggdark::dark_theme_minimal() +
    theme_half_open() +
    theme(legend.position = c(0, 0), 
          legend.justification = c(-0.1, -0.1))+ 
    ggtitle(""))

png(paste0('/Users/', Sys.info()[7], '/Desktop/fig4a.png'), height = 3, width = 18, units = "in", res = 100)
grid.arrange(c1, c2, c3, c4, nrow = 1)
dev.off()

(m1 = trophic %>% 
    subset(sp == "PISCIVORE") %>% 
    group_by(x, y, sp) %>% 
    summarise(est = median(est)) %>%  
    ggplot(aes(x, y, fill = est, color = est)) + 
    geom_tile(height = 0.8, width = 0.8) +
    # geom_point(alpha = 0.5, size = 0.5) +
    coord_fixed() + 
    facet_grid(~sp) +
    ylab("Northings (km)") + 
    xlab("Eastings (km)") + 
    # scale_fill_gradientn(colours = matlab.like(100), "g/353 sq.m") + 
    # scale_color_gradientn(colours = matlab.like(100), "g/353 sq.m") + 
    scale_fill_viridis_c("g/353 sq.m") + 
    scale_color_viridis_c("g/353 sq.m") + 
    # ggdark::dark_theme_minimal() +
    theme_minimal() + 
    theme(legend.position = c(0.15, 0.35)))

(m2 = trophic %>% 
    subset(sp == "PLANKTIVORE") %>% 
    group_by(x, y, sp) %>% 
    summarise(est = median(est)) %>%  
    ggplot(aes(x, y, fill = est, color = est)) + 
    geom_tile(height = 0.8, width = 0.8) +
    # geom_point(alpha = 0.5, size = 0.5) +
    coord_fixed() + 
    facet_grid(~sp) +
    ylab("Northings (km)") + 
    xlab("Eastings (km)") + 
    # scale_fill_gradientn(colours = matlab.like(100), "g/353 sq.m") + 
    # scale_color_gradientn(colours = matlab.like(100), "g/353 sq.m") + 
    scale_fill_viridis_c("g/353 sq.m") + 
    scale_color_viridis_c("g/353 sq.m") + 
    # ggdark::dark_theme_minimal() +
    theme_minimal() + 
    theme(legend.position = c(0.15, 0.35)))

(m3 = trophic %>% 
    subset(sp == "PRIMARY") %>% 
    group_by(x, y, sp) %>% 
    summarise(est = median(est)) %>%  
    ggplot(aes(x, y, fill = est, color = est)) + 
    geom_tile(height = 0.8, width = 0.8) +
    # geom_point(alpha = 0.5, size = 0.5) +
    coord_fixed() + 
    facet_grid(~sp) +
    ylab("Northings (km)") + 
    xlab("Eastings (km)") + 
    # scale_fill_gradientn(colours = matlab.like(100), "g/353 sq.m") + 
    # scale_color_gradientn(colours = matlab.like(100), "g/353 sq.m") + 
    scale_fill_viridis_c("g/353 sq.m") + 
    scale_color_viridis_c("g/353 sq.m") + 
    # ggdark::dark_theme_minimal() +
    theme_minimal() + 
    theme(legend.position = c(0.15, 0.35)))

(m4 = trophic %>% 
    subset(sp == "SECONDARY") %>% 
    group_by(x, y, sp) %>% 
    summarise(est = median(est)) %>%  
    ggplot(aes(x, y, fill = est, color = est)) + 
    geom_tile(height = 0.8, width = 0.8) +
    # geom_point(alpha = 0.5, size = 0.5) +
    coord_fixed() + 
    facet_grid(~sp) +
    ylab("Northings (km)") + 
    xlab("Eastings (km)") + 
    # scale_fill_gradientn(colours = matlab.like(100), "g/353 sq.m") + 
    # scale_color_gradientn(colours = matlab.like(100), "g/353 sq.m") + 
    scale_fill_viridis_c("g/353 sq.m") + 
    scale_color_viridis_c("g/353 sq.m") + 
    # ggdark::dark_theme_minimal() +
    theme_minimal() + 
    theme(legend.position = c(0.15, 0.35)))

png(paste0('/Users/', Sys.info()[7], '/Desktop/biomass.png'), height = 3, width = 18, units = "in", res = 100)
grid.arrange(m1, m2, m3, m4, nrow = 1)
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
    # ggdark::dark_theme_minimal() +
    theme_minimal() + 
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
    # ggdark::dark_theme_minimal() +
    theme_minimal() + 
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
    # ggdark::dark_theme_minimal() +
    theme_minimal() + 
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
    # ggdark::dark_theme_minimal() +
    theme_minimal() + 
    theme(legend.position = "right"))

png(paste0('/Users/', Sys.info()[7], '/Desktop/trend.png'), height = 3, width = 18, units = "in", res = 100)
grid.arrange(t1, t2, t3, t4, nrow = 1)
dev.off()

