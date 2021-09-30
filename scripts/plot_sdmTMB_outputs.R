library(ggplot2)
library(dplyr)
library(colorRamps)
library(ggpubr)
library(gridExtra)

list = list.files(path = "outputs/", pattern = "_biomass"); list

trophic = NULL

for (i in 3:6) {
  
  # i = 4
  
  load(paste0("outputs/", list[i]))
  sp = strsplit(list[i], split = "_")[[1]][3]; sp
  response_scale = strsplit(list[i], split = "_")[[1]][4]; response_scale
  sdm = sdm_output[,c("X", "Y", "longitude", "latitude", "year", "est", "zeta_s")]; rm(sdm_output)
  colnames(sdm)[1:2] = c("x", "y")
  sdm$est = exp(sdm$est); hist(sdm$est);summary(sdm$est)
  sdm$sp = sp
  trophic = rbind(trophic, sdm)
  
}

options(digits = 1)

(change1 = trophic %>% 
    subset(year == 2006) %>% 
    subset(sp == "PISCIVORE") %>% 
    group_by(x, y, sp) %>% 
    summarise(est = mean(zeta_s)) %>%  
    ggplot(aes(x, y, fill = est, color = est)) + 
    # geom_tile(height = 0.8, width = 0.8) +
    geom_point(alpha = 0.5, size = 0.5) +
    coord_fixed() + 
    facet_grid(~ sp) +
    ylab("Northings (km)") + 
    xlab("Eastings (km)") + 
    scale_fill_gradient2("Linear trend") +
    scale_color_gradient2("Linear trend") +
    # ggdark::dark_theme_minimal() +
    theme_minimal() + 
    theme(legend.position = c(0.15, 0.35)))

(change2 = trophic %>% 
    subset(year == 2006) %>% 
    subset(sp == "PLANKTIVORE") %>% 
    group_by(x, y, sp) %>% 
    summarise(est = mean(zeta_s)) %>%  
    ggplot(aes(x, y, fill = est, color = est)) + 
    # geom_tile(height = 0.8, width = 0.8) +
    geom_point(alpha = 0.5, size = 0.5) +
    coord_fixed() + 
    facet_grid(~ sp) + 
    ylab("Northings (km)") + 
    xlab("Eastings (km)") + 
    scale_fill_gradient2("Linear trend") +
    scale_color_gradient2("Linear trend") +
    # ggdark::dark_theme_minimal() +
    theme_minimal() + 
    theme(legend.position = c(0.15, 0.35)))

(change3 = trophic %>% 
    subset(year == 2006) %>% 
    subset(sp == "PRIMARY") %>% 
    group_by(x, y, sp) %>% 
    summarise(est = mean(zeta_s)) %>%  
    ggplot(aes(x, y, fill = est, color = est)) + 
    # geom_tile(height = 0.8, width = 0.8) +
    geom_point(alpha = 0.5, size = 0.5) +
    coord_fixed() + 
    facet_grid(~ sp) + 
    ylab("Northings (km)") + 
    xlab("Eastings (km)") + 
    scale_fill_gradient2("Linear trend") +
    scale_color_gradient2("Linear trend") +
    # ggdark::dark_theme_minimal() +
    theme_minimal() + 
    theme(legend.position = c(0.15, 0.35)))

(change4 = trophic %>% 
    subset(year == 2006) %>% 
    subset(sp == "SECONDARY") %>% 
    group_by(x, y, sp) %>% 
    summarise(est = mean(zeta_s)) %>%  
    ggplot(aes(x, y, fill = est, color = est)) + 
    # geom_tile(height = 0.8, width = 0.8) +
    geom_point(alpha = 0.5, size = 0.5) +
    coord_fixed() + 
    facet_grid(~ sp) + 
    ylab("Northings (km)") + 
    xlab("Eastings (km)") + 
    scale_fill_gradient2("Linear trend") +
    scale_color_gradient2("Linear trend") +
    # ggdark::dark_theme_minimal() +
    theme_minimal() + 
    theme(legend.position = c(0.15, 0.35)))

png('/Users/Kisei/Desktop/change.png', height = 3, width = 18, units = "in", res = 100)
grid.arrange(change1, change2, change3, change4, nrow = 1)
dev.off()

(map1 = trophic %>% 
    subset(sp == "PISCIVORE") %>% 
    group_by(x, y, sp) %>% 
    summarise(est = median(est)) %>%  
    ggplot(aes(x, y, fill = est, color = est)) + 
    # geom_tile(height = 0.8, width = 0.8) +
    geom_point(alpha = 0.5, size = 0.5) +
    coord_fixed() + 
    facet_grid(~sp) +
    ylab("Northings (km)") + 
    xlab("Eastings (km)") + 
    scale_fill_gradientn(colours = matlab.like(100), "g/353 sq.m") + 
    scale_color_gradientn(colours = matlab.like(100), "g/353 sq.m") + 
    # ggdark::dark_theme_minimal() +
    theme_minimal() + 
    theme(legend.position = c(0.15, 0.35)))

(map2 = trophic %>% 
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
    scale_fill_gradientn(colours = matlab.like(100), "g/353 sq.m") + 
    scale_color_gradientn(colours = matlab.like(100), "g/353 sq.m") + 
    # ggdark::dark_theme_minimal() +
    theme_minimal() + 
    theme(legend.position = c(0.15, 0.35)))

(map3 = trophic %>% 
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
    scale_fill_gradientn(colours = matlab.like(100), "g/353 sq.m") + 
    scale_color_gradientn(colours = matlab.like(100), "g/353 sq.m") + 
    # ggdark::dark_theme_minimal() +
    theme_minimal() + 
    theme(legend.position = c(0.15, 0.35)))

(map4 = trophic %>% 
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
    scale_fill_gradientn(colours = matlab.like(100), "g/353 sq.m") + 
    scale_color_gradientn(colours = matlab.like(100), "g/353 sq.m") + 
    # ggdark::dark_theme_minimal() +
    theme_minimal() + 
    theme(legend.position = c(0.15, 0.35)))

png('/Users/Kisei/Desktop/biomass.png', height = 3, width = 18, units = "in", res = 100)
grid.arrange(map1, map2, map3, map4, nrow = 1)
dev.off()




(trend1 = trophic %>% 
    subset(sp == "PISCIVORE") %>% 
    group_by(year, sp) %>% 
    summarise(mean_est = median(est),
              sd = sd(est, na.rm = T)) %>%
    ggplot(aes(year, mean_est, color = mean_est)) + 
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

(trend2 = trophic %>% 
    subset(sp == "PLANKTIVORE") %>% 
    group_by(year, sp) %>% 
    summarise(mean_est = median(est),
              sd = sd(est, na.rm = T)) %>%
    ggplot(aes(year, mean_est, color = mean_est)) + 
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

(trend3 = trophic %>% 
    subset(sp == "PRIMARY") %>% 
    group_by(year, sp) %>% 
    summarise(mean_est = median(est),
              sd = sd(est, na.rm = T)) %>%
    ggplot(aes(year, mean_est, color = mean_est)) + 
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

(trend4 = trophic %>% 
    subset(sp == "SECONDARY") %>% 
    group_by(year, sp) %>% 
    summarise(mean_est = median(est),
              sd = sd(est, na.rm = T)) %>%
    ggplot(aes(year, mean_est, color = mean_est)) + 
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

png('/Users/Kisei/Desktop/trend.png', height = 3, width = 18, units = "in", res = 100)
grid.arrange(trend1, trend2, trend3, trend4, nrow = 1)
dev.off()
