library(ggplot2)
library(dplyr)
library(colorRamps)

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


change = trophic %>% 
  subset(year == 2006) %>% 
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
  theme_pubr() + 
  theme(legend.position = "right")

map = trophic %>% 
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
  theme_pubr() + 
  theme(legend.position = "right")

trend = trophic %>% 
  group_by(year, sp) %>% 
  summarise(mean_est = median(est),
            sd = sd(est, na.rm = T)) %>%
  ggplot(aes(year, mean_est, color = mean_est)) + 
  geom_line(show.legend = F) +
  geom_point(show.legend = F) + 
  geom_errorbar(aes(ymin = ifelse(mean_est - sd < 0, 0, mean_est - sd), 
                    ymax = mean_est+sd), 
                width=.2,
                show.legend = F,
                position = position_dodge(0.05)) + 
  ylab("g/353 sq.m") + 
  xlab("Year") + 
  facet_grid(~sp) + 
  # ggdark::dark_theme_minimal() +
  theme_pubr() + 
  theme(legend.position = "right")

library(patchwork)
change/map
