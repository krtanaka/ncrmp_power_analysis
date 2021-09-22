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
  
  # replace sim$ with sdmTMB outputs  -----------------------------------------
  
  # Need to match sim and sdm
  sdm = sdm_output[,c("X", "Y", "longitude", "latitude", "year", "est" )]; rm(sdm_output)
  colnames(sdm)[1:2] = c("x", "y")
  sdm$est = exp(sdm$est); hist(sdm$est);summary(sdm$est)
  # sdm$est = sdm$est*(res(survey_grid_kt)[1] * res(survey_grid_kt)[2]*1000000); hist(sdm$est); summary(sdm$est) # convert g/sq.m to q/ whatever given cell size
  
  
  sdm$sp = sp
  
  trophic = rbind(trophic, sdm)
  
}

map = trophic %>% 
  group_by(x, y, sp) %>% 
  summarise(est = median(est)) %>%  
  ggplot(aes(x, y, fill = est)) + 
  geom_tile(height = 0.8, width = 0.8) +
  coord_fixed() + 
  facet_grid(~sp) + 
  ylab("Northing (km)") + 
  xlab("Easting (km)") + 
  scale_fill_gradientn(colours = matlab.like(100), "g/353 sq.m") + 
  theme_minimal() 
# facet_wrap(.~sp, ncol = 4) +
# scale_fill_viridis_c() + 
# ggdark::dark_theme_minimal() + 
# theme(legend.position = "top")

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
  theme_minimal() #+ 
# facet_wrap(.~ sp, ncol = 4) +
# ggdark::dark_theme_minimal() + 
# coord_fixed(ratio = 1.5) + 
# theme(legend.position = c(0.9, 0.9))

library(patchwork)
map/trend
