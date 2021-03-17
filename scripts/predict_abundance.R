rm(list = ls())

library(ggplot2)
library(dplyr)
library(maps)
library(ggdark)


load("data/HAW_Grid.RData")
grid = Hawaii_Survey_Grid
grid$DEPTH = ifelse(grid$DEPTH < -2000, grid$DEPTH_e, grid$DEPTH)
# grid$DEPTH[grid$DEPTH < -2000] <- NA

load("data/ALL_REA_FISH_RAW.rdata")
biomass = df %>% subset(ISLAND == "Hawaii"); rm(df)

load("data/SURVEY MASTER.RData")
survey = SURVEY_MASTER %>% subset(ISLAND == "Hawaii"); rm(SURVEY_MASTER)

grid %>% 
  ggplot(aes(x = X, y = Y, fill = DEPTH)) +
  geom_tile(width = 0.005, height = 0.005) + 
  scale_fill_viridis_c() + 
  coord_fixed()

survey %>% 
  ggplot(aes(LONGITUDE_LOV, LATITUDE_LOV, color = OBS_YEAR)) + 
  geom_point() + 
  scale_color_viridis_c()+ 
  coord_fixed()

grid %>% 
  as.data.frame() %>% 
  group_by(DEPTH_BIN) %>% 
  summarise(sq.km = n()*50) %>% 
  ggplot(aes(x = DEPTH_BIN, y = sq.km, fill = DEPTH_BIN)) +
  geom_bar(stat = "identity")+
  theme_minimal() + 
  scale_fill_viridis_d()

library(mgcv)

gam_df = biomass %>% 
  group_by(LONGITUDE, LATITUDE, DEPTH) %>% 
  summarise(cpue = mean(BIOMASS_G_M2)) %>% 
  rename(lon = LONGITUDE,
         lat = LATITUDE,
         depth = DEPTH)

fish_gam = gam(cpue ~ s(lat, lon, k = 30) + s(depth, k = 5), 
               data = gam_df, 
               family = "tw(theta = NULL, link = 'log',a=1.01,b=1.99)",
               gamma = 1.4)

summary(fish_gam)

plot(fish_gam, pages = 1, shade = T, res = F, all.terms = T)
vis.gam(fish_gam, view=c("lon", "lat"), plot.type = "contour", color = "cm", too.far = 0.03, main = "")
map(add = T, resolution = 0, fill = T)
dsm::rqgam.check(fish_gam, pch = ".")
gam.check(fish_gam)

detach("package:raster", unload = TRUE)
predict_env = grid %>% as.data.frame() %>% 
  select(X, Y, DEPTH) %>% 
  rename(lon = X,
         lat = Y,
         depth = DEPTH) %>% 
  mutate(depth = depth*-1)

pred_biomass = data.frame(predict(fish_gam, type = "response", newdata = predict_env, se.fit = T))

predict_env = data.frame(predict_env, pred_biomass) 

predict_env %>% ggplot(aes(lon, lat, fill = fit)) +
  geom_tile(width = 0.005, height = 0.005) + 
  scale_fill_viridis_c("g/sq.m") + 
  dark_theme_minimal() + 
  coord_fixed()

