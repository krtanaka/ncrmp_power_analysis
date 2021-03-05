rm(list = ls())

library(readr)
library(dplyr)
library(rnaturalearth)
library(rgeos)
library(sf)
library(ggplot2)
library(stringr)

df <- read_csv("data/SURVEY MASTER.csv") %>% subset(ISLAND == "Hawaii")

df %>% group_by(
  DEPTH_BIN,
  REEF_ZONE,
  OBS_YEAR) %>% 
  summarise(n = n()) %>% 
  na.omit() %>% 
  ggplot(aes(OBS_YEAR, n, color = n)) + 
  geom_point() + 
  scale_fill_viridis_c() + 
  facet_grid(DEPTH_BIN ~ REEF_ZONE) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

df$LONGITUDE_LOV = ifelse(df$LONGITUDE_LOV < 0, 
                          df$LONGITUDE_LOV + 360, 
                          df$LONGITUDE_LOV)

time_stamp = str_split_fixed(df$DATE_, "/", 3)
colnames(time_stamp) = c("Month", "Day", "Year")

df = cbind(time_stamp, df)

df = df %>% 
  group_by(LATITUDE_LOV, LONGITUDE_LOV, Day, Month, Year) %>% 
  summarise(n = n())

df = df[complete.cases(df), ]

world <- rnaturalearth::ne_countries(scale = 'small', returnclass = "sp")
box_cut <- bbox2SP(n = 90, s = -90, w = -20, e = 20, proj4string = world@proj4string)
world_crop <- gDifference(world, box_cut)

pacific_crop <- world_crop %>% 
  st_as_sf() %>% 
  st_shift_longitude() %>% 
  st_crop(c(
    # xmin = st_bbox(.)[["xmin"]],
    # xmax = st_bbox(.)[["xmax"]],
    xmin = min(df$LONGITUDE_LOV),
    xmax = max(df$LONGITUDE_LOV),
    ymin = min(df$LATITUDE_LOV),
    ymax = max(df$LATITUDE_LOV)))

df$Month = ifelse(df$Month %in% c(1:9), sprintf("%02d", as.numeric(df$Month)), df$Month)

ggplot() +
  geom_sf(data = pacific_crop, size = 0.01, fill = "gray", color = "gray") +
  geom_point(data = df, aes(x = LONGITUDE_LOV, y = LATITUDE_LOV, color = log(n)), alpha = 0.5, size = 3) +
  scale_color_viridis_c() + 
  # facet_grid(Month ~ Year) +
  facet_wrap(~Year) +
  # facet_wrap(~Month) +
  # scale_x_continuous(expand = c(0, 0), "") +
  # scale_y_continuous(expand = c(0, 0), "") + 
  theme_void()
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
