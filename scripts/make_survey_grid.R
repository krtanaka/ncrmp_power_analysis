library(readr)
library(dplyr)
library(ggplot2)

df <- read_csv("data/jplmursst_strm15.csv")
# df$bathymetry_m = ifelse(df$bathymetry_m %in% c(-30:0), df$bathymetry_m, NA)
df$bathymetry_m = ifelse(df$bathymetry_m > 0, NA, df$bathymetry_m)

survey = read_csv("data/SURVEY MASTER.csv")
survey$longitude = survey$LONGITUDE_LOV
survey$latitude = survey$LATITUDE_LOV
survey$longitude = ifelse(survey$longitude < 0, survey$longitude + 360, survey$longitude)
survey = subset(survey, ISLAND == "Hawaii")

ggplot() + 
  geom_tile(data = df, aes(longitude, latitude, fill = bathymetry_m)) +
  geom_point(data = survey, aes(longitude, latitude, color = as.factor(OBS_YEAR))) +
  scale_fill_viridis_c() + 
  scale_color_manual(values = cm.colors(9),"") +
  coord_fixed()

