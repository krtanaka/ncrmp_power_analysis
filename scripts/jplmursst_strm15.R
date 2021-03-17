library(readr)
library(dplyr)
library(ggplot2)

df <- read_csv("data/jplmursst_strm15.csv")
df$bathymetry_m = ifelse(df$bathymetry_m %in% c(-30:0), df$bathymetry_m, NA)
df$bathymetry_m = ifelse(df$bathymetry_m > 0, NA, df$bathymetry_m)

ggplot() + 
  geom_tile(data = df, aes(longitude, latitude, fill = bathymetry_m)) +
  scale_fill_viridis_c() + 
  scale_color_manual(values = cm.colors(9),"") +
  coord_fixed()