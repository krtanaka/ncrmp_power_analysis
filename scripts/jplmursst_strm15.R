library(readr)
library(dplyr)
library(ggplot2)

df <- read_csv("data/jplmursst_strm15.csv")
df$bathymetry_m = ifelse(df$bathymetry_m %in% c(-30:0), df$bathymetry_m, NA)
# df$bathymetry_m = ifelse(df$bathymetry_m > 0, NA, df$bathymetry_m)

df$cell = df$X1
df$division = 1
df$strat = 2
df$strat = ifelse(df$bathymetry_m %in% c(-10:0), 1, df$strat)
df$strat = ifelse(df$bathymetry_m %in% c(-30:-20), 3, df$strat)
df$depth = df$bathymetry_m*-1

df = df %>% drop_na()

ggplot() + 
  geom_raster(data = df, aes(longitude, latitude, fill = depth)) +
  scale_fill_viridis_c() + 
  coord_fixed() + 
  ggdark::dark_theme_void()

cell = rasterFromXYZ(df[,c("longitude", "latitude", "cell")]); plot(cell)
division = rasterFromXYZ(df[,c("longitude", "latitude", "division")]); plot(division)
strat = rasterFromXYZ(df[,c("longitude", "latitude", "strat")]); plot(strat)
depth = rasterFromXYZ(df[,c("longitude", "latitude", "depth")]); plot(depth)

survey_grid_kt = stack(cell, division, strat, depth)

sp::spplot(survey_grid)
sp::spplot(survey_grid_kt)
