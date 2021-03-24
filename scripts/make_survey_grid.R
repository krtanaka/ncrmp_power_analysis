rm(list = ls())

library(SimSurvey)
library(raster)

plot(survey_grid)

load("data/HAW_Grid.RData")

df = Hawaii_Survey_Grid %>% 
  subset(DEPTH_e > -150) %>% 
  group_by(X, Y) %>% 
  summarise(X = mean(X),
            Y = mean(Y),
            depth = mean(DEPTH_e, na.rm = T))

plot(df$X, df$Y, pch = ".")

library(sp)
library(rgdal)
coordinates(df)=~X+Y

proj4string(df)=CRS("+init=epsg:4326") # set it to lat-long
df = spTransform(df,CRS("+proj=longlat +datum=NAD27"))

gridded(df) = TRUE

  
  r = raster(pts)
projection(r) = CRS("insert your proj4 string here")
Now have a look:
  
  plot(r)