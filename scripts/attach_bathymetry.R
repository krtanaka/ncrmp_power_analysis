rm(list = ls())

library(spatial)
library(raster)
library(lubridate)
library(raster)
library(readr)

Hawaii_Survey_Grid = read_csv("T:/Fish/GIS/Projects/Gridding/HAW/HAW_sites.csv")

source("/Users/Kisei.Tanaka/env_data_summary/scripts/HelperCode/ExpandingExtract.R")

SM = as.data.frame(Hawaii_Survey_Grid)
# SM$X = ifelse(SM$X < 0, SM$X + 360, SM$X)

SM = SM[complete.cases(SM[,c("X", "Y")]), ]
SM_sp = SM; SM_sp = as.data.frame(SM_sp)
coordinates(SM_sp) = ~X + Y

this_r = raster("M:/Environmental Data Summary/DataDownload/Bathymetry_SRTM15/Bathymetry_SRTM15_Bathy_M_AllIslands.nc")

e <- extent(range(pretty(SM$X))[1] , 
            range(pretty(SM$X))[2], 
            range(pretty(SM$Y))[1],
            range(pretty(SM$Y))[2])

this_r <- crop(this_r, e)
this_r[this_r > 0] <- NA

crs(SM_sp) = crs(this_r)

this_Ex = ExpandingExtract(this_r, SM_sp, Dists = c(0, 50, 100, 1000, 2000, 4000, 8000))

SM_sp$DEPTH_e = this_Ex$values
Hawaii_Survey_Grid = as.data.frame(SM_sp)
save(Hawaii_Survey_Grid, file = "data/HAW_Grid.RData")
