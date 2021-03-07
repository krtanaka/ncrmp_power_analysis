library(mapplots)
basemap(xlim, ylim, main = "Gadoid landings")
data(coast)
data(landings)
byx = 1
byy = 0.5
xlim <- c(-15.5,0)
ylim <- c(50.25,56)
grd <- make.grid(landings$Lon,landings$Lat,landings$LiveWeight, byx, byy, xlim, ylim)
breaks <- breaks.grid(grd,zero=FALSE)
draw.grid(grd,breaks)
draw.shape(coast, col="darkgreen")
legend.grid("topright", breaks=breaks/1000, type=2, inset=0.02, title="tonnes")

require(HiClimR)

## Load test case data
x <- TestCase$x

## Generate longitude and latitude mesh vectors
xGrid <- grid2D(lon = unique(TestCase$lon), lat = unique(TestCase$lat))
lon <- c(xGrid$lon)
lat <- c(xGrid$lat)

library(maps)
library(mapproj)
m <- map("usa", plot=FALSE)
map("usa", project="albers", par=c(39, 45))
map.grid(m)

# get unprojected world limits
m <- map('world', plot=FALSE)

# center on NYC
map('world', proj='azequalarea', orient=c(41, -74, 0))
map.grid(m, col=2)
points(mapproject(list(y=41, x=-74)), col=3, pch="x", cex=2)

map('world', proj='orth', orient=c(41, -74, 0))
map.grid(m, col=2, nx=6, ny=5, label=FALSE, lty=2)
points(mapproject(list(y=41, x=-74)), col=3, pch="x", cex=2)

# center on Auckland
map('world', proj='orth', orient=c(-36.92, 174.6, 0))
map.grid(m, col=2, label=FALSE, lty=2)
points(mapproject(list(y=-36.92, x=174.6)), col=3, pch="x", cex=2)

m <- map('nz')
# center on Auckland
map('nz', proj='azequalarea', orient=c(-36.92, 174.6, 0))
points(mapproject(list(y=-36.92, x=174.6)), col=3, pch="x", cex=2)
map.grid(m, col=2)

library(dplyr)
library(RCzechia) # a package of Czech administrative areas
library(sf)


mesto <- kraje() %>% # All Czech NUTS3 ...
  filter(NAZ_CZNUTS3 == 'Hlavní město Praha') %>% # ... city of Prague
  st_transform(5514) # a metric CRS 

grid_spacing <- 1000  # size of squares, in units of the CRS (i.e. meters for 5514)

polygony <- st_make_grid(mesto, square = T, cellsize = c(grid_spacing, grid_spacing)) %>% # the grid, covering bounding box
  st_sf() # not really required, but makes the grid nicer to work with later

plot(polygony, col = 'white')
plot(st_geometry(mesto), add = T)
