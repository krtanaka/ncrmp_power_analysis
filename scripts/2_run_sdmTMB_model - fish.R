library(sdmTMB)
library(dplyr)
library(ggplot2)
library(rgdal)
library(colorRamps)
library(patchwork)
library(raster)
library(sf)

rm(list = ls())

load("data/rea/fish_site_data.Rdata")

islands = as.character(unique(wsd$ISLAND))[9]

unit = c("biomass", "abundance")[1]

response = c("PISCIVORE_BIO", "PLANKTIVORE_BIO", "PRIMARY_BIO", "SECONDARY_BIO", "TotFishBio")[5]

wsd$time = as.numeric(as.POSIXct(wsd$DATE_, format="%Y-%m-%d")) # unix time

wsd$response = response

df = wsd %>% 
  subset(ISLAND %in% islands) %>%
  # group_by(LONGITUDE, LATITUDE, time) %>%
  group_by(LONGITUDE, LATITUDE, OBS_YEAR) %>%
  summarise(response = mean(PISCIVORE_BIO, na.rm = T),
            depth = mean(DEPTH, na.rm = T)) %>%  
  subset(response < quantile(response, prob = 0.99))

df %>% ggplot(aes(response)) + geom_histogram() +
  df %>% group_by(OBS_YEAR) %>% summarise(n = mean(response)) %>% ggplot(aes(OBS_YEAR, n)) + geom_point() + geom_line()

zone <- (floor((df$LONGITUDE[1] + 180)/6) %% 60) + 1
xy_utm = as.data.frame(cbind(utm = project(as.matrix(df[, c("LONGITUDE", "LATITUDE")]), paste0("+proj=utm +units=km +zone=", zone))))
colnames(xy_utm) = c("X", "Y")
df = cbind(df, xy_utm)

# Read in Island Boundaries
load('data/misc/MHI_islands_shp.RData')
crs(ISL_bounds) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
ISL_this = ISL_bounds[which(ISL_bounds$ISLAND %in% toupper(islands)),]
ISL_this_utm = spTransform(ISL_this,CRS(paste0("+proj=utm +units=km +zone=",zone)))
ISL_this_sf = st_transform(st_as_sf(ISL_this), crs = paste0("+proj=utm +units=km +zone=",zone))

n_knots = 300
n_knots = 100 # a coarse mesh for speed
rea_spde <- make_mesh(df, c("X", "Y"), n_knots = n_knots, type = "cutoff_search") # search
# rea_spde <- make_mesh(df, c("X", "Y"), cutoff  = n_knots, type = "cutoff") # predefined

#build barrier to mesh
rea_spde_coast = add_barrier_mesh(rea_spde , ISL_this_sf)

plot(rea_spde_coast$mesh, asp = 1, pch = "."); axis(1); axis(2)
plot(ISL_this_utm, add = TRUE)
points(rea_spde_coast$loc_xy,col = "green", pch = ".", cex = 5)
bar_i = rea_spde_coast$barrier_triangles
norm_i = rea_spde_coast$normal_triangles
points(rea_spde_coast$spde$mesh$loc[,1], rea_spde_coast$spde$mesh$loc[,2], pch = ".", col = "black")
points(rea_spde_coast$mesh_sf$V1[bar_i], rea_spde_coast$mesh_sf$V2[bar_i], col = "red", pch = 4, cex = .5)
points(rea_spde_coast$mesh_sf$V1[norm_i], rea_spde_coast$mesh_sf$V2[norm_i], col = "blue", pch = 1, cex = .5)

df$year = df$OBS_YEAR
df$depth_scaled = scale(log(df$depth))
df$depth_scaled2 = df$depth_scaled ^ 2

plot(df$depth, df$response, pch = 20, bty = "n")

density_model <- sdmTMB(
  
  data = df, 
  
  # formula = response ~ as.factor(year) + depth_scaled + depth_scaled2,
  formula = response ~ as.factor(year) + s(depth, k = 5),
  
  silent = F, 
  # extra_time = missing_year,
  spatial_trend = T, 
  spatial_only = F, 
  time = "year", 
  spde = rea_spde, 
  anisotropy = T,
  family = tweedie(link = "log"),
  # family = poisson(link = "log"),
  # family = binomial(link = "logit"), weights = n,
  # family = nbinom2(link = "log"),
  # family = Beta(link = "logit"),
  
  control = sdmTMBcontrol(step.min = 0.01, step.max = 1)
  
); beepr::beep(2)

density_model <- run_extra_optimization(density_model); beepr::beep(2)

# look at gradients
max(density_model$gradients)

df$residuals <- residuals(density_model)
par(pty = "s")
qqnorm(df$residuals, ylim = c(-5, 5), xlim = c(-5, 5), bty = "n", pch = 20); abline(a = 0, b = 1)

m_p <- predict(density_model); m_p = m_p[,c("response", "est")]
m_p$back_abs_res = abs(df$residuals)

ggdark::invert_geom_defaults()

p1 = ggplot(df, aes_string("X", "Y", color = "residuals")) +
  geom_point(alpha = 0.8, size = round(abs(df$residuals), digits = 0)) + 
  xlab("Eastings") +
  ylab("Northings") + 
  scale_color_gradient2() + 
  theme_minimal()

p2 = m_p  %>% 
  ggplot(aes(response, exp(est))) + 
  geom_point(alpha = 0.2, aes(size = back_abs_res)) + 
  coord_fixed(ratio = 1) +
  ylab("prediction") + 
  xlab("observation") + 
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method = "lm", se = T) + 
  theme_minimal()

p1 / p2

#  extract some parameter estimates
sd <- as.data.frame(summary(TMB::sdreport(density_model$tmb_obj)))
r <- density_model$tmb_obj$report()
r

# prediction onto new data grid
load("data/gis_bathymetry/raster/gua.RData") # bathymetry 

utmcoor <- SpatialPoints(cbind(topo$x, topo$y), proj4string = CRS("+proj=utm +units=m +zone=55"))
longlatcoor <- spTransform(utmcoor,CRS("+proj=longlat"))
topo$x <- coordinates(longlatcoor)[,1]
topo$y <- coordinates(longlatcoor)[,2]
# topo$x = ifelse(topo$x > 180, topo$x - 360, topo$x)

grid = topo
grid <- topo %>% subset(x > range(pretty(df$LONGITUDE))[1] 
                        & x < range(pretty(df$LONGITUDE))[2] 
                        & y > range(pretty(df$LATITUDE))[1] 
                        & y < range(pretty(df$LATITUDE))[2])
res = 2
grid$longitude = round(grid$x, digits = res)
grid$latitude = round(grid$y, digits = res)

grid$longitude = grid$x
grid$latitude = grid$y

zone <- (floor((grid$longitude[1] + 180)/6) %% 60) + 1
xy_utm = as.data.frame(cbind(utm = project(as.matrix(grid[, c("longitude", "latitude")]),
                                           paste0("+proj=utm +units=km +zone=", zone))))

colnames(xy_utm) = c("X", "Y"); plot(xy_utm, pch = ".")

grid = cbind(grid, xy_utm)

grid_year = NULL

years = sort(as.vector(unique(df$year)))

for (y in 1:length(years)) {
  
  # y = 1
  
  grid_y = grid[,c("X", "Y", "depth")]
  grid_y$depth = grid_y$depth *-1
  grid_y$year = years[[y]]
  grid_year = rbind(grid_year, grid_y)
  
}

grid_year$depth_scaled = scale(log(grid_year$depth+0.0001))
grid_year$depth_scaled2 = grid_year$depth_scaled ^ 2

# set the area argument to 0.01 km2 since our grid cells are 100 m x 100 m = 0.01 square kilometers
p <- predict(density_model, 
             newdata = grid_year, 
             return_tmb_object = T, area = 0.01)

p$data$response = response
sdm_output = p$data

save(sdm_output, file = paste0("outputs/sdmTMB_results_", response, "_", unit, "_", n_knots, "_", ".RData"))

plot_map_raster <- function(dat, column = "est") {
  
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_tile(aes(height = 0.5, width = 0.5), alpha = 0.5) +
    # geom_raster() +
    facet_wrap(~year) +
    coord_fixed() +
    xlab("Eastings (km)") +
    ylab("Northings (km)") + 
    scale_fill_gradientn(colours = matlab.like(100), "") +
    ggdark::dark_theme_minimal()
  
}

# pick out a single year to plot since they should all be the same for the slopes. Note that these are in log space.
plot_map_raster(filter(p$data, year == unique(p$data$year)[1]), "zeta_s")
plot_map_raster(p$data, "est") + ggtitle("Predicted density (g/sq.m) (fixed effects + all random effects)")
plot_map_raster(p$data, "exp(est_non_rf)") + ggtitle("Prediction (fixed effects only)")
plot_map_raster(p$data, "omega_s") + ggtitle("Spatial random effects only")
plot_map_raster(p$data, "epsilon_st") + ggtitle("Spatiotemporal random effects only")

# look at just the spatiotemporal random effects:
plot_map_raster(p$data, "est_rf") + scale_fill_gradient2()

(density_map = p$data %>% 
    ggplot(aes_string("X", "Y", fill = "exp(est)")) +
    geom_tile(aes(height = 0.1, width = 0.1)) +
    facet_wrap(~year) +
    coord_fixed() +
    xlab("Eastings (km)") +
    ylab("Northings (km)") + 
    scale_fill_gradientn(colours = matlab.like(100), "g/sq.m") +
    ggtitle(paste0("Predicted ", response, " 2015-2019")) +
    ggdark::dark_theme_minimal() +
    theme(legend.position = "right"))

index <- get_index(p, bias_correct = F)

ggdark::invert_geom_defaults()

(relative_biomass = index %>%
    ggplot(aes(year, est)) + 
    geom_line() +
    geom_point(size = 3) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, colour = NA) +
    xlab('Year') + 
    ylab('biomass') +
    ggtitle("Biomass estimate") +
    theme_pubr())

index %>% 
  mutate(cv = sqrt(exp(se^2) - 1)) %>% 
  dplyr::select(-log_est, -max_gradient, -bad_eig, -se) %>%
  knitr::kable(format = "pandoc", digits = c(0, 0, 0, 0, 2))

# Calculate centre of gravity for latitude and longitude

cog <- get_cog(p) # calculate centre of gravity for each data point
cog$utm = ifelse(cog$coord == "X", "Eastings", "Northings")

(density_cog = ggplot(cog, aes(year, est, ymin = lwr, ymax = upr)) +
    geom_ribbon(alpha = 0.2) +
    geom_line() + 
    geom_point(size = 3) +
    facet_wrap(~utm, scales = "free_y") + 
    ggtitle("Center of gravity") + 
    theme_pubr())
