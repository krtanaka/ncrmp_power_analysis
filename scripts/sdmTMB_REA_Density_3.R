library(sdmTMB)
library(dplyr)
library(ggplot2)
library(rgdal)
library(colorRamps)

rm(list = ls())

load("data/ALL_REA_FISH_RAW.rdata")

df = df %>% 
  subset(REGION == "MHI") %>% 
  mutate(density = COUNT)

df %>% 
  group_by(TAXONNAME) %>% 
  # summarise(n = sum(BIOMASS_G_M2, na.rm = T)) %>% 
  summarise(n = sum(density, na.rm = T)) %>% 
  mutate(freq = n/sum(n)) %>% 
  arrange(desc(freq))

# df$density = ifelse(df$TAXONNAME == "Aprion virescens", df$density, 0)
df$density = ifelse(df$TAXONNAME == "Chromis vanderbilti", df$density, 0) # most abundant in MHI

df %>% 
  group_by(ISLAND) %>% 
  summarise(n = mean(density, na.rm = T),
            lat = mean(LATITUDE)) %>% 
  arrange(desc(lat))

islands = c("Kauai", #1
            "Lehua", #2
            "Niihau", #3
            "Kaula", #4
            "Oahu", #5
            "Molokai", #6
            "Maui", #7
            "Lanai", #8
            "Molokini", #9
            "Kahoolawe", #10
            "Hawaii")[11]

df = df %>% 
  subset(ISLAND %in% islands) %>% 
  group_by(LONGITUDE, LATITUDE, OBS_YEAR, DEPTH) %>% 
  summarise(density = mean(density, na.rm = T))

hist(df$density)
summary(df$density)

zone <- (floor((df$LONGITUDE[1] + 180)/6) %% 60) + 1
xy_utm = as.data.frame(cbind(utm = project(as.matrix(df[, c("LONGITUDE", "LATITUDE")]), paste0("+proj=utm +zone=", zone))))
colnames(xy_utm) = c("X", "Y"); plot(xy_utm, pch = ".", bty = 'n')
df = cbind(df, xy_utm)

rea_spde <- make_mesh(df, c("X", "Y"), n_knots = 300, type = "cutoff_search") # a coarse mesh for speed

plot(rea_spde, pch = "."); axis(1); axis(2)

df$year = df$OBS_YEAR
df$depth = df$DEPTH
df$depth_scaled = scale(log(df$depth))
df$depth_scaled2 = df$depth_scaled ^ 2

qplot(df$depth, df$density, color = df$density)
plot(df[9:11])

obs_year = unique(df$year)
full_year = seq(min(df$year), max(df$year), by = 1)
missing_year = setdiff(full_year, obs_year)
missing_year = as.integer(missing_year)

density_model <- sdmTMB(
  
  data = df, 
  formula = density ~ as.factor(year) + depth_scaled + depth_scaled2,
  silent = F, 
  # extra_time = missing_year,
  spatial_trend = T, 
  spatial_only = F, 
  time = "year", 
  spde = rea_spde, 
  anisotropy = T,
  family = tweedie(link = "log"),
  control = sdmTMBcontrol(step.min = 0.01, step.max = 1)
  
)

# look at gradients
max(density_model$gradients)

df$residuals <- residuals(density_model)
qqnorm(df$residuals, ylim = c(-5, 5), xlim = c(-5, 5), bty = "n");abline(a = 0, b = 1)

m_p <- predict(density_model); m_p = m_p[,c("density", "est")]

m_p  %>% 
  ggplot(aes(density, exp(est))) + 
  geom_point(alpha = 0.2) + 
  coord_fixed() +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method = "lm", se = F) + 
  ggdark::dark_theme_minimal()

plot_map_point <- function(dat, column = "est") {
  ggplot(dat, aes_string("X", "Y", color = column)) +
    geom_point(alpha = 0.2, size = 2) +
    coord_fixed() +
    xlab("Eastings") +
    ylab("Northings") #+ ggdark::dark_theme_light()
}

plot_map_point(df, "residuals") + scale_color_gradient2()

#  extract some parameter estimates
sd <- as.data.frame(summary(TMB::sdreport(density_model$tmb_obj)))
r <- density_model$tmb_obj$report()

# prediction onto new data grid
load("data/Topography_NOAA_CRM_vol10.RData")

# topo <- topo %>% subset(x < -157.5 & x > -158.5 & y > 21 & y < 22) #oahu
topo <- topo %>% subset(x < -154.8 & x > -156.2 & y > 18.8 & y < 20.4) #hawaii

grid = topo

res = 5

# grid$longitude = round(grid$x, digits = res)
# grid$latitude = round(grid$y, digits = res)

grid$longitude = grid$x
grid$latitude = grid$y

grid = grid %>% 
  group_by(longitude, latitude) %>% 
  subset(longitude > range(df$LONGITUDE)[1]) %>% 
  subset(longitude < range(df$LONGITUDE)[2]) %>%  
  subset(latitude > range(df$LATITUDE)[1]) %>%
  subset(latitude < range(df$LATITUDE)[2]) %>% 
  summarise(depth = mean(Topography, na.rm = T)*-1) 

zone <- (floor((grid$longitude[1] + 180)/6) %% 60) + 1
xy_utm = as.data.frame(cbind(utm = project(as.matrix(grid[, c("longitude", "latitude")]),
                                           paste0("+proj=utm +zone=", zone))))

colnames(xy_utm) = c("X", "Y"); plot(xy_utm, pch = ".")

grid = cbind(grid, xy_utm)

grid_year = NULL

year = as.vector(unique(df$year))

for (y in 1:length(year)) {
  
  # y = 1
  
  grid_y = grid  
  grid_y$year = year[[y]]  
  
  grid_year = rbind(grid_year, grid_y)
  
}

grid_year$depth_scaled = scale(log(grid_year$depth+0.0001))
grid_year$depth_scaled2 = grid_year$depth_scaled ^ 2

grid_year_missing = NULL

for (y in 1:length(missing_year)) {
  
  # y = 1
  
  # Add missing years to our grid:
  grid_from_missing_yr <- grid_year[grid_year$year ==  missing_year[[y]]-1, ]
  grid_from_missing_yr$year <-  missing_year[[y]] # `L` because `year` is an integer in the data
  grid_year_missing <- rbind(grid_year_missing, grid_from_missing_yr)
  
}

# grid_year = rbind(grid_year, grid_year_missing)

p <- predict(density_model, 
             newdata = grid_year, 
             return_tmb_object = T) # area = 2.5)

plot_map_raster <- function(dat, column = "est") {
  
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_tile(aes(height = 500, width = 500)) +
    # geom_point() +
    facet_wrap(~year) +
    coord_fixed() +
    xlab("Eastings") +
    ylab("Northings") + 
    # scale_fill_viridis_c() +
    scale_fill_gradientn(colours = matlab.like(100)) #+ ggdark::dark_theme_void()
  
}

# pick out a single year to plot since they should all be the same for the slopes. Note that these are in log space.
plot_map_raster(filter(p$data, year == 2015), "zeta_s")
plot_map_raster(p$data, "exp(est)") + ggtitle("Predicted density (g/sq.m) (fixed effects + all random effects)")
plot_map_raster(p$data, "exp(est_non_rf)") + ggtitle("Prediction (fixed effects only)")
plot_map_raster(p$data, "omega_s") + ggtitle("Spatial random effects only")
plot_map_raster(p$data, "epsilon_st") + ggtitle("Spatiotemporal random effects only")

# look at just the spatiotemporal random effects:
plot_map_raster(p$data, "est_rf") + scale_fill_gradient2()

density_map = ggplot(p$data, aes_string("X", "Y", fill = "est", color = "est")) +
  geom_tile(aes(height = 500, width = 500)) +
  # geom_point() +
  facet_wrap(~year) +
  coord_fixed() +
  xlab("Eastings") +
  ylab("Northings") + 
  scale_fill_viridis_c("log(g/sq.m)") + 
  scale_color_viridis_c("log(g/sq.m)") + 
  ggtitle("Predicted density (g/sq.m) (fixed effects + all random effects)") + 
  theme_minimal()

index <- get_index(p, bias_correct = T)

relative_biomass = index %>%
  ggplot(aes(year, est)) + 
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, colour = NA) +
  xlab('Year') + 
  ylab('metric tonnes') + 
  ggtitle("Biomass estimate") + 
  theme_minimal()

index %>% 
  mutate(cv = sqrt(exp(se^2) - 1)) %>% 
  dplyr::select(-log_est, -max_gradient, -bad_eig, -se) %>%
  knitr::kable(format = "pandoc", digits = c(0, 0, 0, 0, 2))

# Calculate centre of gravity for latitude and longitude

cog <- get_cog(p) # calculate centre of gravity for each data point

density_cog = ggplot(cog, aes(year, est, ymin = lwr, ymax = upr)) +
  geom_ribbon(alpha = 0.2) +
  geom_line() + 
  geom_point() + 
  facet_wrap(~coord, scales = "free_y") + 
  ggtitle("center of gravity (lat and lon)") + 
  theme_minimal()

# table of COG by latitude
plot(data.frame(Y = p$data$Y, est = exp(p$data$est), year = p$data$year) %>%
  group_by(year) %>% summarize(cog = sum(Y * est) / sum(est)), type = "b")

library(patchwork)

density_map + (relative_biomass / density_cog)

