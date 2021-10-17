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

if (unit == "biomass") response_variables = c("PISCIVORE_BIO", "PLANKTIVORE_BIO", "PRIMARY_BIO", "SECONDARY_BIO", "TotFishBio")
if (unit == "abundance") response_variables = c("PISCIVORE_ABUN", "PLANKTIVORE_ABUN", "PRIMARY_ABUN", "SECONDARY_ABUN", "TotFishAbund")

knots = c(100, 300, 500, 1000)[3]

for (r in 1:length(response_variables)) {
  
  # r = 5
  
  response = response_variables[r]
  
  wsd$time = as.numeric(as.POSIXct(wsd$DATE_, format="%Y-%m-%d")) # unix time
  
  if (response == "PISCIVORE_BIO") wsd$response = wsd$PISCIVORE_BIO
  if (response == "PLANKTIVORE_BIO") wsd$response = wsd$PLANKTIVORE_BIO
  if (response == "PRIMARY_BIO") wsd$response = wsd$PRIMARY_BIO
  if (response == "SECONDARY_BIO") wsd$response = wsd$SECONDARY_BIO
  if (response == "TotFishBio") wsd$response = wsd$TotFishBio
  
  df = wsd %>% 
    subset(ISLAND %in% islands) %>%
    # group_by(LONGITUDE, LATITUDE, time) %>%
    group_by(LONGITUDE, LATITUDE, OBS_YEAR) %>%
    summarise(response = mean(response, na.rm = T),
              depth = mean(DEPTH, na.rm = T)) %>% 
    subset(response < quantile(response, prob = 0.999))
  
  df %>% ggplot(aes(response)) + geom_histogram() +
    df %>% group_by(OBS_YEAR) %>% summarise(n = mean(response)) %>% ggplot(aes(OBS_YEAR, n)) + geom_point() + geom_line()
  
  zone <- (floor((df$LONGITUDE[1] + 180)/6) %% 60) + 1
  xy_utm = as.data.frame(cbind(utm = project(as.matrix(df[, c("LONGITUDE", "LATITUDE")]), paste0("+proj=utm +units=km +zone=", zone))))
  colnames(xy_utm) = c("X", "Y")
  df = cbind(df, xy_utm)
  
  # Read in coastline boundaries
  load('data/misc/MHI_islands_shp.RData')
  crs(ISL_bounds) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  ISL_this = ISL_bounds[which(ISL_bounds$ISLAND %in% toupper(islands)),]
  ISL_this_utm = spTransform(ISL_this,CRS(paste0("+proj=utm +units=km +zone=",zone)))
  ISL_this_sf = st_transform(st_as_sf(ISL_this), crs = paste0("+proj=utm +units=km +zone=",zone))
  
  rea_spde <- make_mesh(df, c("X", "Y"), n_knots = knots, type = "cutoff_search") # search
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
    
    # formula = response ~ 0 + as.factor(year), #index standardization
    # formula = response ~ 1, #just estimates and intercept and accounts for all other variation through the random effects
    # formula = response ~ as.factor(year) + depth_scaled + depth_scaled2,
    formula = response ~ as.factor(year) + s(depth, k = 5),
    # formula = response ~ as.factor(year) + s(depth, k = 3),
    
    silent = F, 
    spatial_trend = T, 
    spatial_only = F, 
    time = "year", 
    spde = rea_spde, 
    anisotropy = T,
    # family = tweedie(link = "log"),
    # family = poisson(link = "log"),
    family = nbinom2(link = "log"),
    # family = Beta(link = "logit"),
    
    control = sdmTMBcontrol(step.min = 0.01, step.max = 1)
    
  ); beepr::beep(2)
  
  density_model <- run_extra_optimization(density_model); beepr::beep(2)
  
  # look at gradients
  max(density_model$gradients)
  
  df$residuals <- residuals(density_model)
  
  hist(df$residuals, breaks = 30)
  
  par(pty = "s")
  
  png(paste0('outputs/qq_', response, '_', knots, '_knots.png'), height = 4, width = 4, units = "in", res = 100)
  qqnorm(df$residuals, ylim = c(-5, 5), xlim = c(-5, 5), bty = "n", pch = 20, main = response); abline(a = 0, b = 1)
  dev.off()
  
  m_p <- predict(density_model); m_p = m_p[,c("response", "est")]
  m_p$back_abs_res = abs(df$residuals)
  
  p1 = ggplot(df, aes_string("X", "Y", color = "residuals")) +
    geom_point(alpha = 0.8, size = round(abs(df$residuals), digits = 0)) + 
    xlab("Eastings") +
    ylab("Northings") + 
    scale_color_gradient2() + 
    theme_classic()
  
  p2 = m_p  %>% 
    ggplot(aes(response, exp(est))) + 
    geom_point(alpha = 0.2, aes(size = back_abs_res)) + 
    coord_fixed(ratio = 1) +
    ylab("prediction") + 
    xlab("observation") + 
    geom_abline(intercept = 0, slope = 1) +
    geom_smooth(method = "lm", se = T) + 
    theme_classic()
  
  png(paste0('outputs/residuals_', response, '_', knots, '_knots.png'), height = 5, width = 5, units = "in", res = 100)
  print(p1 / p2)
  dev.off()
  
  # #  extract some parameter estimates
  # sd <- as.data.frame(summary(TMB::sdreport(density_model$tmb_obj)))
  # r <- density_model$tmb_obj$report()
  # r
  
  # prediction onto new data grid
  load(paste0("data/gis_bathymetry/raster/", substr(tolower(islands), 1, 3), ".RData")) # bathymetry 
  
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
  
  family = as.character(density_model$family[1])
  
  save(sdm_output, file = paste0("outputs/sdmTMB_results_", response, "_", unit, "_", knots, "knots_", family, ".RData"))
  
  # plot_map_raster <- function(dat, column = "est") {
  # 
  #   ggplot(dat, aes_string("X", "Y", fill = column)) +
  #     geom_tile(aes(height = 0.5, width = 0.5), alpha = 0.5) +
  #     # geom_raster() +
  #     facet_wrap(~year) +
  #     coord_fixed() +
  #     xlab("Eastings (km)") +
  #     ylab("Northings (km)") +
  #     scale_fill_gradientn(colours = matlab.like(100), "") +
  #     ggdark::dark_theme_minimal()
  # 
  # }
  # 
  # # pick out a single year to plot since they should all be the same for the slopes. Note that these are in log space.
  # plot_map_raster(filter(p$data, year == unique(p$data$year)[1]), "zeta_s")
  # plot_map_raster(p$data, "est") + ggtitle("Predicted density (g/sq.m) (fixed effects + all random effects)")
  # plot_map_raster(p$data, "exp(est_non_rf)") + ggtitle("Prediction (fixed effects only)")
  # plot_map_raster(p$data, "omega_s") + ggtitle("Spatial random effects only")
  # plot_map_raster(p$data, "epsilon_st") + ggtitle("Spatiotemporal random effects only")
  # 
  # look at just the spatiotemporal random effects:
  # plot_map_raster(p$data, "est_rf") + scale_fill_gradient2()
  
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
      # theme_pubr() + 
      theme(legend.position = "right"))
  
  ggdark::invert_geom_defaults()
  
  png(paste0('outputs/predicted_', unit, '_', response, '_', knots, '_knots.png'), height = 5, width = 5, units = "in", res = 100)
  print(density_map)
  dev.off()
  
  # index <- get_index(p, bias_correct = F)
  # 
  # ggdark::invert_geom_defaults()
  # 
  # (relative_biomass = index %>%
  #     ggplot(aes(year, est)) +
  #     geom_line() +
  #     geom_point(size = 3) +
  #     geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, colour = NA) +
  #     xlab('Year') +
  #     ylab('biomass') +
  #     ggtitle("Biomass estimate") +
  #     theme_pubr())
  # 
  # index %>% 
  #   mutate(cv = sqrt(exp(se^2) - 1)) %>% 
  #   dplyr::select(-log_est, -max_gradient, -bad_eig, -se) %>%
  #   knitr::kable(format = "pandoc", digits = c(0, 0, 0, 0, 2))
  # 
  # # Calculate centre of gravity for latitude and longitude
  # 
  # cog <- get_cog(p) # calculate centre of gravity for each data point
  # cog$utm = ifelse(cog$coord == "X", "Eastings", "Northings")
  # 
  # (density_cog = ggplot(cog, aes(year, est, ymin = lwr, ymax = upr)) +
  #     geom_ribbon(alpha = 0.2) +
  #     geom_line() + 
  #     geom_point(size = 3) +
  #     facet_wrap(~utm, scales = "free_y") + 
  #     ggtitle("Center of gravity") + 
  #     theme_pubr())
  
}