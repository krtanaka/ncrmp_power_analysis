library(sdmTMB)
library(dplyr)
library(ggplot2)
library(rgdal)
library(colorRamps)

rm(list = ls())

# 4 functional groups
# live coral cover
# adult juvenile density

# for Uku: 
# Total numerical density estimates (individuals per 100 m2) were obtained by dividing fish counts in each survey by the survey area (353 m2 from two 15-m diameter survey cylinders) and multiplying by 100. - Nadon et al. 2020

region = "MHI"
uku_or_not = F

# load("data/ALL_REA_FISH_RAW.rdata")
# df %>% 
#   subset(REGION == region) %>% 
#   group_by(TAXONNAME) %>%
#   summarise(n = sum(BIOMASS_G_M2, na.rm = T)) %>%
#   # summarise(n = sum(COUNT, na.rm = T)) %>% 
#   mutate(freq = n/sum(n)) %>% 
#   arrange(desc(freq)) %>% 
#   top_n(5) 

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

response_variable = "fish_count";      sp = ifelse(uku_or_not == T, "Aprion virescens", "Chromis vanderbilti")
response_variable = "fish_biomass";    sp = ifelse(uku_or_not == T, "Aprion virescens", "Acanthurus olivaceus")
response_variable = "trophic_biomass"; sp = c("PISCIVORE", "PLANKTIVORE", "PRIMARY", "SECONDARY")[2]
response_variable = "coral_cover";     sp = c("CCA", "CORAL", "EMA", "HAL", "I", "MA", "SC", "SED", "TURF")[2]
response_variable = "coral_density";   sp = c("AdColDen", "JuvColDen")[1]

if (response_variable == "fish_count") {
  
  load("data/ALL_REA_FISH_RAW.rdata")
  
  df = df %>% 
    subset(REGION == region & ISLAND %in% islands) %>% 
    mutate(response = ifelse(TAXONNAME == sp, COUNT*100, 0)) %>%  
    group_by(LONGITUDE, LATITUDE, ISLAND, OBS_YEAR, DATE_, DEPTH) %>% 
    summarise(response = sum(response, na.rm = T))
  
  df %>% ggplot(aes(response)) + geom_histogram() + 
    df %>% group_by(OBS_YEAR) %>% summarise(n = median(response)) %>% ggplot(aes(OBS_YEAR, n)) + geom_line()
  
}

if (response_variable == "fish_biomass") {
  
  load("data/ALL_REA_FISH_RAW.rdata")
  
  df = df %>% 
    subset(REGION == region & ISLAND %in% islands) %>% 
    mutate(response = ifelse(TAXONNAME == sp, BIOMASS_G_M2, 0)) %>%  
    # mutate(response = ifelse(TAXONNAME == sp, BIOMASS_G_M2*0.001, 0)) %>%  
    group_by(LONGITUDE, LATITUDE, ISLAND, OBS_YEAR, DATE_, DEPTH) %>% 
    summarise(response = sum(response, na.rm = T))  
  
  df %>% ggplot(aes(response)) + geom_histogram() + 
    df %>% group_by(OBS_YEAR) %>% summarise(n = median(response)) %>% ggplot(aes(OBS_YEAR, n)) + geom_line()
  
} 

if (response_variable == "trophic_biomass") {
  
  load("data/ALL_REA_FISH_RAW.rdata")
  
  df = df %>% 
    subset(REGION == region & ISLAND %in% islands) %>% 
    mutate(response = ifelse(TROPHIC_MONREP == sp, BIOMASS_G_M2, 0)) %>%  
    group_by(LONGITUDE, LATITUDE, ISLAND, OBS_YEAR, DATE_, DEPTH) %>% 
    summarise(response = sum(response, na.rm = T))  
  
  df %>% ggplot(aes(response)) + geom_histogram() + 
    df %>% group_by(OBS_YEAR) %>% summarise(n = median(response)) %>% ggplot(aes(OBS_YEAR, n)) + geom_line()
  
}

if (response_variable == "coral_cover") {
  
  load("data/BenthicCover_2010-2020_Tier1_SITE_MHI_w_CRM.RData") #live coral cover, only for MHI, with CRM_Bathy data
  
  df = df %>% 
    subset(REGION == region & ISLAND %in% islands) %>% 
    mutate(response = CORAL,
           DEPTH = ifelse(DEPTH_e == 0, DEPTH_e*-1 + 0.1, DEPTH_e*-1)) %>%  
    group_by(LONGITUDE, LATITUDE, ISLAND, OBS_YEAR, DATE_, DEPTH) %>% 
    summarise(response = median(response, na.rm = T),
              n = n())
  
  hist(df$response, main = paste0(sp, "_cover"))
  
}

if (response_variable == "coral_density") {
  
  load("data/BenthicREA_sitedata_TAXONCODE.RData_MHI_w_CRM.RData") #live coral cover, only for MHI, with CRM_Bathy data
  
  if (sp == "AdColDen") df = df %>% mutate(response = AdColDen)
  if (sp == "JuvColDen") df = df %>% mutate(response = JuvColDen)
  
  df = df %>% 
    subset(REGION == region & ISLAND %in% islands) %>% 
    mutate( DEPTH = ifelse(DEPTH_e == 0, DEPTH_e*-1 + 0.1, DEPTH_e*-1)) %>%  
    group_by(LONGITUDE, LATITUDE, ISLAND, OBS_YEAR, DATE_, DEPTH) %>% 
    summarise(response = sum(response, na.rm = T))
  
  hist(df$response, main = paste0(sp, "_coral_density"))
  
}

# north-south gradient
df %>% 
  group_by(ISLAND) %>% 
  summarise(n = mean(response, na.rm = T),
            lat = mean(LATITUDE)) %>% 
  arrange(desc(lat)) 

zone <- (floor((df$LONGITUDE[1] + 180)/6) %% 60) + 1
xy_utm = as.data.frame(cbind(utm = project(as.matrix(df[, c("LONGITUDE", "LATITUDE")]), paste0("+proj=utm +units=km +zone=", zone))))
colnames(xy_utm) = c("X", "Y"); plot(xy_utm, pch = ".", bty = 'n')
df = cbind(df, xy_utm)

n_knots = 300
# n_knots = 100 # a coarse mesh for speed
rea_spde <- make_mesh(df, c("X", "Y"), n_knots = n_knots, type = "cutoff_search") 

plot(rea_spde, pch = "."); axis(1); axis(2)

df$year = df$OBS_YEAR
df$depth = df$DEPTH
df$depth_scaled = scale(log(df$depth))
df$depth_scaled2 = df$depth_scaled ^ 2

plot(df$depth, df$response, pch = 20, bty = "n")
plot(df[11:13], pch = ".")

obs_year = unique(df$year)
full_year = seq(min(df$year), max(df$year), by = 1)
missing_year = setdiff(full_year, obs_year)
missing_year = as.integer(missing_year);missing_year

density_model <- sdmTMB(
  
  data = df, 
  
  # formula = response ~ as.factor(year) + depth_scaled + depth_scaled2,
  formula = response ~ as.factor(year) + s(depth, k=3),
  
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
  
  control = sdmTMBcontrol(step.min = 0.01, step.max = 1)
  
)

# look at gradients
max(density_model$gradients)

df$residuals <- residuals(density_model)
qqnorm(df$residuals, ylim = c(-5, 5), xlim = c(-5, 5), bty = "n", pch = 20);abline(a = 0, b = 1)

m_p <- predict(density_model); m_p = m_p[,c("response", "est")]

ggplot(df, aes_string("X", "Y", fill = "residuals")) +
  geom_tile(aes(height = 0.5, width = 0.5)) +
  facet_wrap(.~ISLAND, scales = "free") +
  xlab("Eastings") +
  ylab("Northings") + 
  scale_fill_gradient2() + 
  ggdark::dark_theme_void()

ggdark::invert_geom_defaults()

m_p  %>% 
  ggplot(aes(response, exp(est))) + 
  geom_point(alpha = 0.2) + 
  coord_fixed(ratio = 1) +
  ylab("predicted_density") + 
  xlab("observed_density") + 
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method = "lm", se = T)

#  extract some parameter estimates
sd <- as.data.frame(summary(TMB::sdreport(density_model$tmb_obj)))
r <- density_model$tmb_obj$report()
r

# prediction onto new data grid
load("data/Topography_NOAA_CRM_vol10.RData")

# topo <- topo %>% subset(x < -157.5 & x > -158.5 & y > 21 & y < 22) #oahu
# topo <- topo %>% subset(x < -154.8 & x > -156.2 & y > 18.8 & y < 20.4) #hawaii

grid = topo

# res = 2
# 
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
                                           paste0("+proj=utm +units=km +zone=", zone))))

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

# set the area argument to 0.0081 km2 since our grid cells are 90 m x 90 m = 0.0081 square kilometers
p <- predict(density_model, 
             newdata = grid_year, 
             return_tmb_object = T, area = 0.0081)

p$data$sp = sp
sdm_output = p$data

save(sdm_output, file = paste0("outputs/density_results_", sp, "_", response_variable, "_", n_knots, "_", region, ".RData"))

plot_map_raster <- function(dat, column = "est") {
  
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_tile(aes(height = 0.5, width = 0.5)) +
    facet_wrap(~year) +
    coord_fixed() +
    xlab("Eastings (km)") +
    ylab("Northings (km)") + 
    scale_fill_gradientn(colours = matlab.like(100), "") +
    ggdark::dark_theme_minimal()
  
}

# pick out a single year to plot since they should all be the same for the slopes. Note that these are in log space.
plot_map_raster(filter(p$data, year == 2015), "zeta_s")
plot_map_raster(p$data, "exp(est)") + ggtitle("Predicted density (g/sq.m) (fixed effects + all random effects)")
plot_map_raster(p$data, "exp(est_non_rf)") + ggtitle("Prediction (fixed effects only)")
plot_map_raster(p$data, "omega_s") + ggtitle("Spatial random effects only")
plot_map_raster(p$data, "epsilon_st") + ggtitle("Spatiotemporal random effects only")

# look at just the spatiotemporal random effects:
plot_map_raster(p$data, "est_rf") + scale_fill_gradient2()

density_map = ggplot(p$data, aes_string("X", "Y", fill = "exp(est)")) +
  geom_tile(aes(height = 0.5, width = 0.5)) +
  facet_wrap(~year) +
  coord_fixed() +
  xlab("Eastings (km)") +
  ylab("Northings (km)") + 
  scale_fill_gradientn(colours = matlab.like(100), "g/m^2") +
  ggtitle("Uku predicted density (fixed effects + random effects)") + 
  ggdark::dark_theme_minimal() + 
  theme(legend.position = "right")

index <- get_index(p, bias_correct = F)

ggdark::invert_geom_defaults()

relative_biomass = index %>%
  ggplot(aes(year, est)) + 
  geom_line() +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, colour = NA) +
  xlab('Year') + 
  # ylab('metric tonnes') + 
  ylab('total count (n)') + 
  ggtitle("Biomass estimate") + 
  ggdark::dark_theme_minimal()
# theme_pubr()

index %>% 
  mutate(cv = sqrt(exp(se^2) - 1)) %>% 
  dplyr::select(-log_est, -max_gradient, -bad_eig, -se) %>%
  knitr::kable(format = "pandoc", digits = c(0, 0, 0, 0, 2))

# Calculate centre of gravity for latitude and longitude

cog <- get_cog(p) # calculate centre of gravity for each data point
cog$utm = ifelse(cog$coord == "X", "Eastings", "Northings")

density_cog = ggplot(cog, aes(year, est, ymin = lwr, ymax = upr)) +
  geom_ribbon(alpha = 0.2) +
  geom_line() + 
  geom_point(size = 3) +
  facet_wrap(~utm, scales = "free_y") + 
  ggtitle("Center of gravity") + 
  ggdark::dark_theme_minimal()
# theme_pubr()

# table of COG by latitude
plot(data.frame(Y = p$data$Y, est = exp(p$data$est), year = p$data$year) %>%
       group_by(year) %>% summarize(cog = sum(Y * est) / sum(est)), type = "b", bty = "l", ylab = "Northing")

library(patchwork)

density_map
relative_biomass+density_cog
