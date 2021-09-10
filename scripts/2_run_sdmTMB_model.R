library(sdmTMB)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(rgdal)
library(sp)
library(raster)
library(colorRamps)
library(patchwork)
library(fitdistrplus)
library(sf)
library(lubridate)
library(boot)

rm(list = ls())

SmithsonVerkuilen2006 = function(y){
  n = length(y)
  out = (y * (n-1) + 0.5) / n
  return(out)
}

region = "MHI"

load("data/rea/SURVEY MASTER.RData")

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
            "Hawaii")[5]

response_variable = "coral_cover";     sp = c("CCA", "CORAL", "EMA", "HAL", "I", "MA", "SC", "SED", "TURF")[2]
response_variable = "coral_density";   sp = c("AdColDen", "JuvColDen")[1]

LOAD_EDS = F
if(LOAD_EDS) {
  EDS = read.csv("M:/Environmental Data Summary/Outputs/Survey_Master_Timeseries_2021-02-27.csv")
} 

if (response_variable == "coral_cover") {
  
  load("data/rea/BenthicCover_2010-2020_Tier1_SITE_MHI_w_CRM.RData") #live coral cover, only for MHI, with CRM_Bathy data
  df$LONGITUDE360=df$LONGITUDE
  df$LONGITUDE=df$LONGITUDE360-180
  df$MAX_DEPTH_M=SURVEY_MASTER$new_MAX_DEPTH_M[match(df$SITEVISITID,SURVEY_MASTER$SITEVISITID)]
  df$MIN_DEPTH_M=SURVEY_MASTER$new_MIN_DEPTH_M[match(df$SITEVISITID,SURVEY_MASTER$SITEVISITID)]
  df$MN_DEPTH_M=apply(df[,c("MIN_DEPTH_M","MAX_DEPTH_M")],1,mean,na.rm=T)
  df = df %>% 
    subset(REGION == region & ISLAND %in% islands) %>% 
    mutate(response = CORAL,
           DEPTH = ifelse(DEPTH_e == 0, DEPTH_e*-1 + 0.1, DEPTH_e*-1)) %>%   
    group_by(SITEVISITID,LONGITUDE, LATITUDE, ISLAND, OBS_YEAR, DATE_, DEPTH) %>% 
    summarise(response = median(response, na.rm = T),
              n = n(),
              mn_depth_m = mean(MN_DEPTH_M,na.rm=T)) %>% 
    subset(!is.na(mn_depth_m))
  df$response_beta = SmithsonVerkuilen2006(df$response/100)
  
  if(LOAD_EDS){
    df=left_join(df,
                 EDS[,
                     c(which(names(EDS)=="SITEVISITID"),
                       which(names(EDS)=="DHW.Np10y_Degree_Heating_Weeks_MO03"):length(names(EDS)))
                 ],
                 by="SITEVISITID",)
  }
  # hist(df$response, main = paste0(sp, "_cover"))
  # descdist(df$response)
}

if (response_variable == "coral_density") {
  
  load("data/rea/BenthicREA_sitedata_TAXONCODE.RData_MHI_w_CRM.RData") #live coral cover, only for MHI, with CRM_Bathy data
  
  if (sp == "AdColDen") df = df %>% mutate(response = AdColDen)
  if (sp == "JuvColDen") df = df %>% mutate(response = JuvColDen)
  
  df = df %>% 
    subset(REGION == region & ISLAND %in% islands) %>% 
    mutate( DEPTH = ifelse(DEPTH_e == 0, DEPTH_e*-1 + 0.1, DEPTH_e*-1),  
            MN_DEPTH_M = (MIN_DEPTH_M + MAX_DEPTH_M)/2) %>%  
    group_by(LONGITUDE, LATITUDE, OBS_YEAR) %>% 
    summarise(response = mean(response, na.rm = T), 
              depth = mean(DEPTH, na.rm = T),
              mn_depth_m = mean(MN_DEPTH_M, na.rm = T))
  
  hist(df$response, main = paste0(sp, "_coral_density"),30)
  
}

# north-south gradient
# df %>% 
#   group_by(ISLAND) %>% 
#   summarise(n = mean(response, na.rm = T),
#             lat = mean(LATITUDE)) %>% 
#   arrange(desc(lat)) 

zone <- (floor((df$LONGITUDE[1] + 180)/6) %% 60) + 1
xy_utm = as.data.frame(cbind(utm = project(as.matrix(df[, c("LONGITUDE", "LATITUDE")]), paste0("+proj=utm +units=km +zone=", zone))))
colnames(xy_utm) = c("X", "Y")
df = cbind(df, xy_utm)

#Read in Island Boundaries
#"T:/Common/Maps/Island/islands"
ISL_bounds = readOGR(dsn = "T:/Common/Maps/Island", layer = "islands")
crs(ISL_bounds) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

ISL_this = ISL_bounds[which(ISL_bounds$ISLAND %in% toupper(islands)),]
ISL_this_utm = spTransform(ISL_this,CRS(paste0("+proj=utm +units=km +zone=",zone)))
ISL_this_sf = st_transform(st_as_sf(ISL_this), crs = paste0("+proj=utm +units=km +zone=",zone))

#n_knots = 300 
rea_spde <- make_mesh(df, c("X", "Y"), n_knots = n_knots, type = "cutoff_search")
rea_spde <- make_mesh(df, c("X", "Y"), cutoff  = n_knots, type = "cutoff") # predefined

dists = 2#c(5,1,.25)
PLOT_ALONG = T

for (disti in 1:length(dists)){
  
  # rea_spde <- make_mesh(df, c("X", "Y"), cutoff = dists[disti], type = "cutoff") 
  n_knots = 100
  rea_spde <- make_mesh(df, c("X", "Y"), n_knots = n_knots, type = "cutoff_search")
  (rea_spde$mesh$n)
  
  #build barrier to mesh
  rea_spde_coast = add_barrier_mesh(rea_spde,ISL_this_sf)
  
  #Build potential covariates
  df$year = df$OBS_YEAR
  df$depth_scaled = scale(log(df$DEPTH))
  df$depth_scaled2 = df$depth_scaled ^ 2
  
  if(PLOT_ALONG){
    (rea_spde_coast$mesh$n)
    dev.off()
    #par(mfrow = c(2,2))
    #plot(xy_utm, pch = ".", bty = 'n')
    
    plot(rea_spde_coast$mesh,asp=1); axis(1); axis(2)
    plot(ISL_this_utm,add=TRUE)
    points(rea_spde_coast$loc_xy,col="blue",pch=20,cex=1)
    bar_i=rea_spde_coast$barrier_triangles
    norm_i=rea_spde_coast$normal_triangles
    points(rea_spde_coast$spde$mesh$loc[,1],rea_spde_coast$spde$mesh$loc[,2],pch=20,col="black",cex=1.5)
    points(rea_spde_coast$mesh_sf$V1[bar_i],rea_spde_coast$mesh_sf$V2[bar_i],col="red",pch=4,cex=.5)
    points(rea_spde_coast$mesh_sf$V1[norm_i],rea_spde_coast$mesh_sf$V2[norm_i],col="green",pch=1,cex=.5)
    
    plot(df$mn_depth_m, df$response_beta, pch = 20, bty = "n")
    #plot(df$mean_SST_CRW_Daily_YR03, df$response, pch = 20, bty = "n")
  }
  
  obs_year = unique(df$year)
  full_year = seq(min(df$year), max(df$year), by = 1)
  missing_year = setdiff(full_year, obs_year)
  missing_year = as.integer(missing_year);missing_year
  
  density_model <- sdmTMB(
    
    data = df, 
    
    # formula = response ~ 1,
    formula = response_beta ~ as.factor(year) + s(mn_depth_m,k=5),
    # formula = response ~ as.factor(year) + depth_scaled + depth_scaled2,
    # formula = response ~ as.factor(year) + s(depth, k=3),
    # formula = response ~ as.factor(year) + s(depth, k=5) + s(temp, k=5),
    # formula = response ~ as.factor(year) + s(temp, k=5) + s(depth, k=5) + depth_scaled + depth_scaled2,
    silent = F, 
    # extra_time = missing_year,
    spatial_trend = T, 
    spatial_only = F, 
    time = "year", 
    spde = rea_spde_coast, 
    anisotropy = F,
    family = Beta(link = "logit"),
    # family = tweedie(link = "log"),
    # family = poisson(link = "log"),
    # family = binomial(link = "logit"), weights = n,
    # family = nbinom2(link = "log"),
    control = sdmTMBcontrol(step.min = 0.01, step.max = 1)
    
  )
  
  density_model <- run_extra_optimization(density_model, nlminb_loops = 0)#, newton_steps = 1)
  
  # look at gradients
  max(density_model$gradients)
  
  df$residuals <- residuals(density_model)
  m_p <- predict(density_model);
  m_p = m_p[,c("response_beta", "est")]
  m_p$back_est=inv.logit(m_p$est)
  m_p$back_abs_res=inv.logit(abs(df$residuals))
  
  #  extract some parameter estimates
  # sd <- as.data.frame(summary(TMB::sdreport(density_model$tmb_obj)))
  # r <- density_model$tmb_obj$report()
  
  if(PLOT_ALONG){
    qqnorm(df$residuals, ylim = c(-5, 5), xlim = c(-5, 5), bty = "n", pch = 20);abline(a = 0, b = 1)
    
    p1 = ggplot(df, aes_string("X", "Y", color = "residuals")) +
      geom_point(alpha = 0.8, size = round(abs(df$residuals), digits = 0)) + 
      # facet_wrap(.~ISLAND, scales = "free") +
      xlab("Eastings") +
      ylab("Northings") + 
      scale_color_gradient2() 
    
    p2 = m_p  %>% 
      ggplot(aes(response_beta*100, back_est*100)) + 
      geom_point(alpha = 0.2,aes(size=back_abs_res),color="black") + 
      coord_fixed(ratio = 1) +
      ylab("prediction") + 
      xlab("observation") + 
      geom_abline(intercept = 0, slope = 1) +
      geom_smooth(method = "lm", se = T)
    
    p1 / p2
  }
  
  
  #######################################################################
  # prediction onto new data grid
  #######################################################################
  load("data/crm/Topography_NOAA_CRM_vol10.RData") # depth covariate
  # load("data/crm/Topography_NOAA_CRM_vol10_SST_CRW_Monthly.RData") # depth and SST covariate
  # topo$x = ifelse(topo$x > 180, topo$x - 360, topo$x)
  
  grid = topo
  grid$longitude = grid$x
  grid$latitude = grid$y
  zone <- (floor((grid$longitude[1] + 180)/6) %% 60) + 1
  xy_utm = as.data.frame(cbind(utm = project(as.matrix(grid[, c("longitude", "latitude")]),
                                             paste0("+proj=utm +units=km +zone=", zone))))
  colnames(xy_utm) = c("X", "Y"); plot(xy_utm, pch = ".")
  grid = cbind(grid, xy_utm)
  
  grid_year = NULL
  years = sort(as.vector(unique(df$year)))
  # only depth covariate
  for (y in 1:length(years)) {# y = 1
    grid_y = grid[,c("X", "Y", "Topography")]
    colnames(grid_y)[3] = "mn_depth_m"
    grid_y = grid_y %>% 
      group_by(X,Y) %>% 
      summarise(mn_depth_m = mean(mn_depth_m)*-1)
    grid_y$year = years[[y]]
    grid_year = rbind(grid_year, grid_y)
  }
  dsc=scale(log(grid_year$mn_depth_m+0.0001))
  grid_year$depth_scaled = dsc
  grid_year$depth_scaled2 = grid_year$depth_scaled ^ 2
  
  ggplot(grid_year,aes(X,Y,color=mn_depth_m))+geom_point()+facet_wrap("year")
  
  #explore effects
  select_sites=sample(nrow(df),10)
  effectdata=expand.grid(year=years,mn_depth_m=0:30,ss=select_sites)#,Y=df$Y[select_sites])
  effectdata$X=df$X[effectdata$ss]
  effectdata$Y=df$Y[effectdata$ss]
  head(effectdata)
  effectdata_1000=effectdata
  effectdata_1000$X=effectdata$X+1000*runif(nrow(effectdata),min=-1,max=1)
  effectdata_1000$Y=effectdata$Y+1000*runif(nrow(effectdata),min=-1,max=1)
  effectdata_100=effectdata
  effectdata_100$X=effectdata$X+100*runif(nrow(effectdata),min=-1,max=1)
  effectdata_100$Y=effectdata$Y+100*runif(nrow(effectdata),min=-1,max=1)
  effectdata_10=effectdata
  effectdata_10$X=effectdata$X+10*runif(nrow(effectdata),min=-1,max=1)
  effectdata_10$Y=effectdata$Y+10*runif(nrow(effectdata),min=-1,max=1)
  
  eff <- predict(density_model, newdata=effectdata)
  eff$back_est=inv.logit(eff$est)
  eff1000 <- predict(density_model, newdata=effectdata_1000)
  eff1000$back_est=inv.logit(eff1000$est)
  eff100 <- predict(density_model, newdata=effectdata_100)
  eff100$back_est=inv.logit(eff100$est)
  eff10 <- predict(density_model, newdata=effectdata_10)
  eff10$back_est=inv.logit(eff10$est)
  
  peff=ggplot()+
    geom_point(data=eff,aes(mn_depth_m,back_est),color="white")
  peff10=ggplot()+
    geom_point(data=eff10,aes(mn_depth_m,back_est),color="gray75")
  peff100=ggplot()+
    geom_point(data=eff100,aes(mn_depth_m,back_est),color="gray50")
  peff1000=ggplot()+
    geom_point(data=eff1000,aes(mn_depth_m,back_est),color="black")
  (peff +peff10)/(peff100 + peff1000)
  
  pa=ggplot(eff,aes(back_est*100))+geom_histogram(binwidth=5)+facet_wrap('year')           
  pb=ggplot(eff,aes(x=mn_depth_m,back_est*100))+geom_point()+facet_wrap('year')+stat_smooth(method="loess")           
  pa/pb
  
  
  
  
  pc=ggplot(eff,aes(X,Y,color=100*back_est))+
    geom_point(size=1)+
    scale_color_gradient(low="lightblue",high="red")+
    coord_equal()+theme_bw()+facet_wrap('year')
  
  p <- predict(density_model, 
               newdata = grid_year,
               return_tmb_object = T,
               area = 0.0081)
  p$data$back_est=inv.logit(p$data$est)
  hist(p$data$back_est)
  
  ggplot(p$data,aes(X,Y,color=back_est))+geom_point()+facet_wrap("year")
  
  #Prep for Output
  p$data$sp = sp
  sdm_output = p$data
  n_knots=(rea_spde_coast$mesh$n)
  save(list=c("rea_spde_coast","density_model","sdm_output"),
       file = paste0("outputs/density_results_YEAR_MNDEPTH_", sp, "_", response_variable, "_", n_knots, "_", region, ".RData"))
  
  ###############Extra grid prediction code...
  # grid <- topo %>% subset(x < -157.5 & x > -158.5 & y > 21 & y < 22) #oahu
  # grid <- topo %>% subset(x < -154.8 & x > -156.2 & y > 18.8 & y < 20.4) #hawaii
  # grid <- topo %>% subset(x < -160.0382 & x > -160.262333 & y > 21.77143 & y < 22.03773) #Niihau
  # res = 2
  # grid$longitude = round(grid$x, digits = res)
  # grid$latitude = round(grid$y, digits = res)
  # grid = grid %>% 
  #   group_by(longitude, latitude) %>% 
  #   # subset(longitude > range(df$LONGITUDE)[1]) %>%
  #   # subset(longitude < range(df$LONGITUDE)[2]) %>%
  #   # subset(latitude > range(df$LATITUDE)[1]) %>%
  #   # subset(latitude < range(df$LATITUDE)[2]) %>%
  #   summarise(depth = mean(Topography, na.rm = T)*-1) 
  # # aggregating SST annually
  # for (y in 1:length(years)) {
  #   
  #   # y = 1
  #   
  #   grid_year_sst = grid %>% dplyr::select(contains(as.character(years[[y]])))
  #   
  #   grid_y = NULL
  #   
  #   grid_depth = grid[,3]
  #   grid_xy = grid[,which(names(grid)=="X"):which(names(grid)=="Y")]
  #   
  #   for (m in 1:12) {
  #     
  #     grid_m = cbind(grid_xy, grid_depth, grid_year_sst[,m])
  #     colnames(grid_m)[4] = "temp"
  #     grid_y = rbind(grid_y, grid_m)
  #     
  #   }
  #   
  #   grid_y = grid_y %>% 
  #     group_by(X, Y) %>% 
  #     summarise(temp = mean(temp), 
  #               depth = mean(grid_depth)*-1)
  #   
  #   grid_y$year = years[[y]]
  #   
  #   grid_year = rbind(grid_year, grid_y)
  #   
  #   rm(grid_depth, grid_xy, grid_y, grid_year_sst, grid_m)
  #   
  # }
  # 
  # # aggregating SST at each month-year time step, e.g., 03/1985 - 03/2021
  # for (t in 1:435) {
  #   
  #   # t = 1
  #   
  #   grid_t = grid[,c(442, 443, t+3)]
  #   grid_t$time = colnames(grid_t)[3]
  #   colnames(grid_t)[3] = "temp"
  #   
  #   
  #   grid_year = rbind(grid_year, grid_t)
  #   print(t/435)
  #   
  # }
  # grid_year$year = substr(grid_year$time, 1, 4)
  # grid_year = grid_year %>% subset(year %in% years)
  
  
  # # if you want to interpolate missing years...
  # grid_year_missing = NULL
  # for (y in 1:length(missing_year)) {
  #   
  #   # y = 1
  #   
  #   # Add missing years to our grid:
  #   grid_from_missing_yr <- grid_year[grid_year$year ==  missing_year[[y]]-1, ]
  #   grid_from_missing_yr$year <-  missing_year[[y]] # `L` because `year` is an integer in the data
  #   grid_year_missing <- rbind(grid_year_missing, grid_from_missing_yr)
  #   
  # }
  # grid_year = rbind(grid_year, grid_year_missing)
  
  # set the area argument to 0.0081 km2 since our grid cells are 90 m x 90 m = 0.0081 square kilometers
}
#########################################
#########################################
#########################################
# Model Selection #
#########################################
lf=list.files("outputs/",full.names = T)
lf.C=lf[grep("CORAL",lf)]

mesh_l = vector(mode = "list", length = length(lf.C))
mod_l = vector(mode = "list", length = length(lf.C))
pred_l = vector(mode = "list", length = length(lf.C))
for(i in 1:length(lf.C)){
  load(lf.C[i])
  
  mesh_l[[i]]=rea_spde_coast
  mod_l[[i]]=density_model
  pred_l[[i]]=sdm_output
  print(i)
}
#density_model=mod_l[[5]]
mod_df=data.frame(model_num=1:8,
                  knots=rep(c(2046,450,4829,7421),2),
                  cutoff=rep(c(1,5,.25,.5),2),
                  fixed=sort(rep(c("NULL","YEAR_DEPTH"),4)),
                  AIC=unlist(lapply(mod_l,AIC)),
                  filename=lf.C,range=NA,phi=NA,sigma_O=NA,sigma_E=NA,sigma_Z=NA
)
for(i in 1:nrow(mod_df)){
  tsum=tidy(mod_l[[i]],effects="ran_pars")
  mod_df$range[i]=tsum[which(tsum$term=="range")[1],"estimate"]
  mod_df$phi[i]=tsum[which(tsum$term=="phi")[1],"estimate"]
  mod_df$sigma_O[i]=tsum[which(tsum$term=="sigma_O")[1],"estimate"]
  mod_df$sigma_E[i]=tsum[which(tsum$term=="sigma_E")[1],"estimate"]
  mod_df$sigma_Z[i]=tsum[which(tsum$term=="sigma_Z")[1],"estimate"]
  print(i)
}

ggplot(mod_df,aes(cutoff,AIC,shape=fixed,fill=knots))+
  geom_point(aes(color=knots),size=5)+
  geom_label_repel(aes(label=model_num),color="white")

ggplot(mod_df,aes(AIC,sigma_E,shape=fixed,fill=knots))+
  geom_point(aes(color=knots),size=5)+
  geom_label_repel(aes(label=model_num),color="white")

ggplot(mod_df,aes(AIC,sigma_O,shape=fixed,fill=knots))+
  geom_point(aes(color=knots),size=5)+
  geom_label_repel(aes(label=model_num),color="white")

ggplot(mod_df,aes(AIC,sigma_Z,shape=fixed,fill=knots))+
  geom_point(aes(color=knots),size=5)+
  geom_label_repel(aes(label=model_num),color="white")

ggplot(mod_df,aes(AIC,range,shape=fixed,fill=knots))+
  geom_point(aes(color=knots),size=5)+
  geom_label_repel(aes(label=model_num),color="white")



plot_map_raster <- function(dat, column = "est") {
  
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_tile(aes(height = 0.8, width = 0.8), alpha = 0.8) +
    # geom_raster() +
    facet_wrap(~year) +
    coord_fixed() +
    xlab("Eastings (km)") +
    ylab("Northings (km)") + 
    scale_fill_viridis_c()+#colours = matlab.like(100), "") +
    ggdark::dark_theme_minimal()
  
}

#p$data=pred_l[[5]]
# pick out a single year to plot since they should all be the same for the slopes. Note that these are in log space.


plot_map_raster(filter(p$data, year == 2016), "zeta_s")
plot_map_raster(p$data, "exp(est)") + ggtitle("Predicted Coral Cover (%) (fixed effects + all random effects)")
plot_map_raster(p$data, "exp(est_non_rf)") + ggtitle("Prediction (fixed effects only)")
plot_map_raster(p$data, "omega_s") + ggtitle("Spatial random effects only")
plot_map_raster(p$data, "epsilon_st") + ggtitle("Spatiotemporal random effects only")

# look at just the spatiotemporal random effects:
plot_map_raster(p$data, "est_rf") + scale_fill_gradient2()

trend = plot_map_raster(filter(p$data, year == 2016), "zeta_s") + ggtitle("Linear trend")

density_map = p$data %>% 
  # group_by(X, Y) %>% 
  # summarise(est = mean(est)) %>% 
  ggplot(aes_string("X", "Y", fill = "exp(est)")) +
  geom_tile(aes(height = 1, width = 1)) +
  # geom_point(alpha = 0.5) +
  facet_wrap(~year) +
  coord_fixed() +
  xlab("Eastings (km)") +
  ylab("Northings (km)") + 
  # scale_fill_gradientn(colours = matlab.like(100), "g / m^2") +
  scale_fill_gradientn(colours = matlab.like(100), "kg/sq.m") +
  # scale_color_gradientn(colours = matlab.like(100), "# per 353 m^2") +
  # ggtitle(paste0(sp, " predicted density 2015-2019")) + 
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
  # ylab('biomass') +
  ylab('metric tonnes') +
  # ylab('total count (n)') + 
  ggtitle("Biomass estimate") +
  # ggtitle("Abundance estimate") + 
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

# png("/Users/kisei/Desktop/sdmTMB.png", height = 8, width = 12, units = "in", res = 100)
(density_map + trend )/ (relative_biomass+density_cog)
# dev.off()

