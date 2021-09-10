rm(list=ls())
library(fitdistrplus)
SmithsonVerkuilen2006=function(y){
  n=length(y)
  out=(y * (n-1) + 0.5) / n
  return(out)
}

#Load sector shapefiles
# pick survey design ------------------------------------------------------

island= "oah"
target = c("fish", "ben_den", "ben_cov")[3]

#Load Appropriate Datafiles
load("data/rea/SURVEY MASTER.RData")
# if(target=="fish"){
#   load("data/rea/ALL_REA_FISH_RAW_SST.RData")
#   load(paste0("data/gis_sector/raster/",toupper(island), ".RData"))
# }
# if(target=="ben_den"){
#   load("data/rea/BenthicREA_sitedata_TAXONCODE.RData_MHI_w_CRM.RData") #live coral cover, only for MHI, with CRM_Bathy data
# }
if(target=="ben_cov"){
  load("data/survey_grid_w_sector_reef/survey_grid_Oahu.RData")
  load("data/rea/BenthicCover_2010-2020_Tier1_SITE_MHI_w_CRM.RData") #live coral cover, only for MHI, with CRM_Bathy data
  df$CORAL_b=SmithsonVerkuilen2006(df$CORAL/100)
  df$STRATA_NAME=paste0(df$SEC_NAME,"_",df$DEPTH_BIN)
  df$STRATA_YEAR=paste0(df$SEC_NAME,"_",df$DEPTH_BIN,"_",df$OBS_YEAR)
  
  uSTR=unique(df$STRATA_YEAR)
  n_=function(x){return(length(grep(pattern="_",x)))}
  apply(uSTR,grep,pattern="_")
  SDYmat=matrix(unlist(strsplit(uSTR,"_")),ncol=4,byrow = F)
  #Fit Beta Distribution for Each Strata for Cover
  BetaFits=data.frame(STRATA=uSTR,
                      SECTOR=N=NA,shape1_est=NA,shape1_sd=NA,shape2_est=NA,shape2_sd=NA,LogL=NA,AIC=NA,BIC=NA)
  for(b_i in 1:length(uSTR)){
    thissub=subset(df,STRATA_NAME==uSTR[b_i])
    thisfit=fitdist(data=thissub$CORAL_b,distr="beta")
    BetaFits$N=length(thissub$CORAL_b)
  }
  
}