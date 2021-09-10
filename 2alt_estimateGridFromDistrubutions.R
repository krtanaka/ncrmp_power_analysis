rm(list=ls())
library(fitdistrplus)
library(fields)
library(parallel)
library(raster)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(rgdal)
library(sp)
library(colorRamps)
library(patchwork)
library(sf)
library(lubridate)
library(boot)
library(foreach)
library(tictoc)

n.cores <- detectCores()
clstr <- makePSOCKcluster(n.cores)
print(clstr)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = clstr)
foreach::getDoParRegistered()

SmithsonVerkuilen2006=function(y){
  n=length(y)
  out=(y * (n-1) + 0.5) / n
  return(out)
}

tomdist=function(pt,ref){
  D=sqrt((ref[,1]-pt[1])^2+(ref[,2]-pt[2])^2)
  return(D)
}
EstBeta_pt=function(ref_pt,df_B,NN=10){
  require(fitdistrplus)
  require(fields)
  #setup output
  outdf=data.frame(x=NA,y=NA,shape1_est=NA,shape2_est=NA,shape1_sd=NA,shape2_sd=NA,AIC=NA,BIC=NA,N=NA,LogL=NA,MaxD=NA,MeanD=NA)
  #rename and set up inputs as matrix
  names(df_B)=c("x","y","r")
  ref_ptm=matrix(ref_pt,ncol=2,byrow = F)
  df_pts=matrix(unlist(df_B[,c("x","y")]),ncol=2,byrow = F)
  #calc D and sort for top NN
  D=as.vector(tomdist(ref_ptm,df_pts))
  names(D)=1:nrow(df_B)
  Ds=sort(D)[1:NN]
  Ds_i=as.numeric(names(Ds))
  #subset to NN
  thissub=df_B[Ds_i,]
  #Fit
  thisfit=fitdist(data=thissub$r,distr="beta")
  #output
  outdf$x=ref_pt[1]
  outdf$y=ref_pt[2]
  outdf$shape1_est=thisfit$estimate[1]
  outdf$shape2_est=thisfit$estimate[2]
  outdf$shape1_sd=thisfit$sd[1]
  outdf$shape2_sd=thisfit$sd[2]
  outdf$LogL=thisfit$loglik
  outdf$AIC=thisfit$aic
  outdf$BIC=thisfit$bic
  outdf$N=length(Ds)
  outdf$MaxD=max(Ds,na.rm=T)
  outdf$MeanD=mean(Ds,na.rm=T)
  
  return(outdf)
}

# Prep Point Data Set
island= "Oahu"
load("data/rea/BenthicCover_2010-2020_Tier1_SITE_MHI_w_CRM.RData") #live coral cover, only for MHI, with CRM_Bathy data
df=subset(df,ISLAND==island)
df$CORAL_b=SmithsonVerkuilen2006(df$CORAL/100)
df$LONGITUDE360=df$LONGITUDE
df$LONGITUDE=df$LONGITUDE360-180
zone <- (floor((df$LONGITUDE[1] + 180)/6) %% 60) + 1
xy_utm = as.data.frame(cbind(utm = project(as.matrix(df[, c("LONGITUDE", "LATITUDE")]),
                                           paste0("+proj=utm +units=km +zone=", zone))))
colnames(xy_utm) = c("Xu", "Yu")
df = cbind(df, xy_utm)

#Load Survey Grid Over Which to Estimate Functions
load("data/survey_grid_w_sector_reef/survey_grid_Oahu.RData")

#for each grid point, find the NN nearest points in df, fit a beta dist, and fill in appropriate raster layers in the stack
NN=10
sgk=xyFromCell(survey_grid_kt$cell,cell=1:ncell(survey_grid_kt))#[,c("x","y")]

BlockBreak=round(seq(1,nrow(sgk),length.out = 51),0)
BlockBind=NULL
for(i_b in 1:(length(BlockBreak)-1)){
  sgk_sub=sgk[BlockBreak[i_b]:BlockBreak[i_b+1],]#sample(1:nrow(sgk),Nsubs,replace = F),]
  Nsubs=nrow(sgk_sub)
  # start=EstBeta_pt(ref_pt = sgk_sub[1,],
  #               df_B = df[,c("Xu","Yu","CORAL_b")],
  #               NN = 10)
  #progress bar
  # pb <- txtProgressBar(max = Nsubs, style = 3)
  # progress <- function(n) setTxtProgressBar(pb, n)
  # opts <- list(progress = progress)
  
  tic(paste(Nsubs,"fits"))
  parallelBetaFit <- foreach(
    i_par = 1:nrow(sgk_sub), 
    .combine = 'rbind'#,    .options.snow = opts
  ) %dopar% {
    EstBeta_pt(ref_pt = sgk_sub[i_par,],
               df_B = df[,c("Xu","Yu","CORAL_b")],
               NN = 10)
  }
  BetaGrid=parallelBetaFit#cbind(data.frame(CellNum=1:nrow(sgk_sub),sgk_sub),parallelBetaFit)
  BetaGrid$cfMean=1/(1+(BetaGrid$shape2_est)/(BetaGrid$shape1_est))
  print(paste0("Complete rows ",BlockBreak[i_b]," to ",BlockBreak[i_b+1],", of ",nrow(sgk),"."))
  BlockBind=rbind(BlockBind,BetaGrid)
  toc()
}
  
save(list = "BlockBind",file="data/rea/BenthicCover_BetaGrid.Rdata")
#Read in Island Boundaries
#"T:/Common/Maps/Island/islands"
ISL_bounds=readOGR(dsn="T:/Common/Maps/Island",
                   layer = "islands")
crs(ISL_bounds)="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

ISL_this=ISL_bounds[which(ISL_bounds$ISLAND%in%toupper(island)),]
ISL_this_utm=spTransform(ISL_this,CRS(paste0("+proj=utm +units=km +zone=",zone)))
ISL_this_sf=st_transform(st_as_sf(ISL_this),crs=paste0("+proj=utm +units=km +zone=",zone))

BetaMean=ggplot()+
  geom_point(data=BlockBind,aes(x,y,color=cfMean))+
  geom_point(data=df,aes(Xu,Yu,fill=CORAL_b),shape=21,size=3,color="black")+
  geom_sf(data=ISL_this_sf, color = "darkgreen", fill = "gray75")+
  coord_sf()+
  scale_color_gradient2(name="Coral Cover",low="#c6ffdd",mid="#fbd786",high = "#f7797d",limits=c(0,.5))+
  scale_fill_gradient2(name="Coral Cover",low="#c6ffdd",mid="#fbd786",high = "#f7797d",limits=c(0,.5))

ggsave(filename = "outputs/CoralCoverBetaMean_Oahu.jpg",plot = BetaMean)



repNA=function(x){return(rep(NA,x))}
sgk=survey_grid_kt
k=init(sgk,fun=repNA)
BetaSG=brick(k,k,k,k,k,k,k,k,k,k,values=T)
names(BetaSG)=c("shape1_est","shape2_est","shape1_sd","shape2_sd","AIC","BIC","LogL","N","MeanD","MaxD")





#Select Set of N Closest Point to Random Target
BetaFits=NULL
Nruns=100
Nmax=100
dfC=df[,c("LATITUDE","LONGITUDE","CORAL_b")]
for(i_R in 1:Nruns){
  ThisTarget=dfC[sample(1:nrow(dfC),1),]
  NNd=gcdist(lat1 = ThisTarget$LATITUDE,lon1 = ThisTarget$LONGITUDE,
             lats = dfC$LATITUDE,lons = dfC$LONGITUDE)
  names(NNd)=1:length(NNd)
  NN_Nmax=sort(NNd)[1:Nmax]
  NN_ind=as.numeric(names(NN_Nmax))
  
  BetaFit_Nmax=data.frame(Round=i_R,N=3:Nmax,shape1_est=NA,shape1_sd=NA,shape2_est=NA,shape2_sd=NA,LogL=NA,AIC=NA,BIC=NA)
  for(i_N in 3:Nmax){
    thissub=dfC[NN_ind[1:i_N],]
    thisfit=fitdist(data=thissub$CORAL_b,distr="beta")
    BetaFit_Nmax$shape1_est[i_N-2]=thisfit$estimate[1]
    BetaFit_Nmax$shape2_est[i_N-2]=thisfit$estimate[2]
    BetaFit_Nmax$shape1_sd[i_N-2]=thisfit$sd[1]
    BetaFit_Nmax$shape2_sd[i_N-2]=thisfit$sd[2]
    BetaFit_Nmax$LogL[i_N-2]=thisfit$loglik
    BetaFit_Nmax$AIC[i_N-2]=thisfit$aic
    BetaFit_Nmax$BIC[i_N-2]=thisfit$bic
  }
  BetaFits=rbind(BetaFits,BetaFit_Nmax)
  print(i_R)
}

BFcon=BetaFits %>% 
  group_by(Round) %>% 
  mutate(FinalShape1=shape1_est[N==Nmax],
         FinalShape2=shape2_est[N==Nmax],
         Shape1_RSE=sqrt((FinalShape1-shape1_est)^2),
         Shape2_RSE=sqrt((FinalShape2-shape2_est)^2),
         Shape1_pRSE=Shape1_RSE/max(Shape1_RSE),
         Shape2_pRSE=Shape2_RSE/max(Shape2_RSE))

BFconM=BFcon %>% 
  group_by(N) %>% 
  summarize(Shape1_pRSEmn=mean(Shape1_pRSE),
            Shape2_pRSEmn=mean(Shape2_pRSE))
BFconM$Shape1_pRSEmn_int=c(NA,abs(diff(BFconM$Shape1_pRSEmn)))
BFconM$Shape2_pRSEmn_int=c(NA,abs(diff(BFconM$Shape2_pRSEmn)))

ggplot(data=BFcon,aes(x=N,y=Shape1_pRSE))+
  geom_point()+
  stat_smooth(method="loess",se = T)+
  ylim(0,1.5)


ggplot(data=BFconM,aes(x=N))+
  geom_point(aes(y=Shape1_pRSEmn_int),color="blue")+
  geom_point(aes(y=Shape2_pRSEmn_int),color="red")+
  geom_path(aes(y=Shape1_pRSEmn_int),color="blue")+
  geom_path(aes(y=Shape2_pRSEmn_int),color="red")+
  scale_x_continuous(breaks=seq(0,100,by=5))+
  geom_vline(xintercept = c(10,45),color="purple")
  
ggplot(data=BFconM,aes(x=N))+
  geom_point(aes(y=Shape1_pRSEmn),color="blue")+
  geom_point(aes(y=Shape2_pRSEmn),color="red")+
  geom_path(aes(y=Shape1_pRSEmn),color="blue")+
  geom_path(aes(y=Shape2_pRSEmn),color="red")+
  scale_x_continuous(breaks=seq(0,100,by=5))+
  geom_vline(xintercept = c(10,45),color="purple")




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