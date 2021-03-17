source("C:/Users/Kisei/Google Drive/R/FVCOM/MENH-FVCOM - Fall.R")
library(mgcv)
library(foreign)

# setwd("~/Google Drive/Research/Trawl Survey Data/ME & NH/GAM")
# setwd("E:/Kisei/Research/Trawl Survey Data/ME & NH/GAM")
# setwd("D:/Kisei/GDrive/Research/Trawl Survey Data/ME & NH/GAM")
setwd("C:/Users/Kisei/Google Drive/Research/Trawl Survey Data/ME & NH/GAM")

# load("C:/Users/Kisei/Google Drive/Research/Manuscripts/CJFAS/other materials/TweedieGAMResult.RData") #w/o fvcom
load("C:/Users/Kisei/Google Drive/Research/Manuscripts/CJFAS/other materials/TweedieGAMResult_w_FVCOM_knots.RData") #models fitted w/ fvcom data, k= added

yearlist=c("1978", "1979", "1980", "1981", "1982", "1983", "1984", 
           "1985", "1986", "1987", "1988", "1989", "1990", "1991", 
           "1992", "1993", "1994", "1995", "1996", "1997", "1998", 
           "1999", "2000", "2001", "2002", "2003", "2004", "2005", 
           "2006", "2007", "2008", "2009", "2010", "2011", "2012", 
           "2013")

for(year in yearlist){
  
  predict_env = get(paste0("ts", year))
  predict_env$depth = predict_env$depth *-1
  predict_env = predict_env[which(predict_env$depth > 0),]
  predict_env = predict_env[which(predict_env$salinity_b > 20),]
  colnames(predict_env)[1] = "La"
  colnames(predict_env)[2] = "Lo"
  colnames(predict_env)[3] = "Do"
  colnames(predict_env)[4] = "Se"
  colnames(predict_env)[6] = "D"
  colnames(predict_env)[8] = "Ds"
  colnames(predict_env)[9] = "Te"
  colnames(predict_env)[10] = "S"

  predict_est = c() #overall estimates

  looplist = c("juv.f", "juv.m", "adu.f", "adu.m")
  
  for (x in looplist) #run one group each time for all datasets
  {
    if (x == "juv.f")  {
      gam1 = fl_f_juv
    }
    if (x == "juv.m")  {
      gam1 = fl_m_juv
    }
    if (x == "adu.f")  {
      gam1 = fl_f_adu
    }
    if (x == "adu.m")  {
      gam1 = fl_m_adu
    }
    
    pred = data.frame(predict(gam1, type = "response", newdata = predict_env, se.fit = T))

    predict_est = c(predict_est, (pred[1])) 
    
    pred = NULL
  }
  
  out = data.frame(predict_env, rowSums(data.frame(predict_est)), predict_est)
  
  colnames(out) = c(names(predict_env), 
                    "Total","FL_JUV_F","FL_JUV_M","FL_ADU_F","FL_ADU_M") 
  
  print(year)
  
  write.csv(out, file = paste("predict_FL_", year, ".csv", sep = ''), row.names = FALSE)
  
}

summary(out)
