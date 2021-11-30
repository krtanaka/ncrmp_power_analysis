#######################################################################################################################
### import results from power analysis and create table of allocation per strata                                    ###
#######################################################################################################################

# !!!!!!!!!!!!! need to update this script to run with a loop and pull multiple island names - right now just works for guam !!!!!!!!

library(dplyr)

rm(list = ls()) # clean workspace

load('outputs/sim_results_guaPISCIVORE_100_biomass.RData') # load output from power analysis

# list how many sites are in which strata
a<-sim_results[[1]]
# strip depth and lat/long
b<-a %>% dplyr::select(strat,strat_sets)
# select unique records
b<-b %>% filter(!duplicated(b))
b

# load table
load('outputs/sector_key_gua.RData') # load reference table for strata numbers
names(tab)

# match strata
c<-merge(b,tab,by="strat", all.y=T,all.x=T)
# save file
write.csv(c,file="outputs/sector_allocation_gua.csv")
