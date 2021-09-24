#######################################################################################################################
### Simulate stratified-random surveys in MHI region with reconstructed fish density (count (n) or biomass (g/sq.m) ###
### Simple Power analysis by comparing different survey "efforts"                                                   ###
#######################################################################################################################

library(SimSurvey)
library(raster)
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)
library(ggnewscale)

rm(list = ls())

load("data/modeled_survey_variability.RData") #modeled at grid scale

# pick an island
islands = c("Hawaii", "Kauai", "Lanai", "Maui", "Molokai", "Niihau", "Oahu" )#[sample(1:7, 1)]

isl_power = NULL

for (isl in 1:length(islands)){
  
  # isl = 5
  
  island = islands[isl]
  
  print(island)
  
  design = c("traditional", "downscaled", "downscaled_alt")[1:2]
  
  power = NULL
  
  for (ds in 1:length(design)) {
    
    # ds = 3
    
    if (design[ds] == "traditional") load(paste0("data/survey_grid_w_sector_reef/survey_grid_", island, ".RData")) #survey domain with sector & reef & depth_bins
    if (design[ds] == "downscaled") load(paste0("data/survey_grid_w_zones/fish/survey_grid_", island, ".RData")) #survey domain with tom's downscaled zones
    if (design[ds] == "downscaled_alt") load(paste0("data/survey_grid_w_zones_alt/fish/survey_grid_", island, ".RData")) #survey domain with tom's downscaled zones
    
    # bring in sim$ as a place holder
    sim = sim_abundance(years = 2000:2020, ages = 1:2) %>% sim_distribution(grid = survey_grid_kt)
    
    #population variable that we will replace...
    I <- sim$N
    I
    
    # pick species and response_scale (n or g/sq.m)
    # fish_or_trophic_biomass
    list = list.files(path = "outputs/", pattern = "_biomass"); list
    
    ds_sp = NULL
    
    for (t in 3:6) {
      
      # t = 5
      
      load(paste0("outputs/", list[t]))
      sp = strsplit(list[t], split = "_")[[1]][3]; sp
      response_scale = strsplit(list[t], split = "_")[[1]][4]; response_scale
      
      # replace sim$ with sdmTMB outputs
      sdm = sdm_output[,c("X", "Y", "longitude", "latitude", "year", "est" )]; rm(sdm_output)
      colnames(sdm)[1:2] = c("x", "y")
      sdm$est = exp(sdm$est)*prod(res(survey_grid_kt))*1000^2
      
      sim$years = sort(unique(sdm$year))
      sim$ages = 1
      
      sim_grid = sim$grid_xy
      
      sim_grid = sim_grid %>%
        mutate(x = round(x, 1),
               y = round(y, 1)) %>%
        group_by(x, y) %>%
        summarise(cell = round(median(cell), 0))
      
      sdm_grid = sdm %>%
        mutate(x = round(x, 1),
               y = round(y, 1)) %>%
        group_by(x, y, year) %>%
        summarise(est = sum(est))
      
      sdm_grid = as.data.frame(sdm_grid)
      sim_grid = as.data.frame(sim_grid)
      
      head(sdm_grid)
      head(sim_grid)
      
      df = merge(sim_grid, sdm_grid)
      
      N = df %>% group_by(year) %>% summarise(age = sum(est)) 
      N = matrix(N$age, nrow = 1, ncol = 9)
      rownames(N) <- "1"
      colnames(N) = sort(unique(sdm$year))
      names(dimnames(N)) = c("age", "year")
      N
      sim$N = N
      
      sp_N = df %>% group_by(year, cell) %>% summarise(N = sum(est)) 
      sim$sp_N = sp_N
      
      I <- sim$N
      I
      
      load("data/survey_effort_MHI_2014-2019.RData")
      
      effort = c("high", "median", "low")[1]
      
      # t_sample = survey_effort_MHI %>%
      #   subset(ISLAND == island) %>%
      #   dplyr::select(effort) %>%
      #   as.character() %>%
      #   as.numeric() %>%
      #   round(0)
      
      t_samples = seq(10, 500, by = 10)
      
      df_power = data.frame(N = t_samples, RMSE = NA)
      
      for (i_sampling in 1:length(t_samples)){
        
        t_sample = t_samples[i_sampling]
        n_sims = 10 # number of simulations
        total_sample = t_sample # total sample efforts you want to deploy
        min_sets = 0 # minimum number of sets per strat
        trawl_dim = c(0.01, 0.0353) # 0.000353 sq.km (353 sq.m) from two 15-m diameter survey cylinders
        resample_cells = T
        
        n <- id <- division <- strat <- N <- NULL
        
        strat_sets <- cell_sets <- NULL
        
        cells <- data.table(rasterToPoints(sim$grid))
        cells$sd = predict(g, cells); sd = cells[,c("strat", "sd")]; sd = sd %>% group_by(strat) %>% summarise(sd = mean(sd,na.rm=T)) # add modeled trophic biomass variability
        strat_det <- cells[, list(strat_cells = .N), by = "strat"]; strat_det
        strat_det$tow_area <- prod(trawl_dim); strat_det
        strat_det$cell_area <- prod(res(sim$grid)); strat_det
        strat_det$strat_area <- strat_det$strat_cells * prod(res(sim$grid)); strat_det
        strat_det = right_join(strat_det, sd); strat_det
        strat_det$weight = strat_det$strat_area * strat_det$sd; strat_det
        strat_det$strat_sets = round((total_sample * strat_det$weight) / sum(strat_det$weight), 0); strat_det
        # strat_det$strat_sets = round((total_sample * strat_det$strat_area) / sum(strat_det$strat_area), 0); strat_det
        strat_det$strat_sets[strat_det$strat_sets < min_sets] <- min_sets; strat_det # make sure minimum number of sets per strat is not 0 or 1
        strat_table = strat_det %>% dplyr::select(strat, strat_sets); strat_table
        
        if (design == "downscaled_alt") {
          
          # randomly drop some of tom's zones then allocate sampling units by area
          strat_det$drop = rbinom(length(unique(strat_det$strat)), 1, prob = 2/3); strat_det
          # strat_det$drop = rbinom(length(unique(strat_det$strat)), 1, prob = strat_det$weight); strat_det
          strat_det$weight = strat_det$strat_area * strat_det$sd * strat_det$drop; strat_det
          strat_det$weight = (strat_det$weight - min(strat_det$weight)) / (max(strat_det$weight)-min(strat_det$weight)); strat_det
          strat_det$strat_sets = round((total_sample * strat_det$weight) / sum(strat_det$weight), 0); strat_det
          # strat_det$weight = ifelse(strat_det$drop == 0, 0, strat_det$weight); strat_det
          strat_table = strat_det %>% dplyr::select(strat, strat_sets); strat_table
          
        }
        
        cells <- merge(cells, strat_det, by = c("strat")) # add "strat" "strat_cells" "tow_area" ...
        
        i <- rep(seq(nrow(cells)), times = length(sim$years)) # number of cells * number of years
        y <- rep(sim$years, each = nrow(cells)) # number of years * number of cells
        
        cells <- cells[i, ] # increase the number of rows by number or years
        cells$year <- y
        
        i <- rep(seq(nrow(cells)), times = n_sims) # number of cells * number of simulations
        s <- rep(seq(n_sims), each = nrow(cells)) # number of simulations * number of cells
        
        cells <- cells[i, ] # increase the number of rows by number or simulations
        cells$sim <- s
        
        sets <- cells[, .SD[sample(.N, strat_sets, replace = resample_cells)], 
                      by = c("sim", "year", "strat")]
        
        sets[, `:=`(cell_sets, .N), by = c("sim", "year", "cell")]
        sets$set <- seq(nrow(sets))
        sets
        
        setkeyv(sets, c("sim", "year", "cell"))
        
        sp_I <- data.table(sim$sp_N[, c("cell", "year", "N")])
        
        i <- rep(seq(nrow(sp_I)), times = n_sims)
        s <- rep(seq(n_sims), each = nrow(sp_I))
        
        sp_I <- sp_I[i, ]
        sp_I$sim <- s
        setdet <- merge(sets, sp_I, by = c("sim", "year", "cell"))
        
        setdet$n <- stats::rbinom(rep(1, nrow(setdet)), 
                                  size = round(setdet$N/1),
                                  prob = (setdet$tow_area/setdet$cell_area)*0.5)## DEBUGGGING !!!!!
        
        setdet$detection = ifelse(setdet$N > 0 & setdet$n == 0, 0, 1)
        setdet$detection = ifelse(setdet$N > 0 & setdet$n == 0, 0, 1)
        
        setkeyv(setdet, "set")
        setkeyv(sets, "set")
        
        sim$I <- I
        
        sim$setdet <- setdet
        
        setdet <- sim$setdet
        
        data = list(setdet = setdet)
        
        data$setdet <- data$setdet[, c("sim", 
                                       "year", 
                                       "division",
                                       "strat",
                                       "strat_area", 
                                       "tow_area",
                                       "set",
                                       "n"), 
                                   with = FALSE]
        data = data$setdet
        data
        
        metric = "n"
        strat_groups = c("sim", "year", "division", "strat", "strat_area", "tow_area")
        survey_groups = c("sim", "year")
        
        confidence = 95
        
        Nh <- strat_area <- tow_area <- Wh <- total <- sumYh <- nh <- gh <- meanYh <- varYh <- meanYst_lcl <- meanYst <- varYst <- df <- meanYst_ucl <- sumYst <- N <- sumYst_lcl <- sumYst_ucl <- NULL
        
        lc <- (100 - confidence)/200; lc
        uc <- (100 - confidence)/200 + (confidence/100); uc
        
        d <- copy(data)
        d <- d[, c(strat_groups, metric), with = FALSE]
        
        setnames(d, names(d), c(strat_groups, "metric"))
        setkeyv(d, strat_groups)
        
        strat_tab <- d[, list(sumYh = sum(metric), # sum of samples (n)
                              meanYh = mean(metric, na.rm = T), # mean of samples (n)
                              varYh = stats::var(metric, na.rm = T), # variance of samples (n)
                              nh = .N), # number of strata
                       by = strat_groups]; strat_tab
        
        strat_tab[, `:=`(Nh, strat_area/tow_area)]; strat_tab
        strat_tab[, `:=`(Wh, Nh/sum(Nh)), by = survey_groups]; strat_tab
        strat_tab[, `:=`(total, Nh * sumYh/nh)]; strat_tab
        strat_tab[, `:=`(gh, Nh * (Nh - nh)/nh)]; strat_tab
        
        survey_tab <- strat_tab[, list(n = sum(nh, na.rm = T),
                                       N = sum(Nh, na.rm = T), 
                                       meanYst = sum(Wh * meanYh), 
                                       varYst = (1/((sum(Nh))^2)) * sum(gh * varYh), 
                                       df = ((sum(gh * varYh))^2)/(sum((gh^2 * varYh^2)/(nh - 1)))), by = survey_groups]; survey_tab
        
        survey_tab[, `:=`(meanYst_lcl, (meanYst - (sqrt(varYst)) * abs(stats::qt(lc, df))))]; survey_tab
        survey_tab[, `:=`(meanYst_ucl, (meanYst + (sqrt(varYst)) * abs(stats::qt(lc, df))))]; survey_tab
        survey_tab[, `:=`(sumYst, N * meanYst)]; survey_tab
        survey_tab[, `:=`(sumYst_lcl, (sumYst - abs(stats::qt(lc, df)) * N * sqrt(varYst)))]; survey_tab
        survey_tab[, `:=`(sumYst_ucl, (sumYst + abs(stats::qt(lc, df)) * N * sqrt(varYst)))]; survey_tab
        
        survey_tab[sapply(survey_tab, is.nan)] <- NA
        
        survey_tab <- survey_tab[, c(survey_groups,
                                     "n", 
                                     "N", 
                                     "df",
                                     "varYst", 
                                     "meanYst",
                                     "meanYst_lcl", 
                                     "meanYst_ucl",
                                     "sumYst",
                                     "sumYst_lcl", 
                                     "sumYst_ucl"), 
                                 with = FALSE]
        
        survey_tab$varYst <- sqrt(survey_tab$varYst)
        
        setnames(survey_tab, names(survey_tab), c(survey_groups, 
                                                  "sets", 
                                                  "sampling_units", 
                                                  "df", 
                                                  "sd", 
                                                  "mean",
                                                  "mean_lcl",
                                                  "mean_ucl",
                                                  "total", 
                                                  "total_lcl", 
                                                  "total_ucl"))
        survey_tab
        
        sim$total_strat = survey_tab
        
        total <- NULL
        
        I_hat <- sim$total_strat[, list(sim, year, total)]
        names(I_hat) <- c("sim", "year", "I_hat")
        
        #I <- data.frame(year = sim$years, I = colSums(sim$I))
        Iy <- data.frame(year = sim$years, I = colSums(sim$I))
        
        comp <- merge(I_hat, Iy, by = "year")
        comp$error <- comp$I_hat - comp$I
        means <- error_stats(comp$error)
        sim$total_strat_error <- comp
        sim$total_strat_error_stats <- means
        I_hat <- sim$length_strat[, list(sim, year, length, total)]
        
        sim$total_strat_error_stats
        sim$total_strat_error
        df = sim$total_strat_error
        
        me = formatC(sim$total_strat_error_stats[1], digits = 3)
        mae = formatC(sim$total_strat_error_stats[2], digits = 3)
        mse = formatC(sim$total_strat_error_stats[3], digits = 3)
        rmse = formatC(sim$total_strat_error_stats[4], digits = 3)
        
        df_power$RMSE[i_sampling]=rmse
        df_power$TrueBiomass[i_sampling]=comp$I
        df_power$SurveyBiomass[i_sampling]=comp$I_hat
        print(i_sampling)
      }
      
      df_power$sp = sp
      
      ds_sp = rbind(ds_sp, df_power)
      
      
    }
    
    ds_sp$design = design[ds]
    ds_sp$strata = length(unique(sim$grid_xy$strat))
    
    power = rbind(ds_sp, power)
    
    
  }
  
  power$isl = island
  
  isl_power = rbind(power,isl_power )
  
}

save(isl_power, file = paste0('outputs/rmse_power_results_', Sys.Date(), '.RData'))
load(paste0('outputs/rmse_power_results_2021-09-10.RData'))

load("data/survey_effort_MHI_2014-2019.RData")

# these numbers are from "plot_REA_survey_efforts.R"
survey_effort_MHI_year = data.frame(Year = c("2013", "2016", "2019"),
                                    N = c(498, 400, 487))

low = survey_effort_MHI %>%
  select(ISLAND, low) %>% 
  rename(isl = ISLAND, 
         sites = low) %>% 
  mutate(effort = "low")

mid = survey_effort_MHI %>%
  select(ISLAND, median) %>% 
  rename(isl = ISLAND, 
         sites = median) %>% 
  mutate(effort = "mid")

high = survey_effort_MHI %>%
  select(ISLAND, high) %>% 
  rename(isl = ISLAND, 
         sites = high) %>% 
  mutate(effort = "high")

efforts = rbind(low, mid, high)

traditional_strata = isl_power %>% 
  subset(design == "traditional") %>%
  group_by(design, isl) %>% 
  summarise(traditional = mean(strata)) %>% 
  as.data.frame() %>% 
  select(isl, traditional)
rename(strata = traditional)

downscaled_strata = isl_power %>% 
  subset(design == "downscaled") %>%
  group_by(design, isl) %>% 
  summarise(zone_based = mean(strata)) %>% 
  as.data.frame() %>% 
  select(isl, zone_based)
rename(strata = downscaled)

reduced_strata = isl_power %>% 
  subset(design == "downscaled_alt") %>%
  group_by(design, isl) %>% 
  summarise(zone_triaged = mean(strata)) %>% 
  as.data.frame() %>% 
  select(isl, zone_triaged)
rename(strata = downscaled)

strata_num = merge(downscaled_strata, traditional_strata)
strata_num = merge(strata_num, reduced_strata)

isl_power$design = ifelse(isl_power$design == "downscaled", "zone-based", isl_power$design)
isl_power$design = ifelse(isl_power$design == "downscaled_alt", "zone-triaged", isl_power$design)

isl_power %>%
  # subset(isl == "Oahu") %>% 
  # subset(sp == "PISCIVORE") %>% 
  mutate(RMSE = as.numeric(RMSE)) %>% 
  ggplot() + 
  # geom_smooth(aes(N, RMSE, color = design), show.legend = T, se = T) +
  # geom_point(aes(N, RMSE, color = design), alpha = 0.2) +
  # geom_hex(aes(N, RMSE, color = design, fill = design), alpha = 0.2, bins = 50) +
  facet_wrap(sp ~ isl, scales = "free_y", ncol = 7) +
  # ggnewscale::new_scale_color() +
  geom_vline(aes(xintercept = sites, color = effort), data = efforts) +
  # geom_vline(aes(xintercept = N, color = Year), data = survey_effort_MHI_year) +
  # geom_text_repel(data = survey_effort_MHI_year, 
  #                 aes(x = N, y = 0, label = Year), 
  #                 fontface = "bold",   
  #                 nudge_x = c(0, 0, 0),
  #                 nudge_y = c(6, 3, 3)) +
annotate("segment", 
         x = 498, xend = 498, y = 250000000, yend = 100000000,
         colour = "gray", 
         arrow = arrow()) + 
  annotate("segment",
           x = 400, xend = 400, y = 150000000, yend = 100000000,
           colour = "gray",
           arrow = arrow()) + 
  annotate("segment", 
           x = 487, xend = 487, y = 200000000, yend = 100000000,
           colour = "gray", 
           arrow = arrow()) + 
  annotate("text", 
           x = c(400, 487, 498),
           y = c(160000000, 210000000, 260000000), 
           label = c("2016", "2019", "2013")) + 
  scale_color_discrete("") + 
  xlab("Sampling Efforts") + 
  # scale_y_log10() + 
  # scale_x_log10() + 
  # scale_x_log10(breaks = trans_breaks('log10', function(x) 10^x),
  #               labels = trans_format('log10', math_format(10^.x)), "Sampling Efforts") +
  # scale_y_log10(breaks = trans_breaks('log10', function(x) 10^x),
  #               labels = trans_format('log10', math_format(10^.x)), "RMSE") +
  # geom_text(data = strata_num,
  #           aes(label = paste0("\n d=", downscaled, "\n t=", traditional, "\n r=,", reduced)),
  #           x = -Inf, y = -Inf,
  #           hjust = -0.1,
  #           vjust = -0.2,
  #           size = 3) +
  theme_pubr() +
  guides(color = guide_legend(override.aes = list(fill = NA))) + 
  theme(legend.position = c(0.8, 0.9))
