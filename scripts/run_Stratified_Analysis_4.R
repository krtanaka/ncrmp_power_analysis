#####################################################################
### Simulate stratified-random survey                             ###
### Simple Power analysis by comparing different survey "efforts" ###
#####################################################################

library(SimSurvey)
library(raster)
library(data.table)
library(ggplot2)
library(dplyr)

rm(list = ls())

set.seed(6)

islands = c("Hawaii", "Kahoolawe", "Kauai", "Lanai", "Maui", "Molokai", "Niihau", "Oahu" )[sample(1:8, 1)]
load(paste0("data/survey_grid_", islands, ".RData")) # our own survey grid, just using depth*bottom strata for now
print(islands)

# options(scipen = 999, digits = 2)

sim = sim_abundance(years = 2010:2020, ages = 1:2,
                    R = sim_R(log_mean = log(10),
                              log_sd = 0.8),
                    Z = sim_Z(log_mean = log(0.2))) %>% 
  sim_distribution(grid = survey_grid_kt)

n_sims = 10
min_sets = 2
set_den = 2/1000
trawl_dim = c(0.01, 0.0353)
resample_cells = F

n <- id <- division <- strat <- N <- n_measured <- n_aged <- NULL

sim <- round_sim(sim)

I <- sim$N

sets <- sim_sets(sim, 
                 resample_cells = resample_cells, 
                 n_sims = n_sims, 
                 trawl_dim = trawl_dim, 
                 set_den = set_den, 
                 min_sets = min_sets)

setkeyv(sets, c("sim", "year", "cell"))

sp_I <- data.table(sim$sp_N[, c("cell", "year", "N")])


i <- rep(seq(nrow(sp_I)), times = n_sims)
s <- rep(seq(n_sims), each = nrow(sp_I))

sp_I <- sp_I[i, ]
sp_I$sim <- s
setdet <- merge(sets, sp_I, by = c("sim", "year", "cell"))

setdet$n <- stats::rbinom(rep(1, nrow(setdet)), 
                          size = round(setdet$N/setdet$cell_sets), 
                          # prob = (setdet$tow_area/setdet$cell_area) * q(setdet$age))
                          prob = (setdet$tow_area/setdet$cell_area))

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

I <- data.frame(year = sim$years, I = colSums(sim$I))

comp <- merge(I_hat, I, by = "year")
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

label = paste0("ME = ", me, "\n", "MAE = ", mae, "\n", "MSE = ", mse, "\n", "RMSE = ", rmse)

strata = sim$grid_xy %>%
  mutate(x = round(x/0.5, digits = 0),
         y = round(y/0.5, digits = 0)) %>%
  group_by(x, y) %>% 
  summarise(strat = round(mean(strat), digits = 0),
            depth = mean(depth)) %>% 
  ggplot(aes(x, y)) +
  coord_fixed() + 
  scale_fill_discrete("Strata") + 
  geom_raster(aes(fill = factor(strat))) + 
  theme_minimal() + 
  ylab("Northing (km)") + xlab("Easting (km)") + 
  theme(legend.position = "bottom") + 
  ggtitle(islands)

sim_output = df %>% 
  ggplot() + 
  geom_point(aes(year, I_hat, color = factor(sim), alpha = 0.5), show.legend = F) +
  geom_line(aes(year, I_hat, color = factor(sim), alpha = 0.8), show.legend = F) +
  geom_point(aes(year, I), size = 1, color = "red") + 
  geom_line(aes(year, I), size = 1, color = "red") + 
  # scale_fill_viridis_d() +
  # scale_color_viridis_d() +
  theme_minimal() + 
  ylab("total_abundance (n)")+
  labs(
    title = "",
    subtitle = paste0("Number of simulations = ", n_sims, "\n",
                      "Min # of sets per strat = ", min_sets, "\n",
                      "number of sets per sq.km = ", set_den)
  )+
  annotate(label = label,
           geom = "text",
           x = Inf,
           y = Inf, 
           size = 4, 
           hjust = 1,
           vjust = 1) 

library(patchwork)
strata + sim_output

