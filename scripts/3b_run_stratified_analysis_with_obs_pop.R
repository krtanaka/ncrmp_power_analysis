#######################################################################################################################
### Simulate stratified-random surveys in MHI region with reconstructed fish density (count (n) or biomass (g/sq.m) ###
### Simple Power analysis by comparing different survey "efforts"                                                   ###
#######################################################################################################################

library(SimSurvey)
library(raster)
library(data.table)
library(ggplot2)
library(dplyr)
# library(ggdark)
library(patchwork)

rm(list = ls())

load("data/modeled_survey_variability.RData")
# set.seed(50)
# options(scipen = 999, digits = 2)

# pick an island ----------------------------------------------------------
island = c("Hawaii", "Kauai", "Lanai", "Maui", "Molokai", "Niihau", "Oahu" )[sample(1:7, 1)]
load(paste0("data/survey_grid_w_sector_reef/survey_grid_", island, ".RData")) #survey domain with sector & reef
load(paste0("data/survey_grid_w_zones/fish/survey_grid_", island, ".RData")) #survey domain with tom's downscaled zones
print(island)

# bring in sim$ as a place holder -----------------------------------------
sim = sim_abundance(years = 2000:2020, 
                    ages = 1:2,
                    R = sim_R(log_mean = log(50), log_sd = 0.8),
                    Z = sim_Z(log_mean = log(0.2))) %>% 
  sim_distribution(grid = survey_grid_kt)

I <- sim$N
I

# pick species and response_scale (n or g/sq.m) ---------------------------

# # fish_count
# list = list.files(path = "outputs/", pattern = "_count"); list

# fish_or_trophic_biomass
list = list.files(path = "outputs/", pattern = "_biomass"); list

# # coral_cover
# list = list.files(path = "outputs/", pattern = "_cover"); list

# # adult or juvenile coral density
# list = list.files(path = "outputs/", pattern = "_density"); list

i = 3

load(paste0("outputs/", list[i]))
sp = strsplit(list[i], split = "_")[[1]][3]; sp
response_scale = strsplit(list[i], split = "_")[[1]][4]; response_scale


# replace sim$ with sdmTMB outputs  -----------------------------------------

sdm = sdm_output[,c("X", "Y", "longitude", "latitude", "year", "est" )]; rm(sdm_output)
colnames(sdm)[1:2] = c("x", "y")
sdm$est = exp(sdm$est)

# replace sim$years and sim$ages
sim$years = sort(unique(sdm$year))
sim$ages = 1

# replace sim$N
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

# df %>%
#   group_by(x, y) %>%
#   mutate(x = round(x, 0),
#          y = round(y, 0)) %>% 
#   summarise(est = median(est)) %>%
#   ggplot(aes(x, y, fill = est)) +
#   geom_raster() +
#   scale_fill_gradientn(colours = colorRamps::matlab.like(100)) +
#   coord_fixed() +
#   ggdark::dark_theme_minimal()  

N = df %>% group_by(year) %>% summarise(age = sum(est)) 
N = matrix(N$age, nrow = 1, ncol = 9)
rownames(N) <- "1"
colnames(N) = sort(unique(sdm$year))
names(dimnames(N)) = c("age", "year")
N
sim$N = N

# replace sim$sp_N
sp_N = df %>% group_by(year, cell) %>% summarise(N = sum(est))
sim$sp_N = sp_N

# replace sim$I
I <- sim$N
I

# simulate stratified random surveys --------------------------------------

load("data/survey_effort_MHI_2014-2019.RData")
effort = c("high", "median", "low")[2]
t_sample = survey_effort_MHI %>% subset(ISLAND == island) %>% select(effort) %>% as.character() %>% as.numeric()

n_sims = 100 # number of simulations
total_sample = t_sample # total sample efforts you want to deploy
min_sets = 2 # minimum number of sets per strat
set_den = 2/1000 # number of sets per [grid unit = km] squared)
trawl_dim = c(0.01, 0.0353) # 0.000353 sq.km (353 sq.m) from two 15-m diameter survey cylinders
resample_cells = F

n <- id <- division <- strat <- N <- NULL

# sets <- sim_sets(sim,
#                  resample_cells = resample_cells,
#                  n_sims = n_sims,
#                  trawl_dim = trawl_dim,
#                  set_den = set_den,
#                  min_sets = min_sets)

strat_sets <- cell_sets <- NULL

cells <- data.table(rasterToPoints(sim$grid))
cells$sd = predict(g, cells); sd = cells[,c("strat", "sd")]; sd = sd %>% group_by(strat) %>% summarise(sd = mean(sd)) # add modeled fish density variability
strat_det <- cells[, list(strat_cells = .N), by = "strat"]; strat_det
strat_det$tow_area <- prod(trawl_dim); strat_det
strat_det$cell_area <- prod(res(sim$grid)); strat_det
strat_det$strat_area <- strat_det$strat_cells * prod(res(sim$grid)); strat_det
strat_det = right_join(strat_det, sd); strat_det
strat_det$strat_sets <- round(strat_det$strat_area * set_den); strat_det
strat_det$weight = strat_det$strat_area * strat_det$sd
strat_det$strat_sets = round((total_sample * strat_det$weight) / sum(strat_det$weight), 0); strat_det
# strat_det$strat_sets = round((total_sample * strat_det$strat_area) / sum(strat_det$strat_area), 0); strat_det
strat_det$strat_sets[strat_det$strat_sets < min_sets] <- min_sets; strat_det # make sure minimum number of sets per strat is not 0 or 1

cells <- merge(cells, strat_det, by = c("strat")) # add "strat" "strat_cells" "tow_area" ...

i <- rep(seq(nrow(cells)), times = length(sim$years)) # number of cells * number of years
y <- rep(sim$years, each = nrow(cells)) # number of years * number of cells

cells <- cells[i, ] # increase the number of rows by number or years
cells$year <- y

i <- rep(seq(nrow(cells)), times = n_sims) # number of cells * number of simulations
s <- rep(seq(n_sims), each = nrow(cells)) # number of simulations * number of cells

cells <- cells[i, ] # increase the number of rows by number or simulations
cells$sim <- s

# .SD = "Subset of Data.table"
# .N = number of instances
# strat_sets, see unique(cells$strat_sets)

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
  geom_raster(aes(fill = factor(strat))) + 
  theme_minimal() + 
  ylab("Northing (km)") + xlab("Easting (km)") + 
  theme(legend.position = "none") + 
  labs(
    title = "",
    subtitle = paste0(paste0(island, "\n", "# of strata = ", length(unique(sim$grid_xy$strat)))))

if (response_scale == "biomass") ylab_scale = "biomass (g)"
if (response_scale == "count") ylab_scale = "abundance (n)"

sim_output = df %>% 
  ggplot() + 
  geom_line(aes(year, I_hat, color = factor(sim), alpha = 0.2), show.legend = F) +
  geom_line(aes(year, I), size = 2, color = "red") + 
  # scale_color_viridis_d() + 
  # dark_theme_minimal() + 
  theme_minimal() + 
  ylab(ylab_scale) +
  labs(
    title = "",
    subtitle = paste0("Survey target = ", sp, "\n",
                      "Total # of surveyed sites = ", total_sample, "\n",
                      "Min # of sets per strat = ", min_sets, "\n",
                      "Number of simulations = ", n_sims))+
  annotate(label = label,
           geom = "text",
           x = Inf,
           y = Inf, 
           size = 4, 
           hjust = 1,
           vjust = 1) 

# png(paste0("outputs/", sp, "_", island, ".png"), res = 100, units = "in", height = 4, width = 8)
strata + sim_output
# dev.off()
