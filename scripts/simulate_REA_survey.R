##############################################
### Simulate REA stratified-random survey  ###
##############################################

library(SimSurvey)
library(raster)
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)

rm(list = ls())

islands = c("Hawaii", "Kahoolawe", "Kauai", "Lanai", "Maui", "Molokai", "Niihau", "Oahu" )[sample(1:8, 1)]

load(paste0("data/survey_grid_", islands, ".RData"))
plot(survey_grid_kt, col = topo.colors(100))

# options(scipen = 999, digits = 2)

sim <- sim_abundance(years = 2010:2020, ages = 1:5) %>%
  sim_distribution(grid = survey_grid_kt)

n_sims = 10 #number of simulations
qq = sim_logistic() # simulating catchability at age 
trawl_dim = c(0.01, 0.0353) # 0.000353 sq.km (353 sq.m) from two 15-m diameter survey cylinders
resample_cells = FALSE
binom_error = TRUE
min_sets = 2        # minimum number of sets per strat
set_den = 2/1000    # number of sets per [grid unit] squared)

n <- age <- id <- division <- strat <- N <- NULL

sim <- round_sim(sim)

# abundance-at-age matrix with age specific catchabilities
I <- sim$N * qq(replicate(length(sim$years), sim$ages))
round(I, digits = 0)

# Convert abundance-at-age matrix to abundance-at-length
# abundance at each size length, sum(I) = sum(I_at_length)
# lak = Length-age-key (i.e. probability of being in a specific length group given age)
I_at_length <- convert_N(N_at_age = I, lak = sim$sim_length(age = sim$ages, length_age_key = TRUE))
round(I_at_length, digits = 0)

# Simulate survey sets

#####################################
### whats inside sim_sets() start ###
#####################################

# sets <- sim_sets(sim,
#                  resample_cells = resample_cells,
#                  n_sims = n_sims,
#                  trawl_dim = trawl_dim,
#                  set_den = set_den,
#                  min_sets = min_sets)

strat_sets <- cell_sets <- NULL

cells <- data.table(rasterToPoints(sim$grid))
# cells %>% ggplot(aes(x, y, fill = strat )) + geom_tile(aes(width = 0.5, height = 0.5))
strat_det <- cells[, list(strat_cells = .N), by = "strat"]; strat_det
strat_det$tow_area <- prod(trawl_dim); strat_det
strat_det$cell_area <- prod(res(sim$grid)); strat_det
strat_det$strat_area <- strat_det$strat_cells * prod(res(sim$grid)); strat_det
strat_det$strat_sets <- round(strat_det$strat_area * set_den); strat_det
strat_det$strat_sets[strat_det$strat_sets < min_sets] <- min_sets; strat_det #make sure minimum number of sets per strat is not 0 or 1

cells <- merge(cells, strat_det, by = c("strat")) # add "strat" "strat_cells" "tow_area" ...

i <- rep(seq(nrow(cells)), times = length(sim$years)) # number of cells * number of years
y <- rep(sim$years, each = nrow(cells)) # number of years * number of cells

cells <- cells[i, ] # increase the number of rows by number or years
cells$year <- y

i <- rep(seq(nrow(cells)), times = n_sims) # number of cells * number of simulations
s <- rep(seq(n_sims), each = nrow(cells)) # number of simulations * number of cells

cells <- cells[i, ] # increase the number of rows by number or simulations
cells$sim <- s 

#####################################
### do stratified random sampling ###
#####################################

# .SD = "Subset of Data.table"
# .N = number of instances
# strat_sets, see unique(cells$strat_sets)

sets <- cells[, .SD[sample(.N, strat_sets, replace = resample_cells)], 
              by = c("sim", "year", "strat")]

sets[, `:=`(cell_sets, .N), by = c("sim", "year", "cell")]
sets$set <- seq(nrow(sets))
sets

###################################
### whats inside sim_sets() end ###
###################################

setkeyv(sets, c("sim", "year", "cell"))

# true abundance & age data
sp_I <- data.table(sim$sp_N[, c("cell", "age", "year", "N")])
# sp_I$N = round(sp_I$N/1000, digits = 0)
hist(sp_I$N)

i <- rep(seq(nrow(sp_I)), times = n_sims) # number of rows in true abundance data * number of simulations
s <- rep(seq(n_sims), each = nrow(sp_I))

sp_I <- sp_I[i, ]
sp_I$sim <- s

setdet <- merge(sets, sp_I, by = c("sim", "year", "cell"))

if (binom_error) {
  
  #add 0.1 bc otherwise all you get is 0
  setdet$n <- stats::rbinom(rep(1, nrow(setdet)),
                            size = round(setdet$N/setdet$cell_sets), 
                            prob = ((setdet$tow_area/setdet$cell_area) * qq(setdet$age)) + 0.1) 
  
} else {
  
  setdet$n <- round((setdet$N/setdet$cell_sets) * ((setdet$tow_area/setdet$cell_area) * qq(setdet$age)))
  
}

setkeyv(setdet, "set")
setkeyv(sets, "set")

setdet <- merge(sets, setdet[, list(N = sum(N), n = sum(n)), 
                             by = "set"], by = "set")
sim$I <- I

df1 = cells %>% dplyr::select(strat, x, y, cell, year)
df2 = sp_I %>% dplyr::select(cell, year, N)

true = left_join(df1, df2)  

true = true %>% dplyr::select(cell, year, strat, N)
sample = setdet %>% dplyr::select(cell, year, strat, n)

cv <- function(x) 100*( sd(x, na.rm = T)/mean(x, na.rm = T))

true %>%
  dplyr::select(strat, year, N) %>%
  group_by(strat, year) %>%
  summarise_each(funs(cv))

sample %>%
  dplyr::select(strat, year, n) %>%
  group_by(strat, year) %>%
  summarise_each(cv)

true = true %>% 
  mutate(scale_n = scale(N)) %>% 
  group_by(year) %>% 
  summarise(n = mean(scale_n, na.rm = T))

sample = sample %>% 
  mutate(scale_n = scale(n)) %>% 
  group_by(year) %>% 
  summarise(n = mean(scale_n, na.rm = T))

true$ts = "true"
sample$ts = "survey"

true_pop = setdet %>% 
  group_by(x, y) %>% 
  summarise(N = sum(N))

d = ggplot() +
  geom_tile(data = cells, aes(x, y, fill = depth, width = 1, height = 1),
            alpha = 0.5,
            show.legend = F) + 
  geom_point(data = true_pop, aes(x, y, size = N, color = "red")) + 
  ggdark::dark_theme_void() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank()) + 
  ggtitle("Survey domain + Target population")

m = ggplot() +
  geom_point(data = sets, aes(x, y, color = factor(strat)), size = 1) +
  facet_wrap(.~year, ncol = 3) + 
  scale_color_discrete("strata") + 
  ggdark::dark_theme_void() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank()) + 
  ggtitle("Simulated Sampling Efforts")

t = rbind(true, sample) %>% 
  ggplot(aes(year, n, color = ts)) + 
  geom_line() + 
  geom_point(size = 3) + 
  scale_color_discrete("") + 
  ggdark::dark_theme_minimal() +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1))+ 
  ggtitle("Survey + True")

(d / m ) | t



# -------------------------------------------------------------------------



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
metric = "n"
strat_groups = c("sim", "year", "division", "strat", "strat_area", "tow_area")
survey_groups = c("sim", "year")
confidence = 95

Nh <- strat_area <- tow_area <- Wh <- total <- sumYh <- nh <- gh <- meanYh <- varYh <- meanYst_lcl <- meanYst <- varYst <- df <- meanYst_ucl <- sumYst <- N <- sumYst_lcl <- sumYst_ucl <- NULL

lc <- (100 - confidence)/200
uc <- (100 - confidence)/200 + (confidence/100)
d <- copy(data)
d <- d[, c(strat_groups, metric), with = FALSE]
setnames(d, names(d), c(strat_groups, "metric"))
setkeyv(d, strat_groups)
strat_tab <- d[, list(sumYh = sum(metric),
                      meanYh = mean(metric), 
                      varYh = stats::var(metric), 
                      nh = .N),
               by = strat_groups]
strat_tab[, `:=`(Nh, strat_area/tow_area)]
strat_tab[, `:=`(Wh, Nh/sum(Nh)), by = survey_groups]
strat_tab[, `:=`(total, Nh * sumYh/nh)]
strat_tab[, `:=`(gh, Nh * (Nh - nh)/nh)]

survey_tab <- strat_tab[, list(n = sum(nh),
                               N = sum(Nh), 
                               meanYst = sum(Wh * meanYh), 
                               varYst = (1/((sum(Nh))^2)) * sum(gh * varYh), 
                               df = ((sum(gh * varYh))^2)/(sum((gh^2 * varYh^2)/(nh - 1)))), by = survey_groups]

survey_tab[, `:=`(meanYst_lcl, (meanYst - (sqrt(varYst)) * abs(stats::qt(lc, df))))]
survey_tab[, `:=`(meanYst_ucl, (meanYst + (sqrt(varYst)) * abs(stats::qt(lc, df))))]
survey_tab[, `:=`(sumYst, N * meanYst)]
survey_tab[, `:=`(sumYst_lcl, (sumYst - abs(stats::qt(lc, df)) * N * sqrt(varYst)))]
survey_tab[, `:=`(sumYst_ucl, (sumYst + abs(stats::qt(lc, df)) * N * sqrt(varYst)))]
survey_tab[sapply(survey_tab, is.nan)] <- NA
survey_tab <- survey_tab[, c(survey_groups,
                             "n", "N", 
                             "df", "varYst", "meanYst", "meanYst_lcl", 
                             "meanYst_ucl", "sumYst", "sumYst_lcl", 
                             "sumYst_ucl"), with = FALSE]
survey_tab$varYst <- sqrt(survey_tab$varYst)
setnames(survey_tab, names(survey_tab), c(survey_groups, 
                                          "sets", "sampling_units", "df", "sd", 
                                          "mean", "mean_lcl", "mean_ucl", "total", 
                                          "total_lcl", "total_ucl"))
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

# ggdark::invert_geom_defaults()

p = df %>% 
  ggplot() + 
  # geom_point(aes(year, I_hat, color = factor(sim), alpha = 0.5), show.legend = F) +
  geom_line(aes(year, I_hat, color = factor(sim), alpha = 0.8), show.legend = F) +
  geom_point(aes(year, I), size = 1, color = "red") + 
  geom_line(aes(year, I), size = 1, color = "red") + 
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  ggdark::dark_theme_classic() +
  # theme_minimal()+ 
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

print(p)

