##############################################
### Simulate REA stratified-random survey  ###
##############################################

library(SimSurvey)
library(raster)
library(data.table)
library(ggplot2)

rm(list = ls())

load("data/survey_grid_kt.RData")

options(scipen = 999, digits = 2)

sim <- sim_abundance(ages = 1:2, years = 1:10) %>%
  sim_distribution(grid = survey_grid_kt)

sim = sim
n_sims = 2
qq = sim_logistic() # simulating catchability at age 
trawl_dim = c(0.01, 0.005) # 50 sq.m
resample_cells = FALSE
binom_error = TRUE
min_sets = 2        # minimum number of sets per strat
set_den = 2/1000    # number of sets per [grid unit] squared)

n <- age <- id <- division <- strat <- N <- n_measured <- n_aged <- NULL

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

cells <- data.table(rasterToPoints(sim$grid)); cells %>% ggplot(aes(x, y, fill = strat )) + geom_tile(aes(width = 0.5, height = 0.5))
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
sp_I$N = round(sp_I$N/1000, digits = 0)
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

m = ggplot() +
  geom_tile(data = cells, aes(x, y, fill = depth, 
                              width = 0.5, height = 0.5), alpha = 0.1) +
  geom_point(data = sets, aes(x, y, color = factor(strat))) + 
  facet_wrap(.~year) + 
  # scale_fill_viridis_c() + 
  ggdark::dark_theme_void()

t = rbind(true, sample) %>% 
  ggplot(aes(year, n, color = ts)) + 
  geom_line() + 
  ggdark::dark_theme_minimal()

library(patchwork)
m + t
