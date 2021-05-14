##########################################
### Simulate stratified-random survey  ###
##########################################

library(SimSurvey)
library(raster)
library(data.table)

rm(list = ls())

load("data/survey_grid_kt.RData")

source('scripts/key_functions.R')

options(scipen = 999)
options(digits = 2)

set.seed(438)

plot(survey_grid)
plot(survey_grid_kt)

R_fun <- sim_R(log_mean = log(c(10000000, 2000000, 100000, 10000000, 50000)), 
               log_sd = 0, random_walk = T, plot = TRUE)
R_fun(years = 1:5)
N0_fun <- sim_N0(N0 = c(50000, 10000, 1000, 500, 100), plot = TRUE)
N0_fun(R0 = 100000, ages = 1:6)

pop <- sim_abundance(ages = 1:6, years = 1:5,
                     R = R_fun, N0 = N0_fun)
pop$N

pop <- sim_abundance(ages = 1:6, 
                     years = 2016:2020,
                     R = R_fun,
                     N0 = N0_fun,
                     Z = sim_Z(log_mean = log(0.8),
                               phi_age = 0.5,
                               phi_year = 0.3),
                     growth = sim_vonB(L0 = 1, 
                                       K = 0.1,
                                       Linf = 10))

pop[1:8] <- lapply(pop[1:8], round, 1)
pop[1:8]

sim1a = pop %>%
  sim_distribution_rea(grid = survey_grid_kt,
                       ays_covar = sim_ays_covar(range = 300, 
                                                 phi_year = 0.2,
                                                 phi_age = 0.2),
                       depth_par = sim_parabola(mu = 10, sigma = 10))

sim1b = pop %>%
  sim_distribution_rea(grid = survey_grid,
                       ays_covar = sim_ays_covar(range = 300, 
                                                 phi_year = 0.2,
                                                 phi_age = 0.2),
                       depth_par = sim_parabola(mu = 10, sigma = 10))

sim1a[10:12] <- lapply(sim1a[10:12], round, 2)
sim1b[10:12] <- lapply(sim1b[10:12], round, 2)

sim1a[10:12]
sim1b[10:12]

sim2a = sim1a %>% 
  sim_survey_rea(min_sets = 10, 
                 ages_cap = 500,
                 lengths_cap = 100,
                 set_den = 2/1000)

sim2b = sim1b %>% 
  sim_survey_rea(min_sets = 10, 
                 ages_cap = 500,
                 lengths_cap = 100,
                 set_den = 2/1000)
sim2a$I
sim2b$I

sim2a$I_at_length
sim2b$I_at_length

sim2a$setdet
sim2b$setdet


sim3a = sim2a %>% 
  run_strat()

sim3b = sim2b %>% 
  run_strat()

sim3a$samp
sim3b$samp

sim3a$samp
sim3b$samp

sim3a$samp
sim3b$samp

sim3a$samp
sim3b$samp


sim %>% strat_error()

pop = sim_abundance() %>% sim_distribution()
surveys = expand_surveys(set_den = c(0.5, 1)/1000,
                         lengths_cap = c(5, 10),
                         ages_cap = c(2, 5))
tests = test_surveys(pop, surveys = surveys, n_sims = 1, n_loops = 10)
plot_total_strat_fan(tests)
plot_length_strat_fan(tests, years = 1:5, lengths = 1:100)
plot_survey_rank(tests)
plot_error_surface(tests, plot_by = "rule")
