##########################################
### Simulate stratified-random survey  ###
##########################################

library(SimSurvey)
library(raster)
library(data.table)
library(ggplot2)

rm(list = ls())

load("data/survey_grid_Kauai.RData")

n_sims = 10
min_sets = 3
set_den = 2/1000

options(scipen = 999, digits = 2)
set.seed(300)

sim = sim_abundance(years = 2010:2020, 
                    ages = 1:5) %>% 
  sim_distribution(grid = survey_grid_kt) %>% 
    sim_survey(trawl_dim = c(0.01, 0.0353), 
               n_sims = n_sims, 
               min_sets = min_sets, 
               set_den = set_den)

sim = sim 
length_group = "inherit"
alk_scale = "division"

strat_data_fun = strat_data
strat_means_fun = strat_means

sim_length_group <- get("length_group", envir = environment(sim$sim_length))

if (is.character(length_group) && length_group == "inherit") {
  length_group <- sim_length_group
} else {
  if (length_group != sim_length_group) {
    warning(paste0("length_group value should be set to ", 
                   sim_length_group, " to match the length group defined inside sim_abundance using sim_length", 
                   "; a mismatch in length groupings will cause issues with strat_error", 
                   " as true vs. estimated length groupings will be mismatched."))
  }
}

data <- strat_data_fun(sim, length_group = length_group, 
                       alk_scale = alk_scale)

data$setdet <- data$setdet[, c("sim", "year", 
                               "division", "strat", "strat_area", 
                               "tow_area", "set", "n"), with = FALSE]

data$lf <- merge(data$setdet[, setdiff(names(data$setdet), 
                                       "n"), with = FALSE], data$lf, by = "set", 
                 all = TRUE)

data$af <- merge(data$setdet[, setdiff(names(data$setdet), 
                                       "n"), with = FALSE], data$af, by = "set", 
                 all = TRUE)

strat_args <- list(data = data$setdet, metric = "n", 
                   strat_groups = c("sim", "year", "division", 
                                    "strat", "strat_area", "tow_area"), 
                   survey_groups = c("sim", "year"))

sim$total_strat <- do.call(strat_means_fun, strat_args)
strat_args$data <- data$lf
strat_args$strat_groups <- c(strat_args$strat_groups, "length")
strat_args$survey_groups <- c(strat_args$survey_groups, "length")
sim$length_strat <- do.call(strat_means_fun, strat_args)
strat_args$data <- data$af
strat_args$strat_groups[strat_args$strat_groups == "length"] <- "age"
strat_args$survey_groups[strat_args$survey_groups == "length"] <- "age"
sim$age_strat <- do.call(strat_means_fun, strat_args)

total <- age <- NULL
I_hat <- sim$total_strat[, list(sim, year, total)]
names(I_hat) <- c("sim", "year", "I_hat")
I <- data.frame(year = sim$years, I = colSums(sim$I))
comp <- merge(I_hat, I, by = "year")
comp$error <- comp$I_hat - comp$I
means <- error_stats(comp$error)
sim$total_strat_error <- comp
sim$total_strat_error_stats <- means
I_hat <- sim$length_strat[, list(sim, year, length, total)]
names(I_hat) <- c("sim", "year", "length", 
                  "I_hat")
sly <- expand.grid(sim = seq(max(sim$total_strat$sim)), year = sim$years, 
                   length = sim$lengths)
I_hat <- merge(sly, I_hat, by = c("sim", "year", 
                                  "length"), all = TRUE)
I_hat$I_hat[is.na(I_hat$I_hat)] <- 0
I <- as.data.frame.table(sim$I_at_length, responseName = "I")
I$year <- as.numeric(as.character(I$year))
I$length <- as.numeric(as.character(I$length))
comp <- merge(data.table(I_hat), data.table(I), by = c("year", 
                                                       "length"))
comp$error <- comp$I_hat - comp$I
means <- error_stats(comp$error)
sim$length_strat_error <- comp
sim$length_strat_error_stats <- means
I_hat <- sim$age_strat[, list(sim, year, age, total)]
names(I_hat) <- c("sim", "year", "age", 
                  "I_hat")
say <- expand.grid(sim = seq(max(sim$total_strat$sim)), year = sim$years, 
                   age = sim$ages)
I_hat <- merge(say, I_hat, by = c("sim", "year", 
                                  "age"), all = TRUE)
I_hat$I_hat[is.na(I_hat$I_hat)] <- 0
I <- as.data.frame.table(sim$I, responseName = "I")
I$year <- as.numeric(as.character(I$year))
I$age <- as.numeric(as.character(I$age))
comp <- merge(data.table(I_hat), data.table(I), by = c("year", 
                                                       "age"))
comp$error <- comp$I_hat - comp$I
means <- error_stats(comp$error)
sim$age_strat_error <- comp
sim$age_strat_error_stats <- means
sim


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
           size = 2, 
           hjust = 1,
           vjust = 1) 

print(p)
