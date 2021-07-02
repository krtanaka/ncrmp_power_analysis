##########################################
### Simulate stratified-random survey  ###
##########################################

library(SimSurvey)
library(raster)
library(data.table)
library(ggplot2)

rm(list = ls())

load("data/survey_grid_kt.RData")

n_sims = 10
min_sets = 2
set_den = 2/100

options(scipen = 999, digits = 2)
set.seed(20)

sim = sim_abundance(years = 2010:2020, ages = 1:5) %>% 
  # sim_distribution(grid = make_grid(res = c(100, 100))) %>%
  sim_distribution(grid = survey_grid_kt) 

sim = sim %>%
  sim_survey(trawl_dim = c(0.01, 0.0353), n_sims = n_sims, min_sets = min_sets, set_den = set_den) %>% 
  run_strat() %>%
  strat_error()
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
  geom_point(aes(year, I_hat, color = factor(sim), alpha = 0.5), show.legend = F) +
  geom_line(aes(year, I_hat, color = factor(sim), alpha = 0.5), show.legend = F) +
  geom_point(aes(year, I), size = 5) + 
  geom_line(aes(year, I), size = 2) + 
  # scale_fill_viridis_d() +
  # scale_color_viridis_d() +
  # ggdark::dark_theme_classic() +
  theme_minimal()+ 
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

pdf(paste0("/Users/Kisei.Tanaka/coral_power_analysis/outputs/survey_sim_output_", min_sets, "_", set_den, "_", n_sims, ".pdf"), height = 4, width = 4)
print(p)
dev.off()


pop = sim_abundance(years = 1:5, ages = 1:5) %>% sim_distribution(grid = survey_grid_kt)

surveys = expand_surveys(set_den = c(0.5, 1, 2, 5, 10, 50)/100,
                         lengths_cap = c(5),
                         ages_cap = c(2))

tests = test_surveys(pop, surveys = surveys, n_sims = 2, n_loops = 5, cores = 6)

plot_total_strat_fan(tests, surveys = 6)
plot_length_strat_fan(tests)
plot_age_strat_fan(tests)

plot_survey_rank(tests, which_strat = "total")
