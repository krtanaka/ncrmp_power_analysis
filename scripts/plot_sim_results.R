library(patchwork)
library(ggplot2)
library(dplyr)
library(ggpubr)

rm(list = ls())

areas = NULL
sims = NULL

for (i in 1:3) {
  
  # i = 1
  
  if (i == 1) load('outputs/sim_results_Oahu_traditional_median_PLANKTIVORE_100_biomass.RData')
  if (i == 2) load('outputs/sim_results_Oahu_downscaled_median_PLANKTIVORE_100_biomass.RData')
  if (i == 3) load('outputs/sim_results_Oahu_downscaled_alt_median_PLANKTIVORE_100_biomass.RData')
  
  area = sim_results[[1]]
  sim = sim_results[[2]]
  
  if (i == 1) {
    area$strategy = "Traditional"
    sim$strategy = "Traditional"
  }
  
  if (i == 2) {
    area$strategy = "Zone-based"
    sim$strategy = "Zone-based"
  }
  
  if (i == 3) {
    area$strategy = "Zone-triaged"
    sim$strategy = "Zone-triaged"
  }
  
  areas = rbind(area, areas)
  sims = rbind(sim, sims)
  
}

(rmse = sims %>% 
  group_by(strategy) %>% 
  summarise(rmse = sqrt(mean(error^2))) %>% 
  mutate(rmse = formatC(rmse, digits = 2)))

(p1 = areas %>% 
  ggplot(aes(x, y)) +
  # coord_fixed() + 
  geom_raster(aes(fill = strat_sets)) + 
  theme_minimal() + 
  ylab("Northing (km)") + xlab("Easting (km)") + 
  theme(legend.position = "right") + 
  facet_grid(~ strategy) + 
  # scale_fill_gradient(low = "gray", high = "red", "# of sites")) + 
  scale_fill_viridis_c("# of sites"))
  # scale_fill_gradientn(colours = topo.colors(100)))

(p2 = sims %>% 
    ggplot() + 
    geom_line(aes(year, I_hat, group = sim, color = sim), alpha = 0.5, show.legend = T) +
    scale_color_viridis_c("simulation") + 
    ggnewscale::new_scale_color() +
    geom_line(aes(year, I, color = "True biomass"), size = 2) + 
    scale_color_viridis_d("") + 
    theme_minimal() + 
    # scale_y_log10() + 
    # scale_x_log10() + 
    ylab("biomass (g)") +
    xlab("") + 
    facet_grid(~ strategy) + 
    geom_text(
      data    = rmse,
      mapping = aes(x = Inf, y = Inf, label = paste0("RMSE = ", rmse)),
      hjust   = 1,
      vjust   = 1.5))

p1/p2
