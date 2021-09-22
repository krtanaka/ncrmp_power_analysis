
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

areas = areas %>% 
  ggplot(aes(x, y)) +
  # coord_fixed() + 
  geom_raster(aes(fill = strat_sets)) + 
  theme_minimal() + 
  ylab("Northing (km)") + xlab("Easting (km)") + 
  theme(legend.position = "right") + 
  facet_grid(~ strategy) + 
  scale_fill_gradient(low = "gray", high = "red", "# of sites") 

sims = sims %>% 
  ggplot() + 
  geom_line(aes(year, I_hat, group = sim, color = "gray", alpha = 0.01), show.legend = F) +
  geom_line(aes(year, I), size = 2, color = "gray10") + 
  theme_minimal() + 
  ylab("biomass (g)") +
  xlab("") + 
  facet_grid(~ strategy) 
 
areas/sims
