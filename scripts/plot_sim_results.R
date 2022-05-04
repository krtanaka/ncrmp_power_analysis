library(patchwork)
library(ggplot2)
library(dplyr)
library(cowplot)
library(colorRamps)

rm(list = ls())

areas = NULL
sims = NULL

for (i in 1:3) {
  
  # i = 2
  
  # # 2005-2019
  # if (i == 1) load('outputs/sim_results_Maui_traditional_median_PLANKTIVORE_100_biomass.RData')
  # if (i == 2) load('outputs/sim_results_Maui_downscaled_median_PLANKTIVORE_100_biomass.RData')
  # if (i == 3) load('outputs/sim_results_Maui_downscaled_alt_median_PLANKTIVORE_100_biomass.RData')
  
  # # 2010-2019
  # if (i == 1) load('outputs/sim_results_Maui_traditional_median_PISCIVORE_100_biomass.RData')
  # if (i == 2) load('outputs/sim_results_Maui_downscaled_median_PISCIVORE_100_biomass.RData')
  # if (i == 3) load('outputs/sim_results_Maui_downscaled_alt_median_PISCIVORE_100_biomass.RData')
  
  # 2005-2019
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
    # subset(year >= 2010) %>%
    group_by(strategy) %>% 
    summarise(rmse = sqrt(mean(error^2))) %>% 
    mutate(rmse = formatC(rmse, digits = 2)))

(p1 = areas %>% 
    ggplot(aes(x, y)) +
    coord_fixed() +
    geom_raster(aes(fill = log(strat_sets))) + 
    theme_half_open() + 
    ylab("Northing (km)") + xlab("Easting (km)") + 
    theme(legend.position = "right") +
    facet_grid(~ strategy) + 
    scale_fill_gradientn(colours = matlab.like(100), "log(# of sites)") +
    theme(#axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          # panel.background = element_rect(fill = "gray10", colour = "gray10"),
          # panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "gray20"), 
          # panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "gray20")
          ) + 
    ggtitle("(a)"))

(p2 = sims %>% 
    subset(year >= 2010) %>%
    ggplot() + 
    geom_point(aes(year, I_hat, group = sim, color = sim), 
               alpha = 0.5, 
               size = 2, 
               position = position_jitter(0.3),
               show.legend = T) +
    # geom_line(aes(year, I_hat, group = sim, color = sim), alpha = 0.5, size = 1, show.legend = T) +
    scale_color_viridis_c("Simulation") + 
    scale_color_gradientn(colours = matlab.like(100), "Simulation") +
    ggnewscale::new_scale_color() +
    geom_line(aes(year, I, color = "True biomass"), size = 1) + 
    geom_point(aes(year, I, color = "True biomass"), size = 2) + 
    scale_color_viridis_d("", option = "A") +
    # scale_color_discrete("") +
    theme_half_open() + 
    # scale_y_log10() +
    # scale_x_log10() + 
    scale_x_continuous(breaks = c(2010, 2012, 2013, 2015, 2016, 2019), labels = c(2010, 2012, 2013, 2015, 2016, 2019)) +
    ylab("Biomass (g)") +
    xlab("") + 
    facet_grid(~ strategy) + 
    geom_text(
      data = rmse,
      mapping = aes(x = Inf, y = Inf, label = paste0("RMSE = ", rmse)),
      hjust = 1,
      vjust = 2) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ggtitle("(b)"))

p1/p2

png("outputs/fig5.png", units = "in", height = 8, width = 13, res = 500)
(p1/p2)
dev.off()
