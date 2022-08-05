library(patchwork)
library(ggplot2)
library(dplyr)
library(cowplot)
library(colorRamps)

rm(list = ls())

areas = NULL
sims = NULL

for (i in 1:3) {
  
  # i = 1
  
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
  rm(area, sim, sim_results)
  
}

(rmse = sims %>% 
    subset(year >= 2010) %>%
    group_by(strategy) %>% 
    summarise(rmse = sqrt(mean(error^2))) %>% 
    mutate(rmse = formatC(rmse, digits = 2)))

(num_strata = areas %>% 
    group_by(strategy) %>% 
    summarise(sum = length(unique(strat))))

(f4aa = areas %>% 
    subset(strategy == "Traditional") %>% 
    ggplot(aes(x, y)) +
    coord_fixed() +
    geom_raster(aes(fill = factor(strat)), show.legend = F) + 
    scale_fill_viridis_d("") +
    theme_classic() +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = c(0.01, 0.01)) +
    annotate("text",  x = Inf, y = Inf, label = "n = 21", vjust = 1.2, hjust = 1.2, size = 4) + 
    theme(panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray80"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray80"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(face = "bold")) + 
    ggtitle("Traditional") + 
    labs(tag = "(a)"))

(f4ab = areas %>% 
    subset(strategy == "Zone-based") %>% 
    ggplot(aes(x, y)) +
    coord_fixed() +
    geom_raster(aes(fill = factor(strat)), show.legend = F) + 
    scale_fill_viridis_d("") +
    theme_classic() +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = c(0.01, 0.01)) +
    annotate("text",  x = Inf, y = Inf, label = "n = 9", vjust = 1.2, hjust = 1.2, size = 4) + 
    theme(panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray80"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray80"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(face = "bold")) + 
    ggtitle("Zone-based") + 
    labs(tag = " "))

(f4ac = areas %>% 
    subset(strategy == "Zone-triaged") %>% 
    ggplot(aes(x, y)) +
    coord_fixed() +
    geom_raster(aes(fill = factor(strat)), show.legend = F) + 
    scale_fill_viridis_d("") +
    theme_classic() +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = c(0.01, 0.01)) +
    annotate("text",  x = Inf, y = Inf, label = "n = 3", vjust = 1.2, hjust = 1.2, size = 4) + 
    theme(panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray80"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray80"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(face = "bold")) + 
    ggtitle("Zone-triaged") + 
    labs(tag = " "))

f4a = (f4aa + f4ab + f4ac); rm(f4aa, f4ab, f4ac)

(f4ba = areas %>% 
    subset(strategy == "Traditional") %>% 
    ggplot(aes(x, y)) +
    coord_fixed() +
    geom_raster(aes(fill = strat_sets)) + 
    theme_classic() +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = c(0.01, 0.01)) +
    theme(legend.position = "right") +
    # scale_fill_viridis_c("") +
    scale_fill_gradient2(mid = "blue", high = "red", "") + 
    theme(legend.position = c(1,1),
          legend.justification = c(1.1, 0.85),
          legend.key.size = unit(0.3, "cm"),
          panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray80"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray80"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(face = "bold")) + 
    ggtitle("Traditional") + 
    labs(tag = "(b)"))

(f4bb = areas %>% 
    subset(strategy == "Zone-based") %>% 
    ggplot(aes(x, y)) +
    coord_fixed() +
    geom_raster(aes(fill = strat_sets)) + 
    theme_classic() +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = c(0.01, 0.01)) +
    theme(legend.position = "right") +
    # scale_fill_viridis_c("") +
    scale_fill_gradient2(mid = "blue", high = "red", "") + 
    theme(legend.position = c(1,1),
          legend.justification = c(1.1, 0.85),
          legend.key.size = unit(0.3, "cm"),
          panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray80"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray80"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(face = "bold")) + 
    ggtitle("Zone-based") + 
    labs(tag = " "))

(f4bc = areas %>% 
    subset(strategy == "Zone-triaged") %>% 
    ggplot(aes(x, y)) +
    coord_fixed() +
    geom_raster(aes(fill = strat_sets)) + 
    theme_classic() +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = c(0.01, 0.01)) +
    theme(legend.position = "right") +
    # scale_fill_viridis_c("") +
    scale_fill_gradient2(mid = "blue", high = "red", "") + 
    theme(legend.position = c(1,1),
          legend.justification = c(1.1, 0.85),
          legend.key.size = unit(0.3, "cm"),
          panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray80"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray80"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(face = "bold")) + 
    ggtitle("Zone-triaged") + 
    labs(tag = " "))

f4b = (f4ba + f4bb + f4bc); rm(f4ba, f4bb, f4bc)

f4a / f4b

png("outputs/fig4.png", units = "in", height = 6, width = 10, res = 500)
(f4a / f4b)
dev.off()

(fig5 = sims %>% 
    mutate(I_hat = I_hat/1000,
           I = I/1000) %>% 
    subset(year >= 2010) %>%
    ggplot() + 
    geom_point(aes(year, I_hat, group = sim, color = sim), 
               alpha = 0.5, 
               size = 2, 
               position = position_jitter(0.3),
               show.legend = T) +
    scale_color_viridis_c("Simulation") + 
    ggnewscale::new_scale_color() +
    geom_line(aes(year, I, color = "True biomass"), size = 1) + 
    geom_point(aes(year, I, color = "True biomass"), size = 2) + 
    scale_color_viridis_d("", option = "A") +
    theme_half_open() +
    scale_x_continuous(breaks = c(2010, 2012, 2013, 2015, 2016, 2019), 
                       labels = c(2010, 2012, 2013, 2015, 2016, 2019)) +
    ylab("Total Biomass (kg)") +
    xlab("") + 
    facet_grid(~ strategy) + 
    geom_text(
      data = rmse,
      mapping = aes(x = Inf, y = Inf, label = paste0("RMSE = ", rmse)),
      hjust = 1,
      vjust = 2) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))

png("outputs/fig5.png", units = "in", height = 6, width = 12, res = 500)
(fig5)
dev.off()

