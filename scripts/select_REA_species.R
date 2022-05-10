library(dplyr)
library(cowplot)
library(ggplot2)
library(colorRamps)

options(scipen = 999)

load("data/rea/ALL_REA_FISH_RAW.rdata")

TROPHIC_group = c("PISCIVORE", "PLANKTIVORE", "PRIMARY", "SECONDARY")

TROPHIC_group_table = NULL

for (t in 1:length(TROPHIC_group)) {
  
  # t = 1
  
  df_t = df %>% 
    subset(REGION == "MHI") %>% 
    subset(OBS_YEAR > 2009) %>% 
    subset(TROPHIC_MONREP == TROPHIC_group[t]) %>% 
    group_by(TAXONNAME) %>% 
    summarise(n = sum(BIOMASS_G_M2, na.rm = T),
              Mean_Biomass = mean(BIOMASS_G_M2, na.rm = T)) %>% 
    mutate(freq = n/sum(n)) %>% 
    mutate(Cumulative_sum = cumsum(freq),
           Total_biomass = n/1000,
           Mean_Biomass = round(Mean_Biomass/1000, 4),
           Proportion = freq, 
           Functional_group = TROPHIC_group[t])
  
  TROPHIC_group_table = rbind(TROPHIC_group_table, df_t)

}

png("/Users/kisei.tanaka/Desktop/functional_group.png", height = 10, width = 20, units = "in", res = 300)
(TROPHIC_group_table %>% 
    subset(Proportion > 0.01) %>%
    ggplot(aes(Proportion, reorder(TAXONNAME, Proportion), fill = Proportion)) +  
    geom_bar(stat = "identity") + 
    guides(color = guide_legend(), size = guide_legend()) + 
    scale_fill_gradientn(colours = matlab.like(10)) + 
    facet_wrap(.~Functional_group, ncol = 4, scales = "free") + 
    ylab("") + 
    theme_half_open() + 
    theme(axis.text.y = element_text(face = "italic")))
dev.off()
