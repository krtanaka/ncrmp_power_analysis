library(dplyr)

load("data/ALL_REA_FISH_RAW.rdata")

df %>% 
  group_by(TAXONNAME) %>% 
  summarise(n = sum(BIOMASS_G_M2, na.rm = T)) %>% 
  mutate(freq = n/sum(n)) %>% 
  arrange(desc(freq)) %>% 
  mutate(cumsum = cumsum(freq)) %>% 
  subset(cumsum < 0.5)

