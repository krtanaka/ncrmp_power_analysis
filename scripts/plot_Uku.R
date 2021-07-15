library(dplyr)
library(ggplot2)

rm(list = ls())

# for Uku: 
# Total numerical density estimates (individuals per 100 m2) were obtained by dividing fish counts in each survey by the survey area (353 m2 from two 15-m diameter survey cylinders) and multiplying by 100. - Nadon et al. 2020

region = "MHI"

islands = c("Kauai", #1
            "Lehua", #2
            "Niihau", #3
            "Kaula", #4
            "Oahu", #5
            "Molokai", #6
            "Maui", #7
            "Lanai", #8
            "Molokini", #9
            "Kahoolawe", #10
            "Hawaii")#[11]

load("data/rea/ALL_REA_FISH_RAW.rdata")

sp = "Aprion virescens"

df = df %>% 
  subset(REGION == region & ISLAND %in% islands) %>% 
  mutate(biomass = ifelse(TAXONNAME == sp, BIOMASS_G_M2, 0), 
         count = ifelse(TAXONNAME == sp, COUNT, 0)) %>%  
  group_by(LONGITUDE, LATITUDE, ISLAND, OBS_YEAR, DATE_, DEPTH) %>% 
  summarise(biomass = sum(biomass, na.rm = T),
            count = sum(count, na.rm = T))  

p1a = df %>% 
  ggplot(aes(biomass)) +
  geom_histogram() + 
  xlab("g/sq.m") + 
  ylab("") + 
  theme_minimal()  + 
  annotate("text", 
           y = Inf,
           x = Inf,
           label = "biomass",
           vjust = 2, 
           hjust = 2) 

p1b = df %>% 
  ggplot(aes(count)) +
  geom_histogram() + 
  xlab("n/353 sq.m") + 
  ylab("") + 
  theme_minimal() + 
  annotate("text", 
           y = Inf,
           x = Inf,
           label = "abundance",
           vjust = 2, 
           hjust = 2) 

p1 = p1a + p1b

p2a = df %>% 
  group_by(OBS_YEAR) %>% 
  summarise(n = mean(biomass)) %>% 
  ggplot(aes(OBS_YEAR, n)) +
  geom_line() + 
  geom_point() +
  xlab("") + 
  ylab("mean biomass g/sq.m") + 
  scale_x_continuous(breaks = c(min(df$OBS_YEAR), max(df$OBS_YEAR))) +  
  theme_pubr()

p2b = df %>% 
  group_by(OBS_YEAR) %>% 
  summarise(n = mean(count)) %>% 
  ggplot(aes(OBS_YEAR, n)) +
  geom_line() + 
  geom_point() +
  xlab("") + 
  ylab("mean abundance n/ 353 sq.m") + 
  scale_x_continuous(breaks = c(min(df$OBS_YEAR), max(df$OBS_YEAR))) +  
  theme_pubr()

p2 = p2a + p2b

biomass = df %>% 
  group_by(OBS_YEAR, ISLAND) %>% 
  summarise(n = mean(biomass)) %>% 
  mutate(n = scale(n, center = T),
         gp = "biomass")

abundance = df %>% 
  group_by(OBS_YEAR, ISLAND) %>% 
  summarise(n = mean(count)) %>% 
  mutate(n = scale(n, center = T),
         gp = "abundance")

p3 = rbind(abundance, biomass) %>% 
  ggplot(aes(OBS_YEAR, n, color = gp)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(.~ISLAND, scales = "free_y") + 
  scale_color_discrete("") + 
  xlab("") + 
  ylab("z-score") + 
  scale_x_continuous(breaks = c(min(df$OBS_YEAR), max(df$OBS_YEAR))) +  
  theme_minimal() + 
  theme(legend.position = c(0.9,0.1))

png("/Users/Kisei/Desktop/Uku_MHI_2005_2019.png", height = 5, width = 12, units = "in", res = 100)
(p1/p2) | p3
dev.off()
