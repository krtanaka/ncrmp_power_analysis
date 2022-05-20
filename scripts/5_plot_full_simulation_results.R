#######################################################################################################################
### Simulate stratified-random surveys in MHI region with reconstructed fish density (count (n) or biomass (g/sq.m) ###
### Simple Power analysis by comparing different survey "efforts"                                                   ###
#######################################################################################################################

library(SimSurvey)
library(raster)
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)
library(ggnewscale)
library(cowplot)

rm(list = ls())

select = dplyr::select

load("data/modeled_survey_variability.RData") #modeled at grid scale

# pick an island
islands = c("Hawaii", "Kauai", "Lanai", "Maui", "Molokai", "Niihau", "Oahu" )#[sample(1:7, 1)]

load(paste0('outputs/rmse_power_results_2021-09-10.RData')) # 2005-2019
load(paste0('outputs/rmse_power_results_2021-11-16.RData')) # 2010-2019

load("data/survey_effort_MHI_2014-2019.RData")

# these numbers are from "plot_REA_survey_efforts.R"
survey_effort_MHI_year = data.frame(Year = c("2013", "2016", "2019"),
                                    N = c(498, 400, 487))

low = survey_effort_MHI %>%
  select(ISLAND, low) %>% 
  rename(isl = ISLAND, 
         sites = low) %>% 
  mutate(effort = "low")

mid = survey_effort_MHI %>%
  select(ISLAND, median) %>% 
  rename(isl = ISLAND, 
         sites = median) %>% 
  mutate(effort = "mid")

high = survey_effort_MHI %>%
  select(ISLAND, high) %>% 
  rename(isl = ISLAND, 
         sites = high) %>% 
  mutate(effort = "high")

efforts = rbind(low, mid, high)

traditional_strata = isl_power %>% 
  subset(design == "traditional") %>%
  group_by(design, isl) %>% 
  summarise(traditional = mean(strata)) %>% 
  as.data.frame() %>% 
  select(isl, traditional)

downscaled_strata = isl_power %>% 
  subset(design == "downscaled") %>%
  group_by(design, isl) %>% 
  summarise(zone_based = mean(strata)) %>% 
  as.data.frame() %>% 
  select(isl, zone_based)

reduced_strata = isl_power %>% 
  subset(design == "downscaled_alt") %>%
  group_by(design, isl) %>% 
  summarise(zone_triaged = mean(strata)) %>% 
  as.data.frame() %>% 
  select(isl, zone_triaged)

strata_num = merge(downscaled_strata, traditional_strata)
strata_num = merge(strata_num, reduced_strata)

isl_power$design = ifelse(isl_power$design == "downscaled", "Zone-based", isl_power$design)
isl_power$design = ifelse(isl_power$design == "downscaled_alt", "Zone-triaged", isl_power$design)
isl_power$design = ifelse(isl_power$design == "traditional", "Traditional", isl_power$design)

isl_power = isl_power %>% 
  group_by(sp, isl) %>% 
  mutate(RMSE = as.numeric(RMSE),
         zscore = (RMSE - mean(RMSE))/sd(RMSE))

library(lemon)

(fig8 = isl_power %>%
    ggplot(aes(N, zscore)) + 
    facet_rep_grid(isl ~ sp) +
    scale_x_log10(breaks = trans_breaks('log10', function(x) 10^x),
                  labels = trans_format('log10', math_format(10^.x)),
                  "Sampling Efforts (N sites per island)") +
    geom_vline(data = efforts, aes(xintercept = sites, color = effort)) + 
    scale_color_discrete("") + 
    ggnewscale::new_scale_color() +
    ggnewscale::new_scale_fill() +
    geom_smooth(aes(color = design, fill = design), se = T) +
    scale_color_viridis_d("") + 
    scale_fill_viridis_d("") + 
    guides(color = guide_legend(override.aes = list(fill = NA))) + 
    xlab("Sampling Efforts (N sites per island)") + 
    ylab("Standadized RMSE") +
    theme(panel.background = element_rect(fill = "white"),
          # panel.grid = element_blank(),
          panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray90"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray90"),
          axis.line = element_line(),
          legend.position = "top"))

png("outputs/fig8.png", units = "in", height = 7, width = 7, res = 300)
(fig8)
dev.off()

recent_efforts = data.frame(  x = c(400, 487, 498),
                              y = c(-0.5, -0.5, -0.5),
                              label = c("2016 (n = 400)", "2019 (n = 487)", "2013 (n = 498)"))

library(ggrepel)

(fig6 = isl_power %>%
    # subset(isl == "Oahu") %>% 
    # subset(sp == "PISCIVORE") %>% 
    mutate(RMSE = as.numeric(RMSE)) %>% 
    ggplot() + 
    geom_smooth(aes(N, zscore, color = design, fill = design), show.legend = T, se = T) +
    # geom_point(aes(N, RMSE, color = design), alpha = 0.2) +
    # geom_hex(aes(N, RMSE, color = design, fill = design), alpha = 0.2, bins = 50) +
    # geom_vline(aes(xintercept = N, color = Year), data = survey_effort_MHI_year) +
    geom_label_repel(data = recent_efforts, 
                     aes(x = x, y = y, label = label), 
                     label.size = NA,
                     fontface = "bold",   
                     label.padding = 0.5, 
                     na.rm = T,
                     fill = alpha(c("white"), 0.8),
                     nudge_x = c(-0.34, -0.5, -0.1),
                     nudge_y = c(1, 3, 5)) + 
    scale_fill_viridis_d("") +
    scale_color_viridis_d("") + 
    xlab("Sampling Efforts") + 
    ylab("Standadized RMSE") + 
    # scale_y_log10() +
    # scale_x_log10() +
    scale_x_log10(breaks = trans_breaks('log10', function(x) 10^x),
                  labels = trans_format('log10', math_format(10^.x)),
                  "Sampling Efforts (N sites per year)",
                  expand = c(0, 0)) +
    # scale_y_log10(breaks = trans_breaks('log10', function(x) 10^x),
    #               labels = trans_format('log10', math_format(10^.x)), "RMSE") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_half_open() +
    guides(color = guide_legend(override.aes = list(fill = NA))) + 
    theme(legend.position = c(0, 0),
          legend.justification = c(-0.2, -0.2)))

png("outputs/fig6.png", units = "in", height = 5, width = 5, res = 500)
(fig6)
dev.off()

df = isl_power %>% 
  mutate(RMSE = as.numeric(RMSE)) %>% 
  group_by(isl, sp, design, N) %>% 
  summarise(rmse = mean(RMSE, na.rm = T))

df$rmse = formatC(df$rmse, format = "e", digits = 2)

# df$isl[duplicated(df$isl)] <- NA
# df$sp[duplicated(df$sp)] <- NA

colnames(df) = c("Island", "Trophic_group", "Survey_design", "Survey_efforts", "RMSE")

df$Survey_design = gsub("(^[[:alpha:]])", "\\U\\1", df$Survey_design, perl = T)
df$Trophic_group = tolower(df$Trophic_group)
df$Trophic_group = gsub("(^[[:alpha:]])", "\\U\\1", df$Trophic_group, perl = T)
# df[is.na(df)] <- " "

(df %>% 
    ggplot(aes(x = Island, y = RMSE, fill = Survey_design)) +
    # ggplot(aes(x = Trophic_group, y = RMSE, fill = Survey_design)) + 
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~ Trophic_group, scale = "free_y", nrow = 1) +
    # facet_grid(~ Trophic_group, scale = "free") +
    # facet_grid(~ Island, scale = "free") + 
    # scale_y_log10() +
    # scale_y_continuous(limits = quantile(df$RMSE, c(0.2, 0.9))) + 
    scale_y_log10(breaks = trans_breaks('log10', function(x) 10^x),
                  labels = trans_format('log10', math_format(10^.x)), "RMSE") +
    scale_fill_viridis_d("") + 
    xlab("") + 
    # coord_flip() + 
    # theme_pubr() + 
    theme_half_open() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)))

df$Effort_level = ""
df$Effort_level = ifelse(df$Survey_efforts %in% c(10:50), "Low effort (10-50 sites)", df$Effort_level)
df$Effort_level = ifelse(df$Survey_efforts %in% c(100:150), "Normal effort (100-150 sites)", df$Effort_level)
# df$Effort_level = ifelse(df$Survey_efforts %in% c(60:100), "Mid (60-100 sites)", df$Effort_level)
# df$Effort_level = ifelse(df$Survey_efforts %in% c(110:150), "High (110-150 sites)", df$Effort_level)

df$Effort_level = factor(df$Effort_level, levels = c('High effort (110-150 sites)', 
                                                     'Mid effort (60-100 sites)',
                                                     'Normal effort (100-150 sites)',
                                                     'Low effort (10-50 sites)'))

(fig7a = df %>% 
    subset(Effort_level %in% c("Normal effort (100-150 sites)")) %>% 
    ggplot(aes(x = Island, y = RMSE, fill = Survey_design)) +
    geom_boxplot(outlier.shape = NA, show.legend = F) +
    facet_wrap(~ Trophic_group, scale = "free", ncol = 4) +
    scale_y_log10(breaks = trans_breaks('log10', function(x) 10^x),
                  labels = trans_format('log10', math_format(10^.x)), "RMSE") +
    scale_fill_viridis_d("") + 
    xlab("") + 
    theme_half_open() +
    theme(panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray90"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray90"),
          axis.text.x = element_text(angle = 90, hjust = 1)) + 
    ggtitle("Normal effort (100-150 sites)") + 
    labs(tag = "(a)"))

(fig7b = df %>% 
    subset(Effort_level %in% c("Low effort (10-50 sites)")) %>% 
    ggplot(aes(x = Island, y = RMSE, fill = Survey_design)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~ Trophic_group, scale = "free", ncol = 4) +
    scale_y_log10(breaks = trans_breaks('log10', function(x) 10^x),
                  labels = trans_format('log10', math_format(10^.x)), "RMSE") +
    scale_fill_viridis_d("") + 
    # scale_color_viridis_d("") + 
    xlab("") + 
    theme_half_open() +
    theme(panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "gray90"),
          panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "gray90"),
          legend.position = "bottom", 
          legend.justification = c(1,1),
          axis.text.x = element_text(angle = 90, hjust = 1)) + 
    ggtitle("Low effort (10-50 sites)") + 
    labs(tag = "(b)"))

library(patchwork)
fig7 = fig7a/fig7b

png("outputs/fig7.png", units = "in", height = 10, width = 15, res = 300)
(fig7)
dev.off()

df = isl_power %>% 
  mutate(RMSE = as.numeric(RMSE)) %>% 
  group_by(
    isl,
    design,
    sp
    # N
  ) %>% 
  summarise(rmse = mean(RMSE, na.rm = T), 
            sd = sd(RMSE, na.rm = T))

df$rmse = formatC(df$rmse, format = "e", digits = 2)
df$sd = formatC(df$sd, format = "e", digits = 2)

df

library(data.table)
DT <- df %>% 
  # subset(sp == "PRIMARY") %>% 
  # subset(sp == "SECONDARY") %>% 
  # subset(sp == "PLANKTIVORE") %>% 
  # subset(sp == "PISCIVORE") %>% 
  data.table()

DT[ , .SD[which.min(rmse)], by = isl]

readr::write_csv(df, 'outputs/results.csv')

