library(sdmTMB)
library(dplyr)
library(ggplot2)
library(rgdal)

rm(list = ls())

load("data/ALL_REA_FISH_RAW.rdata")

df %>% 
  subset(REGION == "MHI") %>% 
  group_by(TAXONNAME) %>% 
  summarise(n = sum(BIOMASS_G_M2, na.rm = T)) %>% 
  mutate(freq = n/sum(n)) %>% 
  mutate(cumsum = cumsum(freq)) %>% 
  arrange(desc(freq))

df$BIOMASS_G_M2 = ifelse(df$TAXONNAME == "Aprion virescens", df$BIOMASS_G_M2, 0)

df %>% 
  subset(REGION == "MHI") %>% 
  group_by(ISLAND) %>% 
  summarise(n = sum(BIOMASS_G_M2, na.rm = T))

df = df %>% 
  # subset(REGION == "MHI") %>%
  # subset(OBS_YEAR >= 2015)
  subset(ISLAND == "Oahu")

df = df %>% 
  group_by(LONGITUDE, LATITUDE, OBS_YEAR, DEPTH) %>% 
  summarise(density = sum(BIOMASS_G_M2, na.rm = T))

qplot(df$DEPTH, df$density, color = df$density) + scale_color_viridis_c() + ggdark::dark_theme_minimal()
hist(df$density)
summary(df$density)

zone <- (floor((df$LONGITUDE[1] + 180)/6) %% 60) + 1

xy_utm = as.data.frame(cbind(utm = project(as.matrix(df[, c("LONGITUDE", "LATITUDE")]),
                                           paste0("+proj=utm +zone=", zone))))

colnames(xy_utm) = c("X", "Y"); plot(xy_utm, pch = ".", bty = 'n')

df = cbind(df, xy_utm)

rea_spde <- make_mesh(df, c("X", "Y"), cutoff = 50) # a coarse mesh for speed

plot(rea_spde, pch = "."); axis(1); axis(2)

df$year = df$OBS_YEAR
df$depth = df$DEPTH
df$depth_scaled = as.numeric(scale(log(df$depth+1)))


# default model 
m <- sdmTMB(
  data = df, 
  formula = density ~ 0 + as.factor(year) + depth,
  silent = F, 
  time = "year", 
  spde = rea_spde,
  spatial_only = F,
  family = tweedie(link = "log")
)

m <- sdmTMB(
  data = df,
  formula = density ~ 0 + logistic(depth_scaled) + as.factor(year),
  silent = F, 
  time = "year", 
  spde = rea_spde,
  family = tweedie(link = "log")
)

m <- sdmTMB(
  data = df,
  formula = density ~ 0 + as.factor(year),
  silent = F, 
  time = "year",
  spde = rea_spde, 
  family = tweedie(link = "log")
)

print(m)

predictions <- predict(m)

head(predictions)

plot(predictions$density, exp(predictions$est), xlim = c(0,1), ylim = c(0,1))

df$resids <- residuals(m) # randomized quantile residuals
hist(df$resids)

qqnorm(df$resids)
abline(a = 0, b = 1)

ggplot(df, aes(X, Y, col = resids)) + scale_colour_gradient2() +
  geom_point() + facet_wrap(~year) + coord_fixed()

nd <- data.frame(
  depth_scaled = seq(min(df$depth_scaled) + 0.2, 
                     max(df$depth_scaled) - 0.2, length.out = 100), 
  year = 2015L # a chosen year
)

p <- predict(m, newdata = nd, se_fit = TRUE, re_form = NA, xy_cols = c("X", "Y"))

ggplot(p, aes(depth_scaled, exp(est), 
              ymin = exp(est - 1.96 * est_se), 
              ymax = exp(est + 1.96 * est_se))) +
  geom_line() + geom_ribbon(alpha = 0.4)



predictions$resids <- residuals(m) # randomized quantile residuals

ggplot(predictions, aes(X, Y, col = resids)) +
  scale_colour_gradient2() +
  geom_point() +
  facet_wrap(~year)

hist(predictions$resids)

qqnorm(predictions$resids);abline(a = 0, b = 1)

m1 <- run_extra_optimization(m, nlminb_loops = 0, newton_steps = 1)
# max(m1$gradients)

plot_map <- function(dat, column) {
  
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_tile(aes(height = 500, width = 500)) +
    facet_wrap(~year) +
    coord_fixed() + 
    ggdark::dark_theme_void()
  
}

plot_map(predictions$data, "exp(est)") +
  scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Prediction (fixed effects + all random effects)")

plot_map(predictions$data, "exp(est_non_rf)") +
  ggtitle("Prediction (fixed effects only)") +
  scale_fill_viridis_c(trans = "sqrt")

plot_map(predictions$data, "omega_s") +
  ggtitle("Spatial random effects only") +
  scale_fill_gradient2()

plot_map(predictions$data, "epsilon_st") +
  ggtitle("Spatiotemporal random effects only") +
  scale_fill_gradient2()

# not bias correcting for vignette-building speed:
index <- get_index(predictions, bias_correct = FALSE)

ggplot(index, aes(year, est)) + geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  xlab('Year') + ylab('Biomass estimate (metric tonnes)')

mutate(index, cv = sqrt(exp(se^2) - 1)) %>% 
  select(-log_est, -max_gradient, -bad_eig, -se) %>%
  knitr::kable(format = "pandoc", digits = c(0, 0, 0, 0, 2))
