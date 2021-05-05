d <- subset(pcod, year >= 2011) # subset for example speed
pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30) # a coarse mesh for example speed
plot(pcod_spde)

# Tweedie:
m <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
            data = d,
            time = "year", 
            spde = pcod_spde, 
            silent = F, 
            family = tweedie(link = "log"))
print(m)
tidy(m, conf.int = TRUE)
tidy(m, effects = "ran_par", conf.int = TRUE)

# Run extra optimization steps to help convergence:
m1 <- run_extra_optimization(m, nlminb_loops = 0, newton_steps = 1)
max(m$gradients)
max(m1$gradients)

# Bernoulli:
pcod_binom <- d
pcod_binom$present <- ifelse(pcod_binom$density > 0, 1L, 0L)
m_bin <- sdmTMB(present ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
                data = pcod_binom, 
                time = "year",
                silent = F, 
                spde = pcod_spde,
                family = binomial(link = "logit"))
print(m_bin)

# Gaussian:
pcod_gaus <- subset(d, density > 0 & year >= 2013)
pcod_spde_gaus <- make_mesh(pcod_gaus, c("X", "Y"), cutoff = 30)
m_pos <- sdmTMB(log(density) ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
                data = pcod_gaus, 
                silent = F, 
                time = "year", 
                spde = pcod_spde_gaus)
print(m_pos)

# With splines via mgcv.
# Make sure to pre-specify an appropriate basis dimension (`k`) since
# the smoothers are not penalized in the current implementation.
# See ?mgcv::choose.k
m_gam <- sdmTMB(log(density) ~ 0 + as.factor(year) + s(depth_scaled, k = 4),
                data = pcod_gaus, 
                silent = F, 
                time = "year", 
                spde = pcod_spde_gaus)
print(m_gam)

# With IID random intercepts:
set.seed(1)
x <- runif(500, -1, 1)
y <- runif(500, -1, 1)
loc <- data.frame(x = x, y = y)
spde <- make_mesh(loc, c("x", "y"), n_knots = 50, type = "kmeans"); plot(spde)
s <- sdmTMB_sim(x = x, y = y, 
                betas = 0, 
                time = 1L,
                phi = 0.1, 
                range = 1.4, 
                sigma_O = 0.2, 
                sigma_E = 0, 
                mesh = spde)
s$g <- gl(50, 10)
iid_re_vals <- rnorm(50, 0, 0.3)
s$observed <- s$observed + iid_re_vals[s$g]
m <- sdmTMB(observed ~ 1 + (1 | g), spde = spde, data = s)
print(m)
tidy(m, "ran_pars", conf.int = TRUE) # see tau_G
theta <- as.list(m$sd_report, "Estimate")
plot(iid_re_vals, theta$RE)


# Fit a spatial only model:
m <- sdmTMB(
  density ~ depth_scaled + depth_scaled2, data = d,
  spde = pcod_spde, 
  silent = F, 
  family = tweedie(link = "log"))
print(m)

# Spatial-trend example:
m <- sdmTMB(density ~ depth_scaled, data = d,
            spde = pcod_spde, 
            silent = F, 
            family = tweedie(link = "log"),
            spatial_trend = TRUE, 
            time = "year")
tidy(m, effects = "ran_par")

# Time-varying effects of depth and depth squared:
m <- sdmTMB(density ~ 0 + as.factor(year),
            time_varying = ~ 0 + depth_scaled + depth_scaled2,
            data = d, 
            time = "year", 
            silent = F, 
            spde = pcod_spde, 
            family = tweedie(link = "log"))
print(m)

# See the b_rw_t estimates; these are the time-varying (random walk) effects.
# These could be added to tidy.sdmTMB() eventually.
summary(m$sd_report)[1:19,]

# Linear breakpoint model on depth:
m_pos <- sdmTMB(log(density) ~ 0 + as.factor(year) +
                  breakpt(depth_scaled) + depth_scaled2, data = pcod_gaus,
                time = "year", spde = pcod_spde_gaus)
print(m_pos)

# Linear covariate on log(sigma_epsilon):
# First we will center the years around their mean
# to help with convergence.
d$year_centered <- d$year - mean(d$year)
m <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
            data = d, time = "year", spde = pcod_spde, family = tweedie(link = "log"),
            epsilon_predictor = "year_centered")
print(m) # sigma_E varies with time now
# coefficient is not yet in tidy.sdmTMB:
as.list(m$sd_report, "Estimate", report = TRUE)$b_epsilon
as.list(m$sd_report, "Std. Error", report = TRUE)$b_epsilon