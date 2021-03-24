library(SimSurvey)

#Simulate starting abundance, random recruitment and total mortality
R_fun <- sim_R(log_mean = log(100000), 
               log_sd = 0.1,
               random_walk = TRUE, #Simulate recruitment as a random walk
               plot = TRUE)

R_fun(years = 1:100)


#generate N0
sim_abundance(R = sim_R(log_mean = log(100000), log_sd = 0.5))

sim_abundance(years = 1:20,
              R = sim_R(log_mean = log(c(rep(100000, 10), rep(10000, 10))), plot = TRUE))

#generate Z values
Z_fun <- sim_Z(log_mean = log(0.5), 
               log_sd = 0.1, 
               phi_age = 0.9, #Autoregressive parameter for the age dimension.
               phi_year = 0.9, #Autoregressive parameter for the year dimension.
               plot = TRUE)

Z_fun(years = 1:100, ages = 1:20)

#generate N0
sim_abundance(Z = sim_Z(log_mean = log(0.5), 
                        log_sd = 0.1,
                        plot = TRUE))


Za_dev <- c(-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.3, 0.2, 0.1, 0) #age
Zy_dev <- c(-0.2, -0.2, -0.2, -0.2, -0.2, 2, 2, 2, 2, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0, 0, 0, 0, 0) #year

Z_mat <- outer(Za_dev, Zy_dev, "+") + 0.5

#Simulate basic population dynamics model
sim_abundance(ages = 1:10, 
              years = 1:20,
              Z = sim_Z(log_mean = log(Z_mat), # generating mortality matrix
                        plot = TRUE))

sim_abundance(ages = 1:10, 
              years = 1:20,
              Z = sim_Z(log_mean = log(Z_mat), 
                        log_sd = 0.5,
                        phi_age = 0.5, 
                        phi_year = 0.5,
                        plot = TRUE))

#generate N0
N0_fun <- sim_N0(N0 = "exp", plot = TRUE)
N0_fun(R0 = 1000, Z0 = rep(0.5, 20), ages = 1:20)

sim_abundance(N0 = sim_N0(N0 = "exp", plot = TRUE))

#Closure for simulating length given age using von Bertalanffy notation
growth_fun <- sim_vonB(Linf = 100, 
                       L0 = 5, 
                       K = 0.2, 
                       log_sd = 0.05, 
                       length_group = 1,
                       plot = TRUE)

growth_fun(age = rep(1:15, each = 100))

growth_fun(age = 1:15, length_age_key = TRUE)

sim_abundance(growth = sim_vonB(plot = TRUE))

sim <- sim_abundance()
plot_trend(sim)
plot_surface(sim, mat = "N")
plot_surface(sim, mat = "Z")
plot_surface(sim, mat = "N_at_length", xlab = "Length", zlab = "N")

