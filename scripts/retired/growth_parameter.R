growth_fun <- sim_vonB(Linf = 10, L0 = 1, K = 0.1, log_sd = 0.1, length_group =1, plot = TRUE, digits = 2)
growth_fun(age = rep(1:5, each = 10))
growth_fun(age = 1:15, length_age_key = T)
