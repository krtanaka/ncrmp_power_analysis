sim_sets_rea = function (sim, 
                         n_sims = 1, 
                         trawl_dim = c(0.01, 0.005), 
                         min_sets = 2, 
                         set_den = 2/1000, 
                         resample_cells = FALSE) {
  
  # sim = sim
  # 
  # # sim <- sim_abundance(ages = 1:5, years = 1:5) %>% sim_distribution(grid = make_grid(res = c(10, 10)))
  # 
  # n_sims = 1 
  # q = sim_logistic()
  #  trawl_dim = c(0.01, 0.005)
  # resample_cells = FALSE
  # binom_error = TRUE
  # min_sets = 2 
  # set_den = 2/1000 
  # lengths_cap = 500
  # ages_cap = 10 
  # age_sampling = "stratified" 
  # age_length_group = 1 
  # age_space_group = "division" 
  # light = TRUE
  
  strat_sets <- cell_sets <- NULL
  cells <- data.table(rasterToPoints(sim$grid))
  
  cells$strat = round(cells$strat, 0)
  
  # cells$strat = ifelse(cells$depth >= -1 & cells$depth <= 6, 1, cells$strat)
  # cells$strat = ifelse(cells$depth > 6 & cells$depth <= 18, 2, cells$strat)
  # cells$strat = ifelse(cells$depth > 18 & cells$depth <= 31, 3, cells$strat)
  
  cells
  
  strat_det <- cells[, list(strat_cells = .N), by = "strat"]; strat_det
  strat_det$tow_area <- prod(trawl_dim); strat_det
  strat_det$cell_area <- prod(res(sim$grid)); strat_det
  strat_det$strat_area <- strat_det$strat_cells * prod(res(sim$grid)); strat_det
  strat_det$strat_sets <- round(strat_det$strat_area * set_den); strat_det
  strat_det$strat_sets[strat_det$strat_sets < min_sets] <- min_sets; strat_det
  cells <- merge(cells, strat_det, by = c("strat")); cells
  
  i <- rep(seq(nrow(cells)), times = length(sim$years)); i
  y <- rep(sim$years, each = nrow(cells)); y
  
  cells <- cells[i, ]; cells
  cells$year <- y; cells
  
  i <- rep(seq(nrow(cells)), times = n_sims); i
  s <- rep(seq(n_sims), each = nrow(cells)); s
  
  cells <- cells[i, ]; cells
  cells$sim <- s; cells
  
  sets <- cells[, .SD[sample(.N, strat_sets, replace = resample_cells)], by = c("sim", "year", "strat")]; sets
  
  sets[, `:=`(cell_sets, .N), by = c("sim", "year", "cell")]; sets
  sets$set <- seq(nrow(sets)); sets
  sets
  
}

sim_survey_rea = function (sim, 
                           n_sims = 1, 
                           q = sim_logistic(),
                           trawl_dim = c(0.01, 0.005), 
                           resample_cells = FALSE,
                           binom_error = TRUE,
                           min_sets = 2, 
                           set_den = 2/1000, 
                           lengths_cap = 500,
                           ages_cap = 10, 
                           age_sampling = "stratified", 
                           age_length_group = 1, 
                           age_space_group = "division", 
                           light = TRUE) {
  
  # sim = sim
  # n_sims = 1 
  # q = sim_logistic()
  # trawl_dim = c(0.01, 0.005) 
  # resample_cells = FALSE
  # binom_error = TRUE
  # min_sets = 2 
  # set_den = 2/1000 
  # lengths_cap = 500
  # ages_cap = 10 
  # age_sampling = "stratified" 
  # age_length_group = 1 
  # age_space_group = "division" 
  # light = TRUE
  
  n <- age <- id <- division <- strat <- N <- n_measured <- n_aged <- NULL
  
  if (!age_sampling %in% c("stratified", "random")) {
    stop("age_sampling must be either \"stratified\" or \"random\". Other options have yet to be implemented.")
  }
  
  if (age_sampling == "random" && ages_cap > lengths_cap) {
    stop("When age_sampling = \"random\", ages_cap cannot exceed lengths_cap.")
  }
  
  if (!age_space_group %in% c("division", "strat", "set")) {
    stop("age_space_group must be either \"division\", \"strat\" or \"set\". Other options have yet to be implemented.")
  }
  
  sim <- round_sim(sim)
  
  I <- sim$N * q(replicate(length(sim$years), sim$ages))
  
  I_at_length <- convert_N(N_at_age = I, 
                           lak = sim$sim_length(age = sim$ages, length_age_key = TRUE))
  
  sets <- sim_sets_rea(sim, 
                       resample_cells = resample_cells, 
                       n_sims = n_sims, 
                       trawl_dim = trawl_dim, 
                       set_den = set_den, 
                       min_sets = min_sets)
  
  setkeyv(sets, c("sim", "year", "cell"))
  
  sp_I <- data.table(sim$sp_N[, c("cell", "age", "year", "N")])
  
  i <- rep(seq(nrow(sp_I)), times = n_sims)
  s <- rep(seq(n_sims), each = nrow(sp_I))
  
  sp_I <- sp_I[i, ]
  sp_I$sim <- s
  
  setdet <- merge(sets, sp_I, by = c("sim", "year", "cell"))
  
  if (binom_error) {
    
    setdet$n <- stats::rbinom(rep(1, nrow(setdet)),
                              size = round(setdet$N/setdet$cell_sets), 
                              prob = ((setdet$tow_area/setdet$cell_area) * q(setdet$age)) + 0.1) #add 0.1 bc otherwise all you get is 0
  } else {
    
    setdet$n <- round((setdet$N/setdet$cell_sets) * ((setdet$tow_area/setdet$cell_area) * q(setdet$age)))
  }
  
  setkeyv(setdet, "set")
  setkeyv(sets, "set")
  rm(sp_I)
  
  samp <- setdet[rep(seq(.N), n), list(set, age)]
  samp$id <- seq(nrow(samp))
  samp$length <- sim$sim_length(samp$age)
  measured <- samp[, list(id = id[sample(.N, ifelse(.N > lengths_cap, 
                                                    lengths_cap, .N), replace = FALSE)]), by = "set"]
  samp$measured <- samp$id %in% measured$id
  length_samp <- samp[samp$measured, ]
  rm(measured)
  
  length_samp$length_group <- group_lengths(length_samp$length, 
                                            age_length_group)
  
  length_samp <- merge(sets[, list(set, sim, year, division, strat)], 
                       length_samp, 
                       by = "set")
  
  if (age_sampling == "stratified") {
    
    aged <- length_samp[, 
                        list(id = id[sample(.N, ifelse(.N > ages_cap, ages_cap, .N), replace = FALSE)]), 
                        by = c("sim", "year", age_space_group, "length_group")]
    
  }
  
  if (age_sampling == "random") {
    
    aged <- length_samp[, 
                        list(id = id[sample(.N, ifelse(.N > ages_cap, ages_cap, .N), replace = FALSE)]), 
                        by = c("set")]
    
  }
  
  samp$aged <- samp$id %in% aged$id
  rm(aged)
  rm(length_samp)
  samp <- samp[, list(set, id, length, age, measured, aged)]
  if (light) 
    samp$id <- NULL
  if (!light) 
    full_setdet <- setdet
  setdet <- merge(sets, setdet[, list(N = sum(N), n = sum(n)), 
                               by = "set"], by = "set")
  setdet <- merge(setdet, samp[, list(n_measured = sum(measured), 
                                      n_aged = sum(aged)), by = "set"], by = "set", 
                  all.x = TRUE)
  setdet$n_measured[is.na(setdet$n_measured)] <- 0
  setdet$n_aged[is.na(setdet$n_aged)] <- 0
  samp_totals <- setdet[, list(n_sets = .N, n_caught = sum(n), 
                               n_measured = sum(n_measured), n_aged = sum(n_aged)), 
                        by = c("sim", "year")]
  sim$I <- I
  sim$I_at_length <- I_at_length
  if (!light) 
    sim$full_setdet <- full_setdet
  sim$setdet <- setdet
  sim$samp <- samp
  sim$samp_totals <- samp_totals
  sim
}