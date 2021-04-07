# plot_survey = function (sim, which_year = 1, which_sim = 1) {
#   
#   xax <- list(title = "", 
#               zeroline = FALSE, 
#               showline = FALSE, 
#               showticklabels = FALSE, 
#               showgrid = FALSE,
#               ticks = "")
#   
#   yax <- c(scaleanchor = "x", xax)
#   
#   sp_strat <- raster::rasterToPolygons(sim$grid$strat, dissolve = TRUE)
#   
#   df_strat <- suppressMessages(ggplot2::fortify(sp_strat) %>% group_by(.data$group))
#   
#   setdet <- sim$setdet
#   setdet <- setdet[setdet$year == which_year & setdet$sim == which_sim, ]
#   
#   samp <- sim$samp
#   samp <- samp[samp$set %in% setdet$set, ]
#   samp <- merge(setdet[, c("year", "sim", "set", "x", "y")], samp, by = "set", all = TRUE)
#   
#   l <- sim$I_at_length[, as.character(which_year)]
#   l <- (l/sum(l)) * 100
#   
#   true_length <- data.frame(length = as.numeric(names(l)), 
#                             percent = l)
#   a <- sim$I[, as.character(which_year)]
#   a <- (a/sum(a)) * 100
#   
#   true_age <- data.frame(age = as.numeric(names(a)), percent = a)
#   length_group <- get("length_group", envir = environment(sim$sim_length))
#   
#   d <- crosstalk::SharedData$new(samp, ~set)
#   
#   base <- plot_ly(data = d)
#   
#   sp_p <- base %>% 
#     group_by(set) %>% 
#     summarise(x = unique(.data$x), 
#               y = unique(.data$y), 
#               n = sum(!is.na(.data$measured))) %>% 
#     add_markers(x = ~x, 
#                 y = ~y, 
#                 text = ~n, 
#                 color = ~n, 
#                 name = "n", 
#                 showlegend = FALSE,
#                 marker = list(size = ~.scale_between(n, 2, 600), sizemode = "area")) %>% 
#     add_paths(data = df_strat,
#               x = ~long, 
#               y = ~lat, color = I("black"), 
#               hoverinfo = "none", 
#               size = I(0.5), showlegend = FALSE, alpha = 0.1) %>%
#     layout(xaxis = xax, 
#            yaxis = yax, 
#            margin = list(t = 0, r = 0, l = 0, b = 0, pad = 0))
#   
#   hist_base <- base %>%
#     mutate(length = group_lengths(length, length_group)) %>%
#     group_by(set) %>% 
#     filter(!is.na(.data$measured)) %>% 
#     slice(rep(1:length(.data$measured), each = 2)) %>% 
#     mutate(lab = rep(c("caught", "sampled"), times = length(.data$measured)/2))
#   
#   lf_p <- hist_base %>% 
#     filter(.data$measured & .data$lab == "sampled" | .data$lab == "caught") %>%
#     add_histogram(x = ~length, 
#                   color = ~lab, 
#                   histnorm = "percent", 
#                   colors = c("#FDE725FF", "#21908CFF"), 
#                   legendgroup = ~lab) %>% 
#     add_lines(data = true_length, 
#               x = ~length, 
#               y = ~percent, 
#               color = I("#440154FF"), 
#               fill = "tozeroy", 
#               name = "true", 
#               fillcolor = "#44015433", 
#               legendgroup = "true") %>% 
#     layout(xaxis = list(title = "Length", 
#                         range = grDevices::extendrange(range(samp$length, na.rm = TRUE))), 
#            yaxis = list(title = "Percent", ticksuffix = "%"))
#   
#   af_p <- hist_base %>% 
#     filter(.data$aged & 
#              .data$lab == "sampled" | 
#              .data$lab == "caught") %>% 
#     add_histogram(x = ~age, 
#                   color = ~lab, 
#                   histnorm = "percent",
#                   colors = c("#FDE725FF", "#21908CFF"), 
#                   legendgroup = ~lab, 
#                   showlegend = FALSE) %>% 
#     add_lines(data = true_age, 
#               x = ~age, 
#               y = ~percent, 
#               color = I("#440154FF"), 
#               fill = "tozeroy", 
#               name = "true",
#               fillcolor = "#44015433", 
#               legendgroup = "true", 
#               showlegend = FALSE) %>% 
#     layout(xaxis = list(title = "Age", 
#                         range = grDevices::extendrange(range(samp$age, na.rm = TRUE))), 
#            yaxis = list(title = "Percent",
#                         ticksuffix = "%"))
#   
#   subplot(sp_p, subplot(lf_p, 
#                         af_p, 
#                         nrows = 2, 
#                         titleX = TRUE, 
#                         titleY = TRUE,
#                         margin = 0.1), 
#           nrows = 1, 
#           margin = c(0, 0.2, 0, 0), 
#           widths = c(0.75, 0.25), 
#           titleX = TRUE, 
#           titleY = TRUE) %>% 
#     colorbar(x = 0.52, y = 1) %>% 
#     layout(legend = list(x = 1, y = 0.95))
# }

sim_sets_rea = function (sim, 
                         n_sims = 1, 
                         trawl_dim = c(0.01, 0.005), 
                         min_sets = 2, 
                         set_den = 2/1000, 
                         resample_cells = FALSE) {
  
  strat_sets <- cell_sets <- NULL
  cells <- data.table(rasterToPoints(sim$grid))
  
  cells$strat = ifelse(cells$depth >= 0 & cells$depth <= 6, 1, cells$strat)
  cells$strat = ifelse(cells$depth > 6 & cells$depth <= 18, 2, cells$strat)
  cells$strat = ifelse(cells$depth > 18 & cells$depth <= 30, 3, cells$strat)
  
  strat_det <- cells[, list(strat_cells = .N), by = "strat"]
  strat_det$tow_area <- prod(trawl_dim)
  strat_det$cell_area <- prod(res(sim$grid))
  strat_det$strat_area <- strat_det$strat_cells * prod(res(sim$grid))
  strat_det$strat_sets <- round(strat_det$strat_area * set_den)
  strat_det$strat_sets[strat_det$strat_sets < min_sets] <- min_sets
  cells <- merge(cells, strat_det, by = c("strat"))
  i <- rep(seq(nrow(cells)), times = length(sim$years))
  y <- rep(sim$years, each = nrow(cells))
  cells <- cells[i, ]
  cells$year <- y
  i <- rep(seq(nrow(cells)), times = n_sims)
  s <- rep(seq(n_sims), each = nrow(cells))
  cells <- cells[i, ]
  cells$sim <- s
  sets <- cells[, .SD[sample(.N, strat_sets, replace = resample_cells)], 
                by = c("sim", "year", "strat")]
  
  sets[, `:=`(cell_sets, .N), by = c("sim", "year", "cell")]
  sets$set <- seq(nrow(sets))
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