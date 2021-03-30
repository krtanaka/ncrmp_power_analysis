plot_survey = function (sim, which_year = 1, which_sim = 1) {
  
  xax <- list(title = "", 
              zeroline = FALSE, 
              showline = FALSE, 
              showticklabels = FALSE, 
              showgrid = FALSE,
              ticks = "")
  
  yax <- c(scaleanchor = "x", xax)
  
  sp_strat <- raster::rasterToPolygons(sim$grid$strat, dissolve = TRUE)
  
  df_strat <- suppressMessages(ggplot2::fortify(sp_strat) %>% group_by(.data$group))
  
  setdet <- sim$setdet
  setdet <- setdet[setdet$year == which_year & setdet$sim == which_sim, ]
  
  samp <- sim$samp
  samp <- samp[samp$set %in% setdet$set, ]
  samp <- merge(setdet[, c("year", "sim", "set", "x", "y")], samp, by = "set", all = TRUE)
  
  l <- sim$I_at_length[, as.character(which_year)]
  l <- (l/sum(l)) * 100
  
  true_length <- data.frame(length = as.numeric(names(l)), 
                            percent = l)
  a <- sim$I[, as.character(which_year)]
  a <- (a/sum(a)) * 100
  
  true_age <- data.frame(age = as.numeric(names(a)), percent = a)
  length_group <- get("length_group", envir = environment(sim$sim_length))
  
  d <- crosstalk::SharedData$new(samp, ~set)
  
  base <- plot_ly(data = d)
  
  sp_p <- base %>% 
    group_by(set) %>% 
    summarise(x = unique(.data$x), 
              y = unique(.data$y), 
              n = sum(!is.na(.data$measured))) %>% 
    add_markers(x = ~x, 
                y = ~y, 
                text = ~n, 
                color = ~n, 
                name = "n", 
                showlegend = FALSE,
                marker = list(size = ~.scale_between(n, 2, 600), sizemode = "area")) %>% 
    add_paths(data = df_strat,
              x = ~long, 
              y = ~lat, color = I("black"), 
              hoverinfo = "none", 
              size = I(0.5), showlegend = FALSE, alpha = 0.1) %>%
    layout(xaxis = xax, 
           yaxis = yax, 
           margin = list(t = 0, r = 0, l = 0, b = 0, pad = 0))
  
  hist_base <- base %>%
    mutate(length = group_lengths(length, length_group)) %>%
    group_by(set) %>% 
    filter(!is.na(.data$measured)) %>% 
    slice(rep(1:length(.data$measured), each = 2)) %>% 
    mutate(lab = rep(c("caught", "sampled"), times = length(.data$measured)/2))
  
  lf_p <- hist_base %>% 
    filter(.data$measured & .data$lab == "sampled" | .data$lab == "caught") %>%
    add_histogram(x = ~length, 
                  color = ~lab, 
                  histnorm = "percent", 
                  colors = c("#FDE725FF", "#21908CFF"), 
                  legendgroup = ~lab) %>% 
    add_lines(data = true_length, 
              x = ~length, 
              y = ~percent, 
              color = I("#440154FF"), 
              fill = "tozeroy", 
              name = "true", 
              fillcolor = "#44015433", 
              legendgroup = "true") %>% 
    layout(xaxis = list(title = "Length", 
                        range = grDevices::extendrange(range(samp$length, na.rm = TRUE))), 
           yaxis = list(title = "Percent", ticksuffix = "%"))
  
  af_p <- hist_base %>% 
    filter(.data$aged & 
             .data$lab == "sampled" | 
             .data$lab == "caught") %>% 
    add_histogram(x = ~age, 
                  color = ~lab, 
                  histnorm = "percent",
                  colors = c("#FDE725FF", "#21908CFF"), 
                  legendgroup = ~lab, 
                  showlegend = FALSE) %>% 
    add_lines(data = true_age, 
              x = ~age, 
              y = ~percent, 
              color = I("#440154FF"), 
              fill = "tozeroy", 
              name = "true",
              fillcolor = "#44015433", 
              legendgroup = "true", 
              showlegend = FALSE) %>% 
    layout(xaxis = list(title = "Age", 
                        range = grDevices::extendrange(range(samp$age, na.rm = TRUE))), 
           yaxis = list(title = "Percent",
                        ticksuffix = "%"))
  
  subplot(sp_p, subplot(lf_p, 
                        af_p, 
                        nrows = 2, 
                        titleX = TRUE, 
                        titleY = TRUE,
                        margin = 0.1), 
          nrows = 1, 
          margin = c(0, 0.2, 0, 0), 
          widths = c(0.75, 0.25), 
          titleX = TRUE, 
          titleY = TRUE) %>% 
    colorbar(x = 0.52, y = 1) %>% 
    layout(legend = list(x = 1, y = 0.95))
}
