process_ms <- function(ms_raw, burnin, pop) {
  ms_raw %>%
    mutate(
      EIR = 365 * (EIR_arab + EIR_fun + EIR_gamb) / pop,
      prev_2_10 = n_detect_730_3650  / n_730_3650,
      clin = (p_inc_clinical_0_36500 / n_0_36500) * 365,
      clin_smooth = roll_mean(clin, n = 365, align = "right", fill = NA),
      sev = (p_inc_severe_0_36500 / n_0_36500)  * 365,
      sev_smooth = roll_mean(sev, n = 365, align = "right", fill = NA, na.rm = TRUE),
      model = "malariasimulation") %>%
    filter(timestep > burnin * year) %>%
    mutate(timestep = timestep - burnin * year) %>%
    select(model, timestep, EIR, prev_2_10, clin, clin_smooth, sev, sev_smooth)
}


process_mlgts <- function(mlgts_raw) {
  mlgts_raw$output %>%
    rename(
      EIR = EIR,
      clin = clin_inc_all,
      clin_smooth = clin_inc_all_smooth,
      sev = sev_inc_all,
      sev_smooth = sev_inc_all_smooth) %>%
    mutate(
      clin_smooth = clin_smooth,
      sev_smooth = sev_smooth,
      timestep = year * 365,
      model = "old") %>%
    select(model, timestep, EIR, prev_2_10, clin, clin_smooth, sev, sev_smooth)
}

make_summary_plot <- function(output){
  output <- output %>%
    mutate(model = factor(model, levels = c("old", "malariasimulation")))
  
  eir_plot <- ggplot(output, aes(x = timestep, y = EIR, col = model)) +
    geom_line() +
    ylab("EIR") + 
    theme_bw() +
    facet_wrap( ~ eir, nrow = 1, scales = "free_y")
  prev_plot <- ggplot(output, aes(x = timestep, y = prev_2_10, col = model)) +
    geom_line() +
    ylab("PfPr 2-10") + 
    theme_bw() +
    facet_wrap( ~ eir, nrow = 1, scales = "free_y")
  inc_plot <- ggplot(output, aes(x = timestep, y = clin_smooth, col = model)) +
    geom_line() +
    ylab("Smooth annual\nclinical incidence") + 
    theme_bw() +
    facet_wrap( ~ eir, nrow = 1, scales = "free_y")
  sev_plot <- ggplot(output, aes(x = timestep, y = sev_smooth, col = model)) +
    geom_line() +
    ylab("Smooth annual\nsevere incidence") + 
    theme_bw() +
    facet_wrap( ~ eir, nrow = 1, scales = "free_y")
  
 p1 <- (prev_plot / eir_plot / inc_plot / sev_plot) #+ plot_layout(guide = "collect")
 p2 <- ((prev_plot + ylim(0, 1)) / (eir_plot + ylim(0, 100)) / (inc_plot + ylim(0, 1)) / (sev_plot + ylim(0, 0.02)))# + plot_layout(guide = "collect")

 (p2 | p1)  + plot_layout(guide = "collect")
}
