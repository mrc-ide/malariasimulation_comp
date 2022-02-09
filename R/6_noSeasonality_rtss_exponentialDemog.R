# No seasonality. RTSS. Exponential demography

source("R/0_prepare_all.R")

demog <- exp_demog
site <- no_season_site

output <- list()
for(i in seq_along(eirs)){
  print(eirs[i])
  
  ### Malaria simulation run ###################################################
  p$rtss <- TRUE
  p$rtss_doses <- round(c(0, 1.5 * month, 3 * month))
  
  p <- p %>%
    set_rtss_epi(
      start = 15 * year,
      end = 20 * year,
      age = 5 * month,
      coverage = 0.8,
      min_wait = 0,
      boosters = 18 * month,
      booster_coverage = 0.7)  
  
  p$g0 = site$V2[4]
  p$g = site$V2[c(5, 7, 9)]
  p$h = site$V2[c(6, 8, 10)]
  
  p <- p %>% set_equilibrium(init_EIR = eirs[i])
  
  ms_raw <- run_simulation((burnin + runtime) * year, parameters = p)
  ms <- process_ms(ms_raw, burnin = burnin, pop = pop)
  ##############################################################################
  
  ### MalariaLaunchR run #######################################################
  rtss_options <- mlgts::epi_pev_flexible_input(0:9, c(rep(0, 5), rep(0.8, 5)))
  mlgts_raw <- mlgts::launch(name = "t1",
                             options = paste("num_people", pop,
                                             "prev", mean(ms$prev_2_10[1:1000]),
                                             "final_run", runtime,
                                             "output_per_yr", 365,
                                             "recalculate", 2,
                                             "prev_tol", 0.0001,
                                             "add", rtss_options, "ab_model 1"),
                             demog = demog,
                             site = site)
  
  mlgts <- process_mlgts(mlgts_raw)
  ##############################################################################
  
  # Combine output
  output[[i]] <- bind_rows(ms, mlgts) %>%
    mutate(eir = eirs[i])
}
output <- bind_rows(output)

summary_plot <- make_summary_plot(output)
ggsave("plots/comparison_noSeasonality_rtss_exponentialDemog.png", width = 20, height = 8)

