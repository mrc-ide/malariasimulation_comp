# No seasonality. Bed nets. Exponential demography

source("R/0_prepare_all.R")

demog <- exp_demog
site <- no_season_site

output <- list()
for(i in seq_along(eirs)){
  print(eirs[i])
  
  ### Malaria simulation run ###################################################
  p <- p %>%
    set_bednets(
      timesteps = as.integer(c(0, 0.25 + (15:20)) * year),
      coverages = c(0, rep(0.5, 6)),
      retention = 5 * year,
      dn0 = matrix(c(0.533), nrow = 7, ncol = 3),
      rn = matrix(c(0.56), nrow = 7, ncol = 3),
      rnm = matrix(c(0.24), nrow = 7, ncol = 3),
      gamman = rep(2.64 * 365, 7))
  
  p$g0 = site$V2[4]
  p$g = site$V2[c(5, 7, 9)]
  p$h = site$V2[c(6, 8, 10)]
  
  p <- p %>% set_equilibrium(init_EIR = eirs[i])
  
  ms_raw <- run_simulation((burnin + runtime) * year, parameters = p)
  ms <- process_ms(ms_raw, burnin = burnin, pop = pop)
  ##############################################################################
  
  ### MalariaLaunchR run #######################################################
  # MalariaLaunchR run
  llin_options <- mlgts::itn_flexible_input(years = 0:9,
                                            coverage = c(rep(0, 5), rep(0.5, 5)),
                                            num_runs = 1)
  vector_options <- "itn_kill_fun 0.533 itn_repel_fun 0.56 itn_repel_min_fun 0.24 itn_kill_arab 0.533 itn_repel_arab 0.56 itn_repel_min_arab 0.24 itn_kill_gamb_ss 0.533 itn_repel_gamb_ss 0.56 itn_repel_min_gamb_ss 0.24"
  mlgts_raw <- mlgts::launch(name = "t1",
                             options = paste("num_people", pop,
                                             "prev", mean(ms$prev_2_10[1:1000]),
                                             "final_run", runtime,
                                             "output_per_yr", 365,
                                             "recalculate", 2,
                                             "prev_tol", 0.0001,
                                             "add", vector_options, llin_options),
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
ggsave("plots/comparison_noSeasonality_nets_exponentialDemog.png", summary_plot, width = 20, height = 8)


