# No seasonality. Treatment. Exponential demography

source("R/0_prepare_all.R")

demog <- exp_demog
site <- no_season_site

output <- list()
for(i in seq_along(eirs)){
  print(eirs[i])
  
  ### Malaria simulation run ###################################################
  p <- p %>%
    set_drugs(list(AL_params)) %>%
    set_clinical_treatment(1, c(0, 15) * year, c(0, 0.8))
  
  p$g0 = site$V2[4]
  p$g = site$V2[c(5, 7, 9)]
  p$h = site$V2[c(6, 8, 10)]
  
  p <- p %>% set_equilibrium(init_EIR = eirs[i])
  
  ms_raw <- run_simulation((burnin + runtime) * year, parameters = p)
  ms <- process_ms(ms_raw, burnin = burnin, pop = pop)
  ##############################################################################
  
  ### MalariaLaunchR run #######################################################
  tx_options <- 
    mlgts::treat_flexible_input(0:9, 
                                drug_0_coverage = rep(0, 10),
                                drug_1_coverage = c(rep(0, 5), rep(0.8, 5)),
                                drug_2_coverage = rep(0, 10),
                                drug_3_coverage = rep(0, 10))
  mlgts_raw <- mlgts::launch(name = "t1",
                             options = paste("num_people", pop,
                                             "prev", mean(ms$prev_2_10[1:1000]),
                                             "final_run", runtime,
                                             "output_per_yr", 365,
                                             "recalculate", 2,
                                             "prev_tol", 0.0001,
                                             "add", tx_options),
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
ggsave("plots/comparison_noSeasonality_treatment_exponentialDemog.png", summary_plot, width = 20, height = 8)




