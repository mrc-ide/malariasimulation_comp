# No seasonality. No interventions. Custom demography

source("R/0_prepare_all.R")

demog <- flat_demog
site <- no_season_site

output <- list()
for(i in seq_along(eirs)){
  print(eirs[i])
  
  ### Malaria simulation run ###################################################
  p <- p %>%
    set_demography(
      agegroups = round(demog$V3 * year),
      timesteps = 1,
      deathrates = matrix(demog$V5 / 365, nrow = 1),
      birthrates = find_birthrates(pop, round(demog$V3 * year), demog$V5 / 365)
    )
  
  p$g0 = site$V2[4]
  p$g = site$V2[c(5, 7, 9)]
  p$h = site$V2[c(6, 8, 10)]
  
  p <- p %>% set_equilibrium(init_EIR = eirs[i])
  
  ms_raw <- run_simulation((burnin + runtime) * year, parameters = p)
  ms <- process_ms(ms_raw, burnin = burnin, pop = pop)
  ##############################################################################
  
  ### MalariaLaunchR run #######################################################
  mlgts_raw <- mlgts::launch(name = "t1",
                             options = paste("num_people", pop,
                                             "prev", mean(ms$prev_2_10),
                                             "final_run", runtime,
                                             "output_per_yr", 365,
                                             "recalculate", 2,
                                             "prev_tol", 0.0001),
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
ggsave("plots/comparison_noSeasonality_noInterventions_customDemog.png", summary_plot, width = 20, height = 8)
