# Install severe branch
# remotes::install_github("mrc-ide/individual")
# remotes::install_github("mrc-ide/malariasimulation@feat/simple_severe")

# Load packages
library("malariasimulation")
library(dplyr)
library(ggplot2)
library(patchwork)
library(RcppRoll)
library(mlgts)

# General set up
year <- 365
pop <- 10000
eirs <- c(0.1, 1, 5, 10, 30)
demog <- read.table("inputs/exponential_demog.txt", sep = "\t")

output <- list()
for(i in seq_along(eirs)){
  print(eirs[i])
  ## Malaria simulation run
  p <- get_parameters(list(
    human_population = pop,
    prevalence_rendering_min_ages = 2 * year,
    prevalence_rendering_max_ages = 10 * year,
    clinical_incidence_rendering_min_ages = 0 * year,
    clinical_incidence_rendering_max_ages = 100 * year,
    severe_incidence_rendering_min_ages = 0 * year,
    severe_incidence_rendering_max_ages = 100 * year,
    individual_mosquitoes = FALSE
  )) %>%
    set_bednets(
      timesteps = as.integer(c(0, 0.25 + (15:20)) * year),
      coverages = c(0, rep(0.5, 6)),
      retention = 5 * year,
      dn0 = matrix(c(0.533), nrow = 7, ncol = 1),
      rn = matrix(c(0.56), nrow = 7, ncol = 1),
      rnm = matrix(c(0.24), nrow = 7, ncol = 1),
      gamman = rep(2.64 * 365, 7)) %>%
    set_equilibrium(init_EIR = eirs[i])
  
  ms_raw <- run_simulation((10 + 10) * year, parameters = p)
  
  ms <- ms_raw %>%
    mutate(prev_2_10 = n_detect_730_3650  / n_730_3650,
           clin = n_inc_clinical_0_36500 / n_0_36500,
           clin_smooth = roll_mean(clin, n = 365, align = "right", fill = NA) * 365,
           sev = n_inc_severe_0_36500 / n_0_36500,
           sev_smooth = roll_mean(sev, n = 365, align = "right", fill = NA, na.rm = TRUE) * 365,
           model = "malariasimulation") %>%
    filter(timestep > 10 * year) %>%
    mutate(timestep = timestep - 10 * year) %>%
    select(model, timestep, prev_2_10, clin, clin_smooth, sev, sev_smooth)
  
  # MalariaLaunchR run
  llin_options <- mlgts::itn_flexible_input(years = 0:9,
                                            coverage = c(rep(0, 5), rep(0.5, 5)),
                                            num_runs = 1)
  vector_options <- "itn_kill_fun 0.533 itn_repel_fun 0.56 itn_repel_min_fun 0.24 itn_kill_arab 0.533 itn_repel_arab 0.56 itn_repel_min_arab 0.24 itn_kill_gamb_ss 0.533 itn_repel_gamb_ss 0.56 itn_repel_min_gamb_ss 0.24"
  mlgts_raw <- mlgts::launch(name = "t1",
                             options = paste("num_people", pop, "prev", mean(ms$prev_2_10[1:1000]), "final_run", 10, "output_per_yr", 365, "recalculate", 2,
                                             "add", llin_options, vector_options),
                             demog = demog,
                             site = mlgts::perennial)
  mlgts <- mlgts_raw$output %>%
    rename(clin = clin_inc_all,
           clin_smooth = clin_inc_all_smooth,
           sev = sev_inc_all,
           sev_smooth = sev_inc_all_smooth) %>%
    mutate(
      clin_smooth = clin_smooth,
      sev_smooth = sev_smooth,
      timestep = year * 365,
      model = "old") %>%
    select(model, timestep, prev_2_10, clin, clin_smooth, sev, sev_smooth)
  
  
  # Combine output
  output[[i]] <- bind_rows(ms, mlgts) %>%
    mutate(eir = eirs[i])
}
output <- bind_rows(output)

prev_plot <- ggplot(output, aes(x = timestep, y = prev_2_10, col = model)) +
  geom_line() +
  ylim(0, 0.7) +
  theme_bw() +
  facet_wrap( ~ eir, nrow = 1)
inc_plot <- ggplot(output, aes(x = timestep, y = clin_smooth, col = model)) +
  geom_line() +
  ylim(0, 1) +
  ylab("Smooth annual\nclinical incidence") + 
  theme_bw() +
  facet_wrap( ~ eir, nrow = 1)
sev_plot <- ggplot(output, aes(x = timestep, y = sev_smooth, col = model)) +
  geom_line() +
  ylim(0, 0.025) +
  ylab("Smooth annual\nsevere incidence") + 
  theme_bw() +
  facet_wrap( ~ eir, nrow = 1)

summary_plot <- (prev_plot / inc_plot / sev_plot) + plot_layout(guide = "collect")

ggsave("plots/comparison_noSeasonality_nets_exponentialDemog.png", summary_plot, width = 10, height =8)
