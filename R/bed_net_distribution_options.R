### Comparing bed net distribution options #####################################


# Bednets are commonly distributed via mass campaigns on multi (usual 3-year) cycles (possibly with some top up on off-years).
# A common use-case of the "old" model was to smooth out the impact of these distribution cycles.
# For a 3-yearly cycle this was achieved by splitting the population into 3 and distributing coverage % to each third
# of the population over the 3 years of the cycle. 
# Here I explore methods of capturing this kind of smoothed output using malaria simulation, for a fixed eir.

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
eir <- 5

p <- get_parameters(list(
  human_population = pop,
  prevalence_rendering_min_ages = 2 * year,
  prevalence_rendering_max_ages = 10 * year,
  clinical_incidence_rendering_min_ages = 0 * year,
  clinical_incidence_rendering_max_ages = 100 * year,
  severe_incidence_rendering_min_ages = 0 * year,
  severe_incidence_rendering_max_ages = 100 * year,
  individual_mosquitoes = FALSE
))


### Malaria simulation - 3-yearly cycle ########################################

# 50% distribution every 3 years - modelling a realistic cycle
p_3year <- p %>%
  set_bednets(
    timesteps = as.integer(c(0, 0.10 + (13:20)) * year),
    coverages = c(0, c(0.5, 0, 0, 0.5, 0, 0, 0.5, 0)),
    retention = 5 * year,
    dn0 = matrix(c(0.533), nrow = 9, ncol = 1),
    rn = matrix(c(0.56), nrow = 9, ncol = 1),
    rnm = matrix(c(0.24), nrow = 9, ncol = 1),
    gamman = rep(2.64 * 365, 9)) %>%
  set_equilibrium(init_EIR = eir)

ms_raw_3year <- run_simulation((10 + 10) * year, parameters = p_3year)

ms_3year <- ms_raw_3year %>%
  mutate(prev_2_10 = n_detect_730_3650  / n_730_3650,
         prev_2_10_smooth = roll_mean(prev_2_10, n = 365 * 1, align = "right", fill = NA),
         clin = n_inc_clinical_0_36500 / n_0_36500,
         clin_smooth = roll_mean(clin, n = 365 * 1, align = "right", fill = NA) * 365,
         sev = n_inc_severe_0_36500 / n_0_36500,
         sev_smooth = roll_mean(sev, n = 365 * 1, align = "right", fill = NA, na.rm = TRUE) * 365,
         model = "malariasimulation",
         method = "3 year distribution, 1 year smooth") %>%
  filter(timestep > 10 * year) %>%
  mutate(timestep = timestep - 10 * year) %>%
  select(model, method, timestep, prev_2_10, prev_2_10_smooth, clin, clin_smooth, sev, sev_smooth)

# A second option is to smooth (the same output) over 3 years
ms_3year_3year_smooth <- ms_raw_3year %>%
  mutate(prev_2_10 = n_detect_730_3650  / n_730_3650,
         prev_2_10_smooth = roll_mean(prev_2_10, n = 365 * 3, align = "right", fill = NA),
         clin = n_inc_clinical_0_36500 / n_0_36500,
         clin_smooth = roll_mean(clin, n = 365 * 3, align = "right", fill = NA) * 365,
         sev = n_inc_severe_0_36500 / n_0_36500,
         sev_smooth = roll_mean(sev, n = 365 * 3, align = "right", fill = NA, na.rm = TRUE) * 365,
         model = "malariasimulation",
         method = "3 year distribution, 3 year smooth") %>%
  filter(timestep > 10 * year) %>%
  mutate(timestep = timestep - 10 * year) %>%
  select(model, method, timestep, prev_2_10, prev_2_10_smooth, clin, clin_smooth, sev, sev_smooth)

# 50% / 3 distribution every year (uncorrelated) - modelling a realistic cycle

third_cov <- 0.5 / 3
p_1year <- p %>%
  set_bednets(
    timesteps = as.integer(c(0, 0.10 + (13:20)) * year),
    coverages = c(0, rep(third_cov, 8)),
    retention = 5 * year,
    dn0 = matrix(c(0.533), nrow = 9, ncol = 1),
    rn = matrix(c(0.56), nrow = 9, ncol = 1),
    rnm = matrix(c(0.24), nrow = 9, ncol = 1),
    gamman = rep(2.64 * 365, 9)) %>%
  set_equilibrium(init_EIR = eir)

# Set uncorrelated rounds
cors_1year <- malariasimulation::get_correlation_parameters(p_1year)
cors_1year$inter_round_rho("bednets", 0)

ms_raw_1year <- run_simulation((10 + 10) * year, parameters = p_1year, correlations = cors_1year)

ms_1year <- ms_raw_1year %>%
  mutate(prev_2_10 = n_detect_730_3650  / n_730_3650,
         prev_2_10_smooth = roll_mean(prev_2_10, n = 365 * 1, align = "right", fill = NA),
         clin = n_inc_clinical_0_36500 / n_0_36500,
         clin_smooth = roll_mean(clin, n = 365 * 1, align = "right", fill = NA) * 365,
         sev = n_inc_severe_0_36500 / n_0_36500,
         sev_smooth = roll_mean(sev, n = 365 * 1, align = "right", fill = NA, na.rm = TRUE) * 365,
         model = "malariasimulation",
         method = "Annual distribution (at 1/3 of 3 year coverage)") %>%
  filter(timestep > 10 * year) %>%
  mutate(timestep = timestep - 10 * year) %>%
  select(model, method, timestep, prev_2_10, prev_2_10_smooth, clin, clin_smooth, sev, sev_smooth)
################################################################################


### MalariaLaunchR run #########################################################
llin_options <- mlgts::itn_flexible_input(years = 0:9,
                                          coverage = c(rep(0, 3), rep(0.5, 7)),
                                          num_runs = 3)
vector_options <- "itn_kill_fun 0.533 itn_repel_fun 0.56 itn_repel_min_fun 0.24 itn_kill_arab 0.533 itn_repel_arab 0.56 itn_repel_min_arab 0.24 itn_kill_gamb_ss 0.533 itn_repel_gamb_ss 0.56 itn_repel_min_gamb_ss 0.24"
mlgts_raw <- mlgts::launch(name = "t1",
                           options = paste("num_people", pop, "prev", mean(ms_3year$prev_2_10[1:1000]), "final_run", 10, "output_per_yr", 365, "recalculate", 2,
                                           "num_runs 3 add", llin_options, vector_options),
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
    model = "old",
    method = "Annual distirbution, population split") %>%
  select(model, method, timestep, prev_2_10, prev_2_10_smooth, clin, clin_smooth, sev, sev_smooth)
################################################################################

output <- bind_rows(ms_3year, ms_3year_3year_smooth, ms_1year, mlgts)

prev_plot <- ggplot(output, aes(x = timestep, y = prev_2_10_smooth, col = method, lty = model)) +
  geom_line() +
  ylim(0, 0.7) +
  ylab("Smooth prevalence") + 
  theme_bw() 
inc_plot <- ggplot(output, aes(x = timestep, y = clin_smooth, col = method, lty = model)) +
  geom_line() +
  ylim(0, 1) +
  ylab("Smooth annual\nclinical incidence") + 
  theme_bw() 
sev_plot <- ggplot(output, aes(x = timestep, y = sev_smooth, col = method, lty = model)) +
  geom_line() +
  ylim(0, 0.025) +
  ylab("Smooth annual\nsevere incidence") + 
  theme_bw()

summary_plot <- (prev_plot / inc_plot / sev_plot) + plot_layout(guide = "collect")

ggsave("plots/comparison_net_distribution.png", summary_plot, width = 10, height =8)
