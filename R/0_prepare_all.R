### Set up ###

source("R/functions.R")
library("malariasimulation")
library(dplyr)
library(ggplot2)
library(patchwork)
library(RcppRoll)
library(mlgts)

year <- 365
month <- 30
pop <- 20000
eirs <- c(0.1, 1, 5, 10, 50)
burnin <- 10
runtime <- 10

### Demography ###
exp_demog <- read.table("inputs/exponential_demog.txt", sep = "\t")
flat_demog <- mlgts::flat_demog


### Seasonality ###
no_season_site <- mlgts::perennial
no_season_site[4:10, 2] <- c(1, rep(0, 6))

seasonal_site <- mlgts::highly_seasonal

### Baseline parametrisation ###
p <- get_parameters(list(
  human_population = pop,
  prevalence_rendering_min_ages = 2 * year,
  prevalence_rendering_max_ages = 10 * year,
  clinical_incidence_rendering_min_ages = 0 * year,
  clinical_incidence_rendering_max_ages = 100 * year,
  severe_incidence_rendering_min_ages = 0 * year,
  severe_incidence_rendering_max_ages = 100 * year,
  individual_mosquitoes = FALSE,
  model_seasonality = TRUE,
  rainfall_floor = 0.01
)) %>%
  set_species(species = list(arab_params, fun_params, gamb_params), 
              proportions = c(0.25, 0.25, 0.5))