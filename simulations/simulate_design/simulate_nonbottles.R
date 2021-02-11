# Simulate design 1 (all trials come from same source population)
# SN February 11, 2021

### Setup

rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)

source('base_model_source/sim_functions.R')

pars = expand.grid(
  alpha  = c(0, 0.0035),
  evolve = c(TRUE, FALSE),
  n.pop0 = c(20, 100),
  block  = 1:4,
  rep    = 1:6
) %>%
  mutate(
    init.row = 2e4,
    end.time = 15,
    n.loci = 25,
    w.max = 2,
    theta = 2.775,
    wfitn = sqrt(1 / 0.14 / 2),
    sig.e = sqrt(0.5),
    pos.p = 0.5,
    trial = 1:nrow(.)
)

# Design table
design = distinct(pars, trial, n.pop0, alpha, evolve)

# Define number of trials
n.trials = nrow(pars)

# Initialize object for storing sims
liszt = vector('list', n.trials)

### Run simulations

# Set seed for reproducibility
set.seed(5280)

# Run sims
for (trial in 1:n.trials) {
  liszt[[trial]] = sim(
    params = pars[trial,],
    evolve = pars$evolve[trial]
  )
  print(trial)
}

# Convert simulations to long data frame
experiment = unroller(liszt)
nrow(experiment)

# Summarise (means for generation in each time step)
experiment.summ = experiment %>%
  group_by(trial, gen) %>%
  summarise(n = n(),
            g = mean(g_i),
            z = mean(z_i),
            w = mean(w_i)) %>%
  merge(y = design, by = 'trial')

# Write csv
write.csv(experiment.summ %>% mutate(design = 0), 
          row.names = FALSE,
          'simulations/simulate_design/output/design_nb_summary.csv')

