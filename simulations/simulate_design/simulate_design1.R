# Simulate design 1 (all trials come from same source population)
# SN February 10, 2021

### Setup

rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)

source('base_model_source/sim_functions.R')

pars = data.frame(
  init.row = 2e4,
  end.time = 15,
  n.loci = 25,
  w.max = 2,
  theta = 2.775,
  wfitn = sqrt(1 / 0.14 / 2),
  sig.e = sqrt(0.5),
  pos.p = 0.5
)

# Read in source populations
bottles = read.csv('simulations/prepare_simulations/outputs/bottlenecks_design1.csv') %>%
  # For some reason, export scrambled rows. Think this caused problems. Fix here:
  arrange(sampno, i) %>%
  # I think including the 'trial' column here is causing problems.
  select(-trial)

head(bottles)

# Design table
design = distinct(bottles, sampno, n.pop0, alpha, evolve)

# Grouping variable: sampno
table(bottles$sampno)
# (note: all trials have the same source population)

# Number of simultions to run
n.trials = length(unique(bottles$sampno))

# Initialize object for storing sims
liszt = vector('list', n.trials)

### Run simulations

# Set seed for reproducibility
set.seed(992)

# Run sims
for (trial in 1:n.trials) {
  liszt[[trial]] = sim(
    params = pars %>% 
      mutate(n.pop0 = design$n.pop0[trial],
             alpha  = design$alpha[trial]), 
    init.popn = bottles %>% filter(sampno %in% trial),
    evolve = design$evolve[trial]
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
  merge(y = design, by.x = 'trial', by.y = 'sampno')

# Write csv
write.csv(experiment.summ %>% mutate(design = 1), 
          row.names = FALSE,
          'simulations/simulate_design/output/design1_summary.csv')
