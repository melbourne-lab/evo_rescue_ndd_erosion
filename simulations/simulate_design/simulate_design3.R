# Simulate design 3 (one source population per rep - almost)
# (saying almost because pop size is paired, i.e.,
# for each combintion of block x rep x alpha x evolve,
# the two populations are from the same source)
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
  n.pop0 = 50,
  n.loci = 25,
  w.max = 2,
  theta = 2.775,
  wfitn = sqrt(1 / 0.14 / 2),
  sig.e = sqrt(0.5),
  pos.p = 0.5,
  alpha = 0
)

# Read in source populations
bottles = read.csv('simulations/prepare_simulations/outputs/bottlenecks_design3.csv') %>%
  # Rename 'trial' population to 'sourcepop'
  rename(sourcepop = trial)

head(bottles)

# Design table (linking trial/sample number to experiment attributes)
design = distinct(bottles, sampno, n.pop0, alpha, evolve, sourcepop, block, rep)

# Grouping variable: sampno
table(bottles$sampno)
# Note: sourcepop is source population.

# Number of simultions to run
n.trials = length(unique(bottles$sampno))

# Initialize object for holding simulation trials
liszt = vector('list', n.trials)

### Run simulations
# Set seed
set.seed(700)

# Run simulations
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

# Convert to long-form data frame
experiment = unroller(liszt)
nrow(experiment)

# Summarise for each trial, merge with design
experiment.summ = experiment %>%
  group_by(trial, gen) %>%
  summarise(n = n(),
            g = mean(g_i),
            z = mean(z_i),
            w = mean(w_i)) %>%
  merge(y = design, by.x = 'trial', by.y = 'sampno')

# Write csv
write.csv(experiment.summ %>% mutate(design = 3), 
          row.names = FALSE,
          'simulations/simulate_design/output/design3_summary.csv')
