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
bottles = read.csv('experiment/expmt_bottleneck24/out/source_pops_design.csv') %>%
  # Add trial number (maybe should do in other file?)
  mutate(trialno = cumsum(i %in% 1))

head(bottles)

# Design table
design = bottles %>% distinct(trialno, sourcepops, n.pop0, alpha, evolve, block, repct)

nrow(design)

# Number of simulations to run
n.trials = nrow(design)

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
    init.popn = bottles %>% filter(trialno %in% trial),
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
  merge(y = design, by.x = 'trial', by.y = 'trialno')

experiment.summ
nrow(experiment.summ)

write.csv(experiment.summ, row.names = FALSE,
          file = 'experiment/expmt_bottleneck24/out/sim_results.csv')
