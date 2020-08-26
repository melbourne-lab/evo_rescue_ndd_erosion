# Script for simulating population dynamics of an adapting bottlenecked populations.
# This script runs more trials but saves only population size - used just to get
# population size (and therefore extinction probabilities) but doesn't save
# other information.
# SN - adapted from simulations/sim_thrice_bottled/simulate_evolving_a000.R
# run on August 21 2020

### Clear namespace
rm(list = ls())

### Load packages

library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

### Load source materials

# Get simulation functions
source('base_model_source/sim_functions.R')

# Get number of bottlenecked source populations
# (ensures that there's the same number of trials to compare with the
# bottlenecked sims)
n.neck.popns = read.csv('simulations/prepare_simulations/outputs/thrice_bottlenecked_popns.csv',
                      stringsAsFactors = FALSE) %>%
  distinct(trial) %>%
  nrow()

# Define number of trials per source populations
trials.per = 250

pars = expand.grid(pop.id = 1:n.neck.popns,
                   pop.trial = 1:trials.per,
                   n.pop0 = c(20, 60)) %>%
  mutate(global.trial = 1:nrow(.),
         end.time = 15,
         init.row = 1e4,
         n.loci = 25,
         w.max = 2,
         theta = 2.75,
         wfitn = sqrt(1 / 0.14 / 2),
         sig.e = 0.5,
         pos.p = 0.5,
         alpha = 0)

### Normie populations (i.e., non-bottlenecked)

liszt.vhi = vector('list', nrow(pars))

set.seed(979388)

for (i in 1:nrow(pars)) {
  liszt.vhi[[i]] = sim( params = pars[i,]) %>%
    # Get population size in each step
    group_by(gen) %>%
    summarise(n = n()) %>%
    # Add zeros if needed
    rbind(data.frame(gen = 1:pars$end.time[1], n = 0)) %>%
    group_by(gen) %>%
    summarise(n = sum(n)) %>%
    ungroup()
  print(paste0('non-bottle ', i, ' of ', nrow(pars)))
}

all.vhi.sizes = unroller(liszt.vhi)

### Combine the results.

all.n.summ = all.vhi.sizes %>% 
  # Add paramters, etc.
  mutate(bottleneck = FALSE) %>%
  merge(y = pars, by.x = 'trial', by.y = 'global.trial') %>%
  # Summarise
  group_by(trial) %>% 
  group_by(n.pop0, gen) %>% 
  summarise(nbar = mean(n), 
            nvar = var(n),
            pext = mean(n>0)) %>%
  ungroup()

write.csv(all.n.summ, na = '', row.names = FALSE,
          file = 'simulations/outputs/thrice_bottlenecked_results/a000_hivar_extinctions.csv')
