# Script for simulating population dynamics of an adapting bottlenecked populations.
# This script runs more trials but saves only population size - used just to get
# population size (and therefore extinction probabilities) but doesn't save
# other information.
# SN - adapted from simulations/sim_thrice_bottled/simulate_evolving_a035.R
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

# Load source populations
# (do this now to get the population numbers)
neck.popns = read.csv('simulations/prepare_simulations/outputs/thrice_bottlenecked_popns.csv',
                      stringsAsFactors = FALSE)

# Define number of trials per source populations
trials.per = 250

pars = expand.grid(pop.id = unique(neck.popns$trial),
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
         alpha = 0.0035)

liszt.vlo = vector('list', nrow(pars))

set.seed(158709)

for (i in 1:nrow(pars)) {
  popid = pars$pop.id[i]
  n0 = pars$n.pop0[i]
  liszt.vlo[[i]] = sim( a = c(1/2, -1/2), 
                        params = pars[i,],
                        init.popn = neck.popns %>%
                          filter(trial %in% popid ) %>% # [neck.popns$trial %in% pars$pop.id[i],] %>%
                          sample_n(n0, replace = FALSE) %>%
                          select(-trial)
  ) %>%
    # Get population size in each step
    group_by(gen) %>%
    summarise(n = n()) %>%
    # Add zeros if needed
    rbind(data.frame(gen = 1:pars$end.time[1], n = 0)) %>%
    group_by(gen) %>%
    summarise(n = sum(n)) %>%
    ungroup()
  print(paste0('ndd bottle ', i, ' of ', nrow(pars)))
}

all.vlo.sizes = unroller(liszt.vlo)

head(all.vlo.trials)

### Normie populations (i.e., non-bottlenecked)

liszt.vhi = vector('list', nrow(pars))

set.seed(940312)

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
  print(paste0('ndd non-bottle ', i, ' of ', nrow(pars)))
}

all.vhi.sizes = unroller(liszt.vhi)

### Combine the results.

all.n.summ = rbind(
  all.vlo.sizes %>% 
    mutate(bottleneck = TRUE) %>%
    merge(y = pars, by.x = 'trial', by.y = 'global.trial'),
  all.vhi.sizes %>% 
    mutate(bottleneck = FALSE) %>%
    merge(y = pars, by.x = 'trial', by.y = 'global.trial') %>%
    mutate(trial = trial + max(all.vlo.sizes$trial))
) %>%
  group_by(trial, bottleneck) %>% 
  mutate(extinct = any(!n)) %>%
  group_by(n.pop0, gen, bottleneck) %>% 
  summarise(nbar = mean(n), 
            nvar = var(n),
            pext = mean(n>0)) %>%
  ungroup()

write.csv(all.n.summ, na = '', row.names = FALSE,
          file = 'simulations/outputs/thrice_bottlenecked_results/a035_evotrue_all_n_bigsamps.csv')
