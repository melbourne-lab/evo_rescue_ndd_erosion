# Script for simulating population dynamics of an adapting populations.
# Here, I am extracting mean population size and extinctions.
# In this script: low genetic variation, density independent.
# SN - adapted from simulations/run_final_sims/sim_a000_hivar_extinctions.R
# adapted and run on September 1, 2020

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


# Define number of trials per source populations
trials = 8000

pars = expand.grid(local.trial = 1:trials,
                   n.pop0 = c(20, 60)) %>%
  mutate(global.trial = 1:nrow(.),
         end.time = 15,
         init.row = 1e4,
         n.loci = 25,
         w.max = 2,
         theta = 2.75,
         wfitn = sqrt(1 / 0.14 / 2),
         sig.e = sqrt(0.5),
         pos.p = 0.5,
         alpha = 0)

liszt.vlo = vector('list', nrow(pars))

set.seed(132190)

for (i in 1:nrow(pars)) {
  pop.init = init.sim(a = c(1/2, -1/2),
                      params = pars[i,]) %>%
    mutate_at(paste0('a', 1:6), function(x) -1/2) %>%
    mutate_at(paste0('b', 1:6), function(x) -1/2) %>%
    mutate_at(paste0('a', 7:12), function(x) 1/2) %>%
    mutate_at(paste0('b', 7:12), function(x) 1/2)
  liszt.vlo[[i]] = sim( a = c(1/2, -1/2), 
                        params = pars[i,],
                        init.popn = pop.init
  ) %>%
    # Get population size in each step
    group_by(gen) %>%
    summarise(n = n()) %>%
    # Add zeros if needed
    rbind(data.frame(gen = 1:pars$end.time[1], n = 0)) %>%
    group_by(gen) %>%
    summarise(n = sum(n)) %>%
    ungroup()
  print(paste0('lo var ', i, ' of ', nrow(pars)))
}

all.vlo.sizes = unroller(liszt.vlo)

all.n.summ = all.vlo.sizes %>% 
  # Add paramters, etc.
  mutate(bottleneck = TRUE) %>%
  merge(y = pars, by.x = 'trial', by.y = 'global.trial') %>%
  # Summarise
  group_by(trial) %>% 
  group_by(n.pop0, gen) %>% 
  summarise(nbar = mean(n), 
            nvar = var(n),
            pext = mean(n>0)) %>%
  ungroup()

write.csv(all.n.summ, na = '', row.names = FALSE,
          file = 'simulations/outputs/final_results/a000_lovar_extinctions.csv')
