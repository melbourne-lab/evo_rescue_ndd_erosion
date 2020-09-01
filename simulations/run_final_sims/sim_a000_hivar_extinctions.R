# Script for simulating population dynamics of an adapting populations.
# Here, I am extracting mean genotype, phenotype, and fitness for all eight combos.
# In this script: low genetic variation, density independent.
# SN - adapted from simulations/sim_thrce_bottlenecked/sim_a000_hivar_extinction.R
# adapted and run on August 31 2020

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

liszt.vhi = vector('list', nrow(pars))

set.seed(19935)

for (i in 1:nrow(pars)) {
  liszt.vhi[[i]] = sim( a = c(1/2, -1/2), params = pars[i,]) %>%
    # Get population size in each step
    group_by(gen) %>%
    summarise(n = n()) %>%
    # Add zeros if needed
    rbind(data.frame(gen = 1:pars$end.time[1], n = 0)) %>%
    group_by(gen) %>%
    summarise(n = sum(n)) %>%
    ungroup()
  print(paste0('hi var ', i, ' of ', nrow(pars)))
}

all.vhi.sizes = unroller(liszt.vhi)

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
          file = 'simulations/outputs/final_results/a000_hivar_extinctions.csv')
