# Script for simulating population dynamics of an adapting bottlenecked populations.
# Here, I am extracting mean genotype, phenotype, and fitness for all eight combos.
# In this script: high genetic variation, density independent.
# SN - adapted from simulations/simulate_evolving_popns.R
# adapted and run on August 26 2020

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
neck.popns = read.csv('simulations/prepare_simulations/outputs/thrice_bottlenecked_600s.csv',
                      stringsAsFactors = FALSE)

# Define number of trials per source populations
trials.per = 75

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
         sig.e = sqrt(0.5),
         pos.p = 0.5,
         alpha = 0)

liszt.vhi = vector('list', nrow(pars))

set.seed(853443)

for (i in 1:nrow(pars)) {
  popid = pars$pop.id[i]
  n0 = pars$n.pop0[i]
  liszt.vhi[[i]] = sim( params = pars[i,]) %>%
    group_by(gen) %>%
    summarise(n = n(),
              g_t = mean(g_i),
              z_t = mean(z_i),
              w_t = mean(w_i))
  print(paste0('hi var ', i, ' of ', nrow(pars)))
}

all.vhi.trials = unroller(liszt.vhi) %>%
  group_by(trial) %>%
  mutate(n.pop0 = n[1],
         ext.gen = max(gen),
         extinct = !ext.gen %in% pars$end.time[1]) %>%
  ungroup()

### Combine the results.

type.summ.all = all.vhi.trials %>%
  group_by(n.pop0, gen) %>%
  summarise(gbar = mean(g_t),
            gvar = var(g_t),
            zbar = mean(z_t),
            zvar = var(z_t),
            wbar = mean(w_t),
            wvar = var(w_t),
            n = n())

write.csv(type.summ.all, na = '', row.names = FALSE,
          file = 'simulations/outputs/thrice_bottlenecked_results/a000_hivar_type_all.csv')

type.summ.ext = all.vhi.trials %>%
  group_by(n.pop0, extinct, gen) %>%
  summarise(gbar = mean(g_t),
            gvar = var(g_t),
            zbar = mean(z_t),
            zvar = var(z_t),
            wbar = mean(w_t),
            wvar = var(w_t),
            n = n())

write.csv(type.summ.ext, na = '', row.names = FALSE,
          file = 'simulations/outputs/thrice_bottlenecked_results/a000_hivar_type_ext.csv')


# Save all simulation output
# save(liszt.vhi, liszt.vhi, file = 'out/bottleneck_adapted_outputs.RData')

type.summ.gen = all.vhi.trials %>%
  group_by(n.pop0, ext.gen, gen) %>%
  summarise(gbar = mean(g_t),
            gvar = var(g_t),
            zbar = mean(z_t),
            zvar = var(z_t),
            wbar = mean(w_t),
            wvar = var(w_t),
            n = n())

write.csv(type.summ.gen, na = '', row.names = FALSE,
          file = 'simulations/outputs/thrice_bottlenecked_results/a000_hivar_type_gen.csv')
