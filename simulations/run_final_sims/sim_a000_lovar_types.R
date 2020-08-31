# Script for simulating population dynamics of an adapting populations.
# Here, I am extracting mean genotype, phenotype, and fitness for all eight combos.
# In this script: low genetic variation, density independent.
# SN - adapted from simulations/sim_thrce_bottlenecked/sim_a000_lowvar_types.R
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
trials = 4000

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

set.seed(853443)

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
    group_by(gen) %>%
    summarise(n = n(),
              g_t = mean(g_i),
              z_t = mean(z_i),
              w_t = mean(w_i))
  print(paste0('lo var ', i, ' of ', nrow(pars)))
}

all.vlo.trials = unroller(liszt.vlo) %>%
  group_by(trial) %>%
  mutate(n.pop0 = n[1],
         ext.gen = max(gen),
         extinct = !ext.gen %in% pars$end.time[1]) %>%
  ungroup()

### Combine the results.

type.summ.all = all.vlo.trials %>%
  group_by(n.pop0, gen) %>%
  summarise(gbar = mean(g_t),
            gvar = var(g_t),
            zbar = mean(z_t),
            zvar = var(z_t),
            wbar = mean(w_t),
            wvar = var(w_t),
            n = n())

write.csv(type.summ.all, na = '', row.names = FALSE,
          file = 'simulations/outputs/final_results/a000_lovar_type_all.csv')

type.summ.ext = all.vlo.trials %>%
  group_by(n.pop0, extinct, gen) %>%
  summarise(gbar = mean(g_t),
            gvar = var(g_t),
            zbar = mean(z_t),
            zvar = var(z_t),
            wbar = mean(w_t),
            wvar = var(w_t),
            n = n())

write.csv(type.summ.ext, na = '', row.names = FALSE,
          file = 'simulations/outputs/final_results/a000_lovar_type_ext.csv')


# Save all simulation output
# save(liszt.vhi, liszt.vlo, file = 'out/bottleneck_adapted_outputs.RData')

type.summ.gen = all.vlo.trials %>%
  group_by(n.pop0, ext.gen, gen) %>%
  summarise(gbar = mean(g_t),
            gvar = var(g_t),
            zbar = mean(z_t),
            zvar = var(z_t),
            wbar = mean(w_t),
            wvar = var(w_t),
            n = n())

write.csv(type.summ.gen, na = '', row.names = FALSE,
          file = 'simulations/outputs/final_results/a000_lovar_type_gen.csv')
