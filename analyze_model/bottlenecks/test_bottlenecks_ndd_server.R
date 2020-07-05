# Script for simulating population dynamics of an adapting bottlenecked
# populations (alpha > 0)
# SN - adopted from the old bottleneck server script
# "written" and run on server june 30 2020

### Clear namespace
rm(list = ls())

### Load packages

library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

### Set working directory
setwd('~/demall/')

### Load source materials

# Get simulation functions
source('base_model_source/sim_functions.R')

# Load source populations
# (do this now to get the population numbers)
neck.popns = read.csv('bottlenecks/out/all_k_pops_min_gen.csv',
                      stringsAsFactors = FALSE)

# Define number of trials per source populations
trials.per = 100

pars = expand.grid(pop.id = unique(neck.popns$trial),
                   pop.trial = 1:trials.per,
                   n.pop0 = c(20, 100)) %>%
  mutate(global.trial = 1:nrow(.),
         end.time = 15,
         init.row = 1e4,
         n.loci = 25,
         w.max = 2,
         theta = 3,
         wfitn = sqrt(1 / 0.14 / 2),
         sig.e = 0,
         pos.p = 0.5,
         alpha = 0.0035)

liszt.vlo = vector('list', nrow(pars))

set.seed(7756758)

for (i in 1:nrow(pars)) {
  popid = pars$pop.id[i]
  n0 = pars$n.pop0[i]
  liszt.vlo[[i]] = sim( a = c(1/2, -1/2), 
                        params = pars[i,],
                        init.popn = neck.popns %>%
                          filter(trial %in% popid ) %>% # [neck.popns$trial %in% pars$pop.id[i],] %>%
                          sample_n(n0, replace = FALSE) %>%
                          select(-trial)
  )
  print(i)
}

all.vlo.trials = unroller(liszt.vlo)

head(all.vlo.trials)

all.vlo.trials %>%
  group_by(trial, gen) %>%
  summarise(n = n())

all.vlo.n = all.vlo.trials %>%
  group_by(trial, gen) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  rbind(expand.grid(trial = 1:nrow(pars), gen = 1:pars$end.time[1], n = 0)) %>%
  group_by(trial, gen) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  merge(y = pars, by.x = 'trial', by.y = 'global.trial')

### normie populations

liszt.vhi = vector('list', nrow(pars))

set.seed(69594)

for (i in 1:nrow(pars)) {
  liszt.vhi[[i]] = sim( params = pars[i,])
  print(i)
}

all.vhi.trials = unroller(liszt.vhi)

head(all.vlo.trials)

all.vhi.trials %>%
  group_by(trial, gen) %>%
  summarise(n = n())

all.vhi.n = all.vhi.trials %>%
  group_by(trial, gen) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  rbind(expand.grid(trial = 1:nrow(pars), gen = 1:pars$end.time[1], n = 0)) %>%
  group_by(trial, gen) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  merge(y = pars, by.x = 'trial', by.y = 'global.trial')

### Combine the results.

all.n = rbind(
  all.vlo.n %>% mutate(bottleneck = TRUE),
  all.vhi.n %>% mutate(bottleneck = FALSE, trial = trial + max(all.vlo.n$trial))
) %>%
  group_by(trial, bottleneck) %>% 
  mutate(extinct = any(!n)) %>%
  group_by(n.pop0, gen, bottleneck) %>% 
  mutate(nbar = mean(n), 
         nvar = var(n),
         pext = mean(n>0)) %>%
  ungroup() %>%
  filter(n > 1)

write.csv(all.n, na = '', row.names = FALSE,
          file = 'out/bottleneck_ndd_test_pop_sizes.csv')

# Save all simulation output
save(liszt.vhi, liszt.vlo, file = 'out/bottleneck_ndd_adapted_outputs.RData')

