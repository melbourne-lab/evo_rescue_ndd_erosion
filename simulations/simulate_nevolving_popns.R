# Script for simulating population dynamics of an adapting bottlenecked populations.
# SN - adapted from analyze_model/bottlenecks/test_bottlenecks_server
# adapted and run on July 6 2020

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
neck.popns = read.csv('simulations/prepare_simulations/outputs/all_k_pops_min_gen.csv',
                      stringsAsFactors = FALSE)

# Define number of trials per source populations
trials.per = 26

pars = expand.grid(pop.id = unique(neck.popns$trial),
                   pop.trial = 1:trials.per,
                   n.pop0 = c(20, 100)) %>%
  mutate(global.trial = 1:nrow(.),
         end.time = 15,
         init.row = 1e4,
         n.loci = 25,
         w.max = 2,
         theta = 2.6,
         wfitn = sqrt(1 / 0.14 / 2),
         sig.e = 0.5,
         pos.p = 0.5,
         alpha = 0)

liszt.vlo = vector('list', nrow(pars))

set.seed(450332)

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
  print(paste0('non-evolving bottle ', i, ' of ', nrow(pars)))
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

### Normie populations (i.e., non-bottlenecked)

liszt.vhi = vector('list', nrow(pars))

set.seed(48686)

for (i in 1:nrow(pars)) {
  liszt.vhi[[i]] = sim( params = pars[i,])
  print(paste0('non-evolving non-bottle ', i, ' of ', nrow(pars)))
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

all.n.summ = rbind(
  all.vlo.n %>% mutate(bottleneck = TRUE),
  all.vhi.n %>% mutate(bottleneck = FALSE, trial = trial + max(all.vlo.n$trial))
) %>%
  group_by(trial, bottleneck) %>% 
  mutate(extinct = any(!n)) %>%
  group_by(n.pop0, gen, bottleneck) %>% 
  summarise(nbar = mean(n), 
            nvar = var(n),
            pext = mean(n>0)) %>%
  ungroup()

write.csv(all.n.summ, na = '', row.names = FALSE,
          file = 'simulations/outputs/a000_evotrue_all_n.csv')

all.demo.summs = rbind(
  all.vlo.trials %>% 
    mutate(bottleneck = TRUE) %>%
    merge(y = pars %>% select(global.trial, n.pop0),
          by.x = 'trial', by.y = 'global.trial') %>%
    mutate(bottleneck = TRUE),
  all.vhi.trials %>%
    merge(y = pars %>% select(global.trial, n.pop0),
          by.x = 'trial', by.y = 'global.trial') %>%
    mutate(bottleneck = FALSE, trial = trial + max(all.vlo.n$trial))
) %>%
  group_by(gen, bottleneck, n.pop0) %>%
  summarise(gbar = mean(g_i),
            gvar = var(g_i),
            zbar = mean(z_i),
            zvar = var(z_i),
            wbar = mean(w_i),
            wvar = var(w_i),
            n.trials = length(unique(trial)))

write.csv(all.demo.summs, na = '', row.names = FALSE,
          file = 'simulations/outputs/a000_evotrue_all_demo.csv')

all.gene.summs = rbind(
  all.vlo.trials %>% 
    mutate(bottleneck = TRUE) %>%
    merge(y = pars %>% select(global.trial, n.pop0),
          by.x = 'trial', by.y = 'global.trial') %>%
    mutate(bottleneck = TRUE),
  all.vhi.trials %>%
    merge(y = pars %>% select(global.trial, n.pop0),
          by.x = 'trial', by.y = 'global.trial') %>%
    mutate(bottleneck = FALSE, trial = trial + max(all.vlo.n$trial))
) %>%
  select(-c(i, g_i, z_i, w_i, r_i, fem)) %>%
  gather(key = loc.copy, value = val, -c(gen, trial, bottleneck, n.pop0)) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(trial, bottleneck, n.pop0, gen, locus) %>%
  summarise(p = mean(val > 0)) %>%
  group_by(trial, bottleneck, n.pop0, gen) %>%
  summarise(p.fix.pos = mean(p == 1),
            p.fix.neg = mean(p == 0),
            v = sum(2 * p * (1 - p)) / pars$n.loci[1]) %>%
  group_by(bottleneck, n.pop0, gen) %>%
  summarise(p.fix.pos = mean(p.fix.pos),
            p.fix.neg = mean(p.fix.neg),
            vbar = mean(v),
            vvar = var(v))

write.csv(all.gene.summs, na = '', row.names = FALSE,
          file = 'simulations/outputs/a000_evotrue_all_gene.csv')

# Save all simulation output
# save(liszt.vhi, liszt.vlo, file = 'out/bottleneck_adapted_outputs.RData')

