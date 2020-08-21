# Script for simulating population dynamics of an adapting bottlenecked populations.
# This is scraping out genetic varation, but trying to avoid the
# censoring/biasing occurring due to extinction.
# I realized after running on a small number (10-per) of trials, though, that
# I'm most interested in genotypes, not genetic variance. Lmao. This script has
# not been run in full.
# Still would be interesting to run this! It would show genetic variation loss
# in *large* populations, i.e., with less genetic stochasticity.
# If re-running, make sure you change the number of iterations in the loop my guy!
# SN - August 21 2020 

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
                      stringsAsFactors = FALSE) %>%
  group_by(trial) %>%
  filter(n() > 120) %>%
  ungroup()

# Define number of trials per source populations
trials.per = 25

pars = expand.grid(pop.id = unique(neck.popns$trial),
                   pop.trial = 1:trials.per) %>%
  mutate(global.trial = 1:nrow(.),
         end.time = 10,
         init.row = 1e4,
         n.loci = 25,
         n.pop0 = 120,
         w.max = 2,
         theta = 2.75,
         wfitn = sqrt(1 / 0.14 / 2),
         sig.e = 0.5,
         pos.p = 0.5,
         alpha = 0.0035)

liszt.vlo = vector('list', nrow(pars))

set.seed(62180)

for (i in 1:nrow(pars)) {
  popid = pars$pop.id[i]
  n0 = pars$n.pop0[i]
  # Run simulation then extract the genetic summaries.
  liszt.vlo[[i]] = sim( a = c(1/2, -1/2), 
                        params = pars[i,],
                        init.popn = neck.popns %>%
                          filter(trial %in% popid ) %>% # [neck.popns$trial %in% pars$pop.id[i],] %>%
                          sample_n(n0, replace = FALSE) %>%
                          select(-trial)
  ) %>%
    # Summarise results.
    group_by(gen) %>%
    summarise(g_t = mean(g_i),
              z_t = mean(z_i),
              w_t = mean(w_i))
  print(paste0('ndd big N0 bottle ', i, ' of ', nrow(pars)))
}

all.vlo.trials = unroller(liszt.vlo)

### Normie populations (i.e., non-bottlenecked)

liszt.vhi = vector('list', nrow(pars))

set.seed(249290)

for (i in 1:nrow(pars)) {
  liszt.vhi[[i]] = sim( params = pars[i,]) %>%
    # Summarise results.
    group_by(gen) %>%
    summarise(g_t = mean(g_i),
              z_t = mean(z_i),
              w_t = mean(w_i))
  print(paste0('ndd big N0 non-bottle ', i, ' of ', nrow(pars)))
}

all.vhi.trials = unroller(liszt.vhi)

### Combine the results.

# Commenting out population size and demo, which I don't think we need (right now).

# all.n.summ = rbind(
#   all.vlo.n %>% mutate(bottleneck = TRUE),
#   all.vhi.n %>% mutate(bottleneck = FALSE, trial = trial + max(all.vlo.n$trial))
# ) %>%
#   group_by(trial, bottleneck) %>% 
#   mutate(extinct = any(!n)) %>%
#   group_by(n.pop0, gen, bottleneck) %>% 
#   summarise(nbar = mean(n), 
#             nvar = var(n),
#             pext = mean(n>0)) %>%
#   ungroup()
# 
# write.csv(all.n.summ, na = '', row.names = FALSE,
#           file = 'simulations/outputs/thrice_bottlenecked_results/a000_evotrue_all_n.csv')
# 
all.demo.summs = rbind(
  all.vlo.trials %>%
    mutate(bottleneck = TRUE) %>%
    merge(y = pars %>% select(global.trial, n.pop0),
          by.x = 'trial', by.y = 'global.trial') %>%
    mutate(bottleneck = TRUE),
  all.vhi.trials %>%
    merge(y = pars %>% select(global.trial, n.pop0),
          by.x = 'trial', by.y = 'global.trial') %>%
    mutate(bottleneck = FALSE, trial = trial + max(all.vlo.trials$trial))
) %>%
  group_by(gen, bottleneck, n.pop0) %>%
  summarise(gbar = mean(g_t),
            gvar = var(g_t),
            zbar = mean(z_t),
            zvar = var(z_t),
            wbar = mean(w_t),
            wvar = var(w_t),
            n.trials = length(unique(trial)))

write.csv(all.demo.summs, na = '', row.names = FALSE,
          file = 'simulations/outputs/thrice_bottlenecked_results/a035_gt_big_init.csv')

# all.gene.summs = rbind(
#   all.vlo.trials %>% 
#     mutate(bottleneck = TRUE) %>%
#     merge(y = pars %>% select(global.trial),
#           by.x = 'trial', by.y = 'global.trial') %>%
#     mutate(bottleneck = TRUE),
#   all.vhi.trials %>%
#     merge(y = pars %>% select(global.trial),
#           by.x = 'trial', by.y = 'global.trial') %>%
#     mutate(bottleneck = FALSE, trial = trial + nrow(pars))
# ) %>%
#   # select(-c(i, g_i, z_i, w_i, r_i, fem)) %>%
#   # gather(key = loc.copy, value = val, -c(gen, trial, bottleneck, n.pop0)) %>%
#   # mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
#   # group_by(trial, bottleneck, n.pop0, gen, locus) %>%
#   # summarise(p = mean(val > 0)) %>%
#   # group_by(trial, bottleneck, n.pop0, gen) %>%
#   # summarise(p.fix.pos = mean(p == 1),
#   #           p.fix.neg = mean(p == 0),
#   #           v = sum(2 * p * (1 - p)) / pars$n.loci[1]) %>%
#   group_by(bottleneck, gen) %>%
#   summarise(p.fix.pos = mean(p.fix.pos),
#             p.fix.neg = mean(p.fix.neg),
#             vbar = mean(v),
#             vvar = var(v),
#             ntrials = n())
# 
# write.csv(all.gene.summs, na = '', row.names = FALSE,
#           file = 'simulations/outputs/thrice_bottlenecked_results/a000_evotrue_big_init_gene.csv')

# Save all simulation output
# save(liszt.vhi, liszt.vlo, file = 'out/bottleneck_adapted_outputs.RData')

