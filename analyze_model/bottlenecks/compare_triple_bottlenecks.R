# Script for simulating population dynamics of an adapting bottlenecked populations.
# Unlike previous script, these populations have gone through three bottlenecks
# (further lowering their genetic variance)
# The expecttions here are
#   (1) populations experiencing the bottleneck will adapt more slowly
#   (i.e., the growth rate, and the growth rate increases, will  be lower
#   in the bottlenecked populations)
#   (2) (is there an influence of population size here?)
# I also use this script to look at how sigma_e (non-genetic phenotypic variance)
# influences adaptation and population dynamics
# (the expectation here is: it will also slow down adaptation - see Lande, 1976)
# SN - coded july 14 2020, run on melb2 on july 15 2020

### Clear namespace
rm(list = ls())

### Load packages

library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

### Load source materials

# Working directory (for server)
setwd('~/demgen_allelic')

# Get simulation functions
source('base_model_source/sim_functions.R')

# Load source populations
# (do this now to get the population numbers)
neck.popns = read.csv('analyze_model/bottlenecks/thrice_bottlenecked_popns.csv',
                      stringsAsFactors = FALSE)

# Define number of trials per source populations
trials.per = 12

pars = expand.grid(pop.id = unique(neck.popns$trial),
                   pop.trial = 1:trials.per,
                   n.pop0 = c(20, 50)) %>%
  mutate(global.trial = 1:nrow(.),
         end.time = 15,
         init.row = 1e4,
         n.loci = 25,
         w.max = 2,
         theta = 3,
         wfitn = sqrt(1 / 0.14 / 2),
         sig.e = 0,
         pos.p = 0.5,
         alpha = 0)

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

all.vlo.n %>%
  group_by(trial) %>% mutate(extinct = any(!n)) %>%
  group_by(n.pop0, gen) %>% mutate(nbar = mean(n), pext = mean(n>0)) %>% ungroup() %>%
  filter(n > 1) %>%
  ggplot() +
  geom_line(aes(x = gen, y = n, group = trial, linetype = factor(n.pop0), colour = extinct)) +
  scale_color_manual(values = c('black', 'red')) +
  geom_line(aes(x = gen, y = nbar, group = n.pop0), size = 2) +
  geom_point(aes(x = gen, y = nbar, fill = pext), shape = 21, stroke = 1, size = 3) +
  scale_fill_gradient(high = 'black', low = 'red') +
  scale_y_log10()

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

all.vhi.n %>%
  group_by(trial) %>% mutate(extinct = any(!n)) %>%
  group_by(n.pop0, gen) %>% mutate(nbar = mean(n), pext = mean(n>0)) %>% ungroup() %>%
  filter(n > 1) %>%
  ggplot() +
  geom_line(aes(x = gen, y = n, group = trial, linetype = factor(n.pop0), colour = extinct)) +
  scale_color_manual(values = c('black', 'red')) +
  geom_line(aes(x = gen, y = nbar, group = n.pop0), size = 2) +
  geom_point(aes(x = gen, y = nbar, fill = pext), shape = 21, stroke = 1, size = 3) +
  scale_fill_gradient(high = 'black', low = 'red') +
  scale_y_log10()

###

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

all.n %>%
  ggplot() +
  geom_ribbon(aes(x = gen, 
                  ymin = nbar - 2*sqrt(nvar / nrow(pars)/2), 
                  ymax = nbar + 2*sqrt(nvar / nrow(pars)/2),
                  group = interaction(n.pop0, bottleneck)),
              alpha = 0.5) +
  geom_line(aes(x = gen, y = n, group = trial, colour = extinct), size = 0.15) +
  scale_color_manual(values = c('black', 'red')) +
  geom_line(aes(x = gen, y = nbar, group = n.pop0), size = 2) +
  geom_point(aes(x = gen, y = nbar, fill = 1-pext), shape = 21, stroke = 1, size = 3) +
  scale_fill_gradient(high = 'black', low = 'red') +
  scale_y_log10() +
  facet_wrap(factor(n.pop0) ~ bottleneck) +
  guides(colour = 'none') +
  theme(legend.position = 'bottom')

all.n %>%
  ggplot() +
  geom_ribbon(aes(x = gen, 
                  ymin = nbar - 2*sqrt(nvar / nrow(pars)/2), 
                  ymax = nbar + 2*sqrt(nvar / nrow(pars)/2),
                  group = interaction(n.pop0, bottleneck)),
              alpha = 0.5) +
  scale_color_manual(values = c('black', 'red')) +
  geom_line(aes(x = gen, y = nbar, 
                group = interaction(n.pop0, bottleneck), 
                linetype = bottleneck), size = 2) +
  geom_point(aes(x = gen, y = nbar, fill = 1-pext), shape = 21, stroke = 1, size = 3) +
  scale_fill_gradient(high = 'black', low = 'red') +
  scale_y_log10() +
  facet_wrap(~ factor(n.pop0)) +
  guides(colour = 'none') +
  theme(legend.position = 'bottom')

# Summary data frame
n.all.n = all.n %>%
  distinct(n.pop0, bottleneck, gen, nbar, nvar, pext) %>%
  mutate(pext = 1 - pext)

### Try again but this time with environmental phenotypic noise...

epars = expand.grid(pop.id = unique(neck.popns$trial),
                    pop.trial = 1:trials.per,
                    n.pop0 = c(20, 50)) %>%
  mutate(global.trial = 1:nrow(.),
         end.time = 15,
         init.row = 1e4,
         n.loci = 25,
         w.max = 2,
         theta = 3,
         wfitn = sqrt(1 / 0.14 / 2),
         sig.e = 0.2,
         pos.p = 0.5,
         alpha = 0)

eliszt.vlo = vector('list', nrow(pars))

set.seed(756758)

for (i in 1:nrow(pars)) {
  popid = epars$pop.id[i]
  n0 = epars$n.pop0[i]
  eliszt.vlo[[i]] = sim( a = c(1/2, -1/2), 
                         params = epars[i,],
                         init.popn = neck.popns %>%
                           filter(trial %in% popid ) %>% # [neck.popns$trial %in% pars$pop.id[i],] %>%
                           sample_n(n0, replace = FALSE) %>%
                           select(-trial)
  )
  print(i)
}

e.vlo.trials = unroller(eliszt.vlo)


eliszt.vhi = vector('list', nrow(pars))

set.seed(9594)

for (i in 1:nrow(epars)) {
  eliszt.vhi[[i]] = sim( params = epars[i,])
  print(i)
}

e.vhi.trials = unroller(eliszt.vhi)

e.all.n = rbind(
  e.vlo.trials %>%
    group_by(trial, gen) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    rbind(expand.grid(trial = 1:nrow(pars), gen = 1:pars$end.time[1], n = 0)) %>%
    group_by(trial, gen) %>%
    summarise(n = sum(n)) %>%
    ungroup() %>%
    merge(y = pars, by.x = 'trial', by.y = 'global.trial') %>%
    mutate(bottleneck = TRUE),
  e.vhi.trials %>%
    group_by(trial, gen) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    rbind(expand.grid(trial = 1:nrow(pars), gen = 1:pars$end.time[1], n = 0)) %>%
    group_by(trial, gen) %>%
    summarise(n = sum(n)) %>%
    ungroup() %>%
    merge(y = pars, by.x = 'trial', by.y = 'global.trial') %>%
    mutate(bottleneck = FALSE,
           trial = trial + max(e.vlo.trials$trial))
) %>%
  group_by(n.pop0, bottleneck, gen) %>% 
  summarise(nbar = mean(n), 
            nvar = var(n),
            pext = mean(!n))

e.all.n %>%
  ggplot() +
  geom_ribbon(aes(x = gen, 
                  ymin = nbar - 2*sqrt(nvar / nrow(pars)/2), 
                  ymax = nbar + 2*sqrt(nvar / nrow(pars)/2),
                  group = interaction(n.pop0, bottleneck)),
              alpha = 0.5) +
  scale_color_manual(values = c('black', 'red')) +
  geom_line(aes(x = gen, y = nbar, 
                group = interaction(n.pop0, bottleneck), 
                linetype = bottleneck), size = 2) +
  geom_point(aes(x = gen, y = nbar, fill = pext), shape = 21, stroke = 1, size = 3) +
  scale_fill_gradient(high = 'red', low = 'black') +
  scale_y_log10() +
  facet_wrap(~ factor(n.pop0)) +
  guides(colour = 'none') +
  theme(legend.position = 'bottom')

rbind(n.all.n %>% ungroup() %>% select(names(e.all.n)) %>% mutate(sig.e = 0),
      e.all.n %>% ungroup() %>% mutate(sig.e = 0.2)) %>%
  ggplot() +
  geom_ribbon(aes(x = gen, 
                  ymin = nbar - 2*sqrt(nvar / nrow(pars)/2), 
                  ymax = nbar + 2*sqrt(nvar / nrow(pars)/2),
                  group = interaction(n.pop0, bottleneck, sig.e)),
              alpha = 0.5) +
  geom_line(aes(x = gen, y = nbar, 
                group = interaction(n.pop0, bottleneck, sig.e), 
                linetype = bottleneck), 
            size = 1) +
  geom_point(aes(x = gen, y = nbar, fill = pext, shape = factor(sig.e)), 
             stroke = 1, size = 3) +
  scale_fill_gradient(high = 'red', low = 'black') +
  scale_shape_manual(values = c(21, 22)) +
  scale_y_log10() +
  facet_wrap( ~ factor(n.pop0)) +
  guides(colour = 'none') +
  theme(legend.position = 'bottom')

# Export for use later
rbind(n.all.n %>% ungroup() %>% select(names(e.all.n)) %>% mutate(sig.e = 0),
      e.all.n %>% ungroup() %>% mutate(sig.e = 0.2)) %>%
  write.csv(file = 'analyze_model/bottlenecks/triple_bottlenecks_summstats.csv',
            row.names = FALSE, na = '')
