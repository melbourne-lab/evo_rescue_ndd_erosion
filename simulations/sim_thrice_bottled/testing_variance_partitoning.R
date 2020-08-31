# Script for simulating population dynamics of an adapting bottlenecked populations.
# Here, I am extracting mean genotype, phenotype, and fitness for all eight combos.
# In this script: low genetic variation (bottlenecks), density dependent.
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
trials.per = 20 # 75

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
         alpha = 0.0035)

liszt.vlo = vector('list', nrow(pars))

set.seed(853443)

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
    group_by(gen) %>%
    summarise(n = n(),
              g_t = mean(g_i),
              z_t = mean(z_i),
              w_t = mean(w_i))
  print(paste0('test ndd lo var ', i, ' of ', nrow(pars)))
}

all.vlo.trials = unroller(liszt.vlo) %>%
  group_by(trial) %>%
  mutate(n.pop0 = n[1],
         ext.gen = max(gen),
         extinct = !ext.gen %in% pars$end.time[1]) %>%
  ungroup() %>%
  merge(y = pars %>% select(pop.id, global.trial),
        by.x = 'trial', by.y = 'global.trial')

all.vlo.trials %>%
  group_by(gen, pop.id, n.pop0) %>%
  summarize(g_win_pop = mean(g_t),
            v_win_pop = var(g_t)) %>%
  group_by(gen, n.pop0) %>%
  summarise(v_amg_pop = var(g_win_pop),
            v = mean(v_win_pop, na.rm = TRUE) + v_amg_pop)

# 

# compare with...

all.vlo.trials %>%
  group_by(gen, n.pop0) %>%
  summarise(gvar = var(g_t))
  
### Combine the results.

type.summ.all.means =  all.vlo.trials %>%
  group_by(n.pop0, gen) %>%
  summarise(gbar = mean(g_t),
            zbar = mean(z_t),
            wbar = mean(w_t),
            n = n())

type.summ.all.vars = all.vlo.trials %>%
  group_by(gen, n.pop0, pop.id) %>%
  summarize(gbar_win_pop = mean(g_t),
            gvar_win_pop = var(g_t),
            zbar_win_pop = mean(z_t),
            zvar_win_pop = var(z_t),
            wbar_win_pop = mean(w_t),
            wvar_win_pop = var(w_t)) %>%
  group_by(n.pop0, gen) %>%
  summarise(gvar_win = mean(gvar_win_pop, na.rm = TRUE),
            gvar_amg = var(gbar_win_pop, na.rm = TRUE),
            gvar_var = var(gvar_win_pop, na.rm = TRUE),
            zvar_win = mean(zvar_win_pop, na.rm = TRUE),
            zvar_amg = var(zbar_win_pop, na.rm = TRUE),
            zvar_var = var(zvar_win_pop, na.rm = TRUE),
            wvar_win = mean(wvar_win_pop, na.rm = TRUE),
            wvar_amg = var(wbar_win_pop, na.rm = TRUE),
            wvar_var = var(wvar_win_pop, na.rm = TRUE))

type.summ.all = merge(type.summ.all.means, type.summ.all.vars,
                      by = c('n.pop0', 'gen'))

write.csv(type.summ.all, na = '', row.names = FALSE,
          file = 'simulations/outputs/thrice_bottlenecked_results/a035_lovar_type_all_test.csv')

type.summ.ext.means = all.vlo.trials %>%
  group_by(n.pop0, extinct, gen) %>%
  summarise(gbar = mean(g_t),
            zbar = mean(z_t),
            wbar = mean(w_t),
            n = n())

type.summ.ext.vars = all.vlo.trials %>%
  group_by(gen, n.pop0, pop.id, extinct) %>%
  summarize(gbar_win_pop = mean(g_t),
            gvar_win_pop = var(g_t),
            zbar_win_pop = mean(z_t),
            zvar_win_pop = var(z_t),
            wbar_win_pop = mean(w_t),
            wvar_win_pop = var(w_t)) %>%
  group_by(n.pop0, extinct, gen) %>%
  summarise(gvar_win = mean(gvar_win_pop, na.rm = TRUE),
            gvar_amg = var(gbar_win_pop, na.rm = TRUE),
            gvar_var = var(gvar_win_pop, na.rm = TRUE),
            zvar_win = mean(zvar_win_pop, na.rm = TRUE),
            zvar_amg = var(zbar_win_pop, na.rm = TRUE),
            zvar_var = var(zvar_win_pop, na.rm = TRUE),
            wvar_win = mean(wvar_win_pop, na.rm = TRUE),
            wvar_amg = var(wbar_win_pop, na.rm = TRUE),
            wvar_var = var(wvar_win_pop, na.rm = TRUE))

type.summ.ext = merge(type.summ.ext.means, type.summ.ext.vars,
                      by = c('n.pop0', 'gen', 'extinct'))

# what to do about NaNs?

write.csv(type.summ.ext, na = '', row.names = FALSE,
          file = 'simulations/outputs/thrice_bottlenecked_results/a035_lovar_type_ext_test.csv')


# Save all simulation output
# save(liszt.vhi, liszt.vlo, file = 'out/bottleneck_adapted_outputs.RData')

type.summ.gen.means = all.vlo.trials %>%
  group_by(n.pop0, ext.gen, gen) %>%
  summarise(gbar = mean(g_t),
            gvar = var(g_t),
            zbar = mean(z_t),
            zvar = var(z_t),
            wbar = mean(w_t),
            wvar = var(w_t),
            n = n())

type.summ.gen.vars = all.vlo.trials %>%
  group_by(gen, n.pop0, pop.id, ext.gen) %>%
  summarize(gbar_win_pop = mean(g_t),
            gvar_win_pop = var(g_t),
            zbar_win_pop = mean(z_t),
            zvar_win_pop = var(z_t),
            wbar_win_pop = mean(w_t),
            wvar_win_pop = var(w_t)) %>%
  group_by(n.pop0, ext.gen, gen) %>%
  summarise(gvar_win = mean(gvar_win_pop, na.rm = TRUE),
            gvar_amg = var(gbar_win_pop, na.rm = TRUE),
            gvar_var = var(gvar_win_pop, na.rm = TRUE),
            zvar_win = mean(zvar_win_pop, na.rm = TRUE),
            zvar_amg = var(zbar_win_pop, na.rm = TRUE),
            zvar_var = var(zvar_win_pop, na.rm = TRUE),
            wvar_win = mean(wvar_win_pop, na.rm = TRUE),
            wvar_amg = var(wbar_win_pop, na.rm = TRUE),
            wvar_var = var(wvar_win_pop, na.rm = TRUE))

type.summ.gen = merge(type.summ.gen.means, type.summ.gen.vars,
                      by = c('n.pop0', 'gen', 'ext.gen')) %>%
  arrange(n.pop0, gen, ext.gen)

write.csv(type.summ.gen, na = '', row.names = FALSE,
          file = 'simulations/outputs/thrice_bottlenecked_results/a035_lovar_type_gen_test.csv')

### Calculating errors with lme4 (because I can't understand Pinheiro and Bates)

all.vlo.trials %>% 
  filter(n.pop0 %in% 60) %>% 
  lmer(formula = z_t ~ factor(gen) + (1 | pop.id)) %>%
  predict(newdata = expand.grid(gen = 1:15),
          re.form = NA)
