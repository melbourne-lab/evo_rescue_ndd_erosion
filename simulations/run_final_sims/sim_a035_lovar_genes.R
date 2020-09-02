# Script for simulating population dynamics of an adapting bottlenecked populations.
# Here, I am extracting mean genetic variance and fixation probabilities.
# In this script: low genetic variation, density dependent.
# SN - adapted from simulations/run_final_simulations/sim_a000_hivar_genes.R
# adapted and run on September 1 2020

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

# Define number of trials
trials = 4000

# Define parameters
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
         alpha = 0.0035)

liszt.vlo = vector('list', nrow(pars))

set.seed(726153)

for (i in 1:nrow(pars)) {
  pop.init = init.sim(a = c(1/2, -1/2),
                      params = pars[i,]) %>%
    mutate_at(paste0('a', 1:6), function(x) -1/2) %>%
    mutate_at(paste0('b', 1:6), function(x) -1/2) %>%
    mutate_at(paste0('a', 7:12), function(x) 1/2) %>%
    mutate_at(paste0('b', 7:12), function(x) 1/2)
  liszt.vlo[[i]] = sim( a = c(1/2, -1/2),
                        params = pars[i,],
                        init.popn = pop.init) %>%
    select(-c(i, g_i, z_i, w_i, r_i, fem)) %>%
    gather(key = loc.copy, value = val, -gen) %>%
    mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
    group_by(gen, locus) %>%
    summarise(p = mean(val > 0)) %>%
    group_by(gen) %>%
    summarise(p.fix.pos = mean(p == 1),
              p.fix.neg = mean(p == 0),
              v = sum(2 * p * (1 - p)) / pars$n.loci[1])
  print(paste0('ndd lo var ', i, ' of ', nrow(pars)))
}

all.vlo.genes = unroller(liszt.vlo) %>%
  merge(y = pars %>% select(global.trial, n.pop0), 
        by.x = 'trial', by.y = 'global.trial') %>%
  group_by(trial)  %>%
  mutate(ext.gen = max(gen),
         extinct = !ext.gen %in% pars$end.time[1]) %>%
  ungroup()

### Combine the results.

gene.summ.all = all.vlo.genes %>%
  group_by(n.pop0, gen) %>%
  summarise(p.fix.pos = mean(p.fix.pos),
            p.fix.neg = mean(p.fix.neg),
            vbar = mean(v),
            vvar = var(v),
            n = n())

write.csv(gene.summ.all, na = '', row.names = FALSE,
          file = 'simulations/outputs/final_results/a035_lowvar_gene_all.csv')

gene.summ.ext = all.vlo.genes %>%
  group_by(n.pop0, gen, extinct) %>%
  summarise(p.fix.pos = mean(p.fix.pos),
            p.fix.neg = mean(p.fix.neg),
            vbar = mean(v),
            vvar = var(v),
            n = n())

write.csv(gene.summ.ext, na = '', row.names = FALSE,
          file = 'simulations/outputs/final_results/a035_lowvar_gene_ext.csv')

gene.summ.gen = all.vlo.genes %>%
  group_by(n.pop0, gen, ext.gen) %>%
  summarise(p.fix.pos = mean(p.fix.pos),
            p.fix.neg = mean(p.fix.neg),
            vbar = mean(v),
            vvar = var(v),
            n = n())

write.csv(gene.summ.gen, na = '', row.names = FALSE,
          file = 'simulations/outputs/final_results/a035_lowvar_gene_gen.csv')
