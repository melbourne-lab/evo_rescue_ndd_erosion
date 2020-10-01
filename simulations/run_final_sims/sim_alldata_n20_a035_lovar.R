# Script for simulating population dynamics of an adapting bottlenecked populations.
# Here, I am extracting all statistics.
# In this script: low genetic variation, initially small, density dependent.
# SN - adapted from simulations/run_final_simulations/test_alldata_n20_a000_lovar.R
# adapted and run on September 22 2020

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
pars = data.frame(
  trial = 1:trials,
  n.pop0 = 20,
  end.time = 15,
  init.row = 1e4,
  n.loci = 25,
  w.max = 2,
  theta = 2.75,
  wfitn = sqrt(1 / 0.14 / 2),
  sig.e = sqrt(0.5),
  pos.p = 0.5,
  alpha = 0.0035
)

liszt = vector('list', nrow(pars))

set.seed(9228324)

for (i in 1:nrow(pars)) {
  pop.init = init.sim(a = c(1/2, -1/2),
                      params = pars[i,]) %>%
    mutate_at(paste0('a', 1:6), function(x) -1/2) %>%
    mutate_at(paste0('b', 1:6), function(x) -1/2) %>%
    mutate_at(paste0('a', 7:12), function(x) 1/2) %>%
    mutate_at(paste0('b', 7:12), function(x) 1/2)
  sim.output = sim( a = c(1/2, -1/2),
                    params = pars[i,],
                    init.popn = pop.init)
  
  demo.summ = sim.output %>%
    group_by(gen) %>%
    summarise(n = n(),
              gbar = mean(g_i),
              zbar = mean(z_i),
              wbar = mean(w_i))
  
  gene.summ = sim.output %>%
    select(-c(i, g_i, z_i, w_i, r_i, fem)) %>%
    gather(key = loc.copy, value = val, -gen) %>%
    mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
    group_by(gen, locus) %>%
    summarise(p = mean(val > 0)) %>%
    group_by(gen) %>%
    summarise(p.fix.pos = mean(p == 1),
              p.fix.neg = mean(p == 0),
              v = sum(2 * p * (1 - p)) / pars$n.loci[1])
  
  liszt[[i]] = cbind(demo.summ, gene.summ %>% select(-gen)) 
  
  print(paste0('ndd 20 lo var ', i, ' of ', nrow(pars)))
}

all.sims = unroller(liszt) %>%
  merge(y = pars %>% select(trial, n.pop0, alpha) %>% mutate(low.var = TRUE), by = 'trial') %>%
  group_by(trial)  %>%
  mutate(ext.gen = max(gen),
         extinct = !ext.gen %in% pars$end.time[1]) %>%
  ungroup()

write.csv(all.sims, row.names = FALSE,
          file = "simulations/outputs/alldata_n20_a035_lowvar.csv")
