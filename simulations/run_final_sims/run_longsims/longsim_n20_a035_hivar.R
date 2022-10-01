# Script for simulating population dynamics of an adapting bottlenecked populations.
# I am running these simulations for 50 generations to get extincton data.
# In this script: high genetic variation, initially small, density dependent.
# SN - adapted from simulations/run_final_simulations/longsim/longsim_alldata_n20_a035_lowvar.R
# adapted and run on March 3 2021

### Clear namespace
rm(list = ls())

### Load packages

library(ggplot2)
library(dplyr)
library(tidyr)

### Load source materials

# Get simulation functions
source('sim_functions.R')

# Define number of trials
trials = 1000

# Define parameters
pars = data.frame(
  size.thresh = 10000,
  init.row = 1e6,
  end.time = 51,
  n.pop0 = 20,
  n.loci = 25,
  trial = 1:trials,
  w.max = 2,
  theta = 2.8,
  wfitn = sqrt(3.5),
  sig.e = sqrt(0.5),
  pos.p = 0.5,
  alpha = 0.0035
)

liszt = vector('list', nrow(pars))

set.seed(420820)

for (i in 1:nrow(pars)) {
  
  sim.output = sim(a = c(1/2, -1/2), params = pars[i,])
  
  demo.summ = sim.output %>%
    group_by(gen) %>%
    summarise(
      n = n(),
      gbar = mean(g_i),
      zbar = mean(z_i),
      wbar = mean(w_i),
      pfem = mean(fem)
    )  %>%
    ungroup()
  
  gene.summ = sim.output %>%
    select(-c(i, g_i, z_i, w_i, r_i, fem)) %>%
    gather(key = loc.copy, value = val, -gen) %>%
    mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
    group_by(gen, locus) %>%
    summarise(p = mean(val > 0)) %>%
    group_by(gen) %>%
    summarise(
      p.fix.pos = mean(p == 1),
      p.fix.neg = mean(p == 0),
      v = sum(2 * p * (1 - p)) / pars$n.loci[1]
    ) %>%
    ungroup()
  
  liszt[[i]] = cbind(demo.summ, gene.summ %>% select(-gen)) %>%
    mutate(trial = i)
  
  print(paste0('ndd 20 hi var ', i, ' of ', nrow(pars)))
}

merge(
  x = do.call(rbind, liszt),
  y = pars %>% select(trial, n.pop0, alpha) %>% mutate(low.var = FALSE), 
  by = 'trial'
) %>%
  write.csv(file = "simulations/outputs/longsims/longsims_n20_a035_hivar.csv",
            row.names = FALSE)
