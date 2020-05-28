# Generate 1000 simulation runs on server.
# Run May 28, 2020
# (note: will be obsolete as density dependence is added)

# This script is designed to evaluate simulation phenotypic changes with the
# expectations of Lande, 1976.

library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

rm(list = ls())

setwd('~/demall')

source('base_model_source/sim_functions.R')

# Global params used in all sims.
theta = 2
wfitn = sqrt(1 / 0.14 / 2)

# For initializing and running simualtion
n.trials = 10000

mlist = rep(c(1, 3, 5, 7, 9), each = n.trials / 5)
liszt = vector(mode = 'list', length = n.trials)

set.seed(404303)

for (trial in 1:n.trials) {
  liszt[[trial]] = sim(
    a = c(-1/2, 1/2),
    params = data.frame(end.time = 2,
                        init.row = 1e2,
                        n.loci = mlist[trial]^2, 
                        n.pop0 = 20,
                        w.max = 2, 
                        theta = theta, 
                        wfitn = wfitn,
                        sig.e = 0)
  ) %>%
    mutate(n.loci = mlist[trial]^2)
  print(trial)
}

for (m in 1:5) {
  assign(paste0('sims.m', m), 
         unroller(liszt[(m-1)*(n.trials/5) + 1:(n.trials/5)]) %>%
           select(trial, names(.)[!grepl('trial', names(.))]) %>%
           select(-c(g_i, w_i, r_i, fem)) %>%
           gather(key = loc.copy, value = val, 
                  -c(trial, i, z_i, gen, n.loci)) %>%
           mutate(loc = gsub('[ab]', '', loc.copy) %>% as.numeric()) %>%
           group_by(trial, gen, loc) %>%
           summarise(n_t = length(unique(i)),
                     p_j = mean(val > 0),
                     z_i = mean(z_i),
                     n.loci = n.loci[1]) %>%
           group_by(trial, gen) %>%
           summarise(n_t = n_t[1],
                     z_t = mean(z_i),
                     v_t = sum(2 * p_j * (1-p_j) / n.loci),
                     n_loci = n.loci[1]) %>%
           mutate(d_t = theta - z_t,
                  lande_d_t = d_t * (wfitn^2 / (wfitn^2 + v_t))) %>%
           select(trial, gen, n_t, z_t, v_t, d_t, lande_d_t, n_loci)
  )
}

mget(grep('sims.m', ls(), value = TRUE)) %>%
  do.call(what = rbind) %>%
  write.csv('out/first_step_lande.csv',
            row.names = FALSE)

