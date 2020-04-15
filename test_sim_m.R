# Generate 1000 simulation runs on server.
# Run April 8, 2020
# (note: will be obsolete as density dependence is added)

library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

rm(list = ls())

setwd('~/demall')

source('base_model_source/sim_functions.R')

n.trials = 8000

mlist = rep(1:8, each = n.trials / 8)
liszt = vector(mode = 'list', length = n.trials)

set.seed(603185)

for (trial in 1:n.trials) {
  liszt[[trial]] = sim(
    a = c(-1/2, 1/2),
    params = data.frame(end.time = 15,
                        init.row = 1e4,
                        n.loci = mlist[trial]^2, 
                        n.pop0 = 20,
                        w.max = 2, 
                        theta = 2.6, 
                        wfitn = sqrt(1 / 0.14 / 2),
                        sig.e = 0)
  ) %>%
    mutate(n.loci = mlist[trial]^2)
  print(trial)
}

for (m in 1:8) {
  assign(paste0('sim.sum.m', m), 
         unroller(liszt[(m-1)*(n.trials/8) + 1:(n.trials/8)]) %>%
           select(trial, names(.)[!grepl('trial', names(.))]) %>%
           select(-c(r_i, fem)) %>%
           gather(key = loc.copy, value = val, -c(trial, i, g_i, z_i, w_i, gen, n.loci)) %>%
           mutate(loc = gsub('[ab]', '', loc.copy) %>% as.numeric()) %>%
           group_by(trial, gen, loc, n.loci) %>%
           summarise(n = length(unique(i)),
                     g = mean(g_i),
                     z = mean(z_i),
                     w = mean(w_i),
                     p = mean(val > 0),
                     v = 2 * p * (1 - p) / m^2,
                     fixp = all(val > 0),
                     fixn = all(val < 0)) %>%
           group_by(trial, gen) %>%
           summarise(pop.size = mean(n),
                     n.fixp = sum(fixp),
                     n.fixn = sum(fixn),
                     g.varn = sum(v),
                     mean.p = mean(p),
                     mean.g = mean(g),
                     mean.z = mean(z),
                     mean.w = mean(w)) %>%
           group_by(gen) %>%
           summarise(n.trials = n(),
                     n.bar = mean(pop.size),
                     n.var = var(pop.size, na.rm = TRUE),
                     g.bar = mean(mean.g),
                     g.var = var(mean.g),
                     w.bar = mean(mean.w),
                     w.var = var(mean.w, na.rm = TRUE),
                     p.bar = mean(mean.p),
                     v.bar = mean(g.varn),
                     v.var = var(g.varn, na.rm = TRUE),
                     p.fixp = mean(n.fixp) / m^2,
                     p.fixn = mean(n.fixn) / m^2) %>%
           mutate(m = m^2)
  )
}

mget(grep('sim.sum.m', ls(), value = TRUE)) %>%
  do.call(what = rbind) %>%
  write.csv('out/test_m_k.csv',
            row.names = FALSE)

