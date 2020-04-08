# Generate 1000 simulation runs on server.
# Run April 8, 2020
# (note: will be obsolete as density dependence is added)

library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

rm(list = ls())

source('base_model_source/sim_functions.R')

n.trials = 1000

liszt = vector(mode = 'list', length = n.trials)

set.seed(1581823)

for (trial in 1:n.trials) {
    liszt[[trial]] = sim(
        a = c(-1/2, 1/2),
        params = data.frame(end.time = 15,
                            init.row = 1e4,
                            n.loci = 25, 
                            n.pop0 = 20,
                            w.max = 2, 
                            theta = 2.6, 
                            wfitn = sqrt(1 / 0.14 / 2),
                            sig.e = 0)
    )
    print(trial)
}
    

all.sims = unroller(liszt) %>%
    select(trial, names(.)[!grepl('trial', names(.))])

sumry = all.sims %>%
    select(-c(r_i, fem)) %>%
    gather(key = loc.copy, value = val, -c(trial, i, g_i, z_i, w_i, gen)) %>%
    mutate(loc = gsub('[ab]', '', loc.copy) %>% as.numeric()) %>%
    group_by(trial, gen, loc) %>%
    summarise(n = length(unique(i)),
              g = mean(g_i),
              z = mean(z_i),
              w = mean(w_i),
              p = mean(val > 0),
              v = 2 * p * (1 - p) / 25,
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
              p.fixp = mean(n.fixp) / 25,
              p.fixn = mean(n.fixn) / 25)

sumry

sumry %>%
    write.csv('out/eg_1k_n0_20.csv',
              row.names = FALSE)

