library(ggplot2)
library(dplyr)
library(tidyr)

kpops = read.csv('simulations/prepare_simulations/outputs/all_k_pops_min_gen.csv')

# How many rows (inidividuals total) are there?
nrow(kpops)
# How many different generations represented?
unique(kpops$gen)
# How many individual/unique trials are there?
length(unique(kpops$trial))
# How large is each population?
table(kpops$trial)
# Most of these are ~1000

# Here - let me assess the genetic variation in each population.

kvars = kpops %>%
  select(-c(g_i, z_i, w_i, r_i, fem, gen)) %>%
  gather(key = loc.copy, value = val, -c(i, trial)) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(trial, locus) %>%
  summarise(n = n() / 2,
            p = mean(val > 0)) %>%
  group_by(trial) %>%
  summarise(n = n[1],
            v = sum(2 * p * (1 - p) / 25))

kvars  
ggplot(kvars) + 
  geom_histogram(aes(x = v), binwidth = 0.025) + 
  geom_segment(aes(x = 0.48, xend = 0.48, y = 0, yend = 3), colour = 'red') +
  xlim(c(0, 0.5))

ggplot(kvars) + 
  geom_point(aes(x = n, y = v))

# Get long-form kpop data frame
klong = kpops %>%
  select(-c(g_i, z_i, w_i, r_i, fem, gen)) %>%
  gather(key = loc.copy, value = val, -c(i, trial)) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy))

datafram = expand.grid(popn = unique(kvars$trial),
                       rep = 1:100,
                       n0 = c(20, 100),
                       v = NA)

set.seed(32300)

for (i in 1:nrow(datafram)) {
  
  pop.init = klong[klong$trial %in% datafram$popn[i],]
  
  datafram$v[i] = pop.init %>%
    filter(i %in% sample(unique(pop.init$i), 
                         size = datafram$popn[i], 
                         replace = FALSE)) %>%
    group_by(locus) %>%
    summarise(p = mean(val > 0)) %>%
    ungroup() %>%
    summarise(v = sum(2 * p * (1 - p) / 25)) %>%
    unlist()
  
  print(i)
  
}

datafram = datafram %>%
  mutate(v = unlist(v)) %>%
  rename(v.samp = v)
  merge(y = kvars, by.x = 'popn', by.y = 'trial')

datafram %>%
  ggplot(aes(x = factor(popn), y = v.samp)) +
  geom_jitter() +
  geom_point(aes(x = factor(popn), y = v), 
             colour = 'red',
             shape = '+') +
  facet_wrap(~ n0)

datafram %>%
  ggplot() +
  geom_violin(aes(x = factor(popn), y = v.samp)) +
  geom_point(aes(x = factor(popn), y = v), 
             colour = 'blue',
             shape = '+',
             size = 5) +
  facet_wrap(~ n0)

datafram %>%
  mutate(v = unlist(v)) %>%
  group_by(popn, n0) %>%
  summarise(vlo = quantile(v, 0.025),
            vhi = quantile(v, 0.975)) %>%
  ggplot() +
  geom_segment(aes(x = factor(popn), xend = factor(popn),
                   y = vlo, yend = vhi, group = n0),
               position = position_dodge(width = 0.25))
# Stuff is messed up here. Who cares.

