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

# grid.y.res = .005
# 
# gridded.v = all.vlo.genes %>% 
#   mutate(v_rounded = round((1/grid.y.res) * v) / (1/grid.y.res)) %>%
#   group_by(v_rounded, gen, n.pop0) %>%
#   summarise(p.extinct = mean(extinct),
#             n = n())

write.csv(all.vlo.genes, na = '', row.names = FALSE,
          file = 'simulations/outputs/test_gridded_v.csv')

# ### Test plot
# 
# testv = read.csv('simulations/outputs/test_gridded_v.csv')
# 
# testv %>%
#   ggplot() +
#   geom_tile(
#     aes(
#       x = gen,
#       y = v_rounded,
#       fill = p.extinct,
#       alpha = log(n + 1, base = max(n) + 1)
#     )
#   )

bigboy = read.csv('simulations/outputs/test_gridded_v.csv')

nrow(bigboy)

bigboy %>%
  ggplot(aes(x = gen)) +
  geom_point(
    aes(
      y = v,
      colour = extinct
    ),
    position = position_jitter(width = 0.5),
    alpha = 0.1
  ) +
  facet_wrap(~ n.pop0)

bigboy %>%
  ggplot(aes(x = gen)) +
  geom_point(
    aes(
      y = v,
      colour = extinct
    ),
    position = position_jitter(width = 0.5),
    alpha = 0.15,
    size = 3
  ) +
  scale_color_manual(values = c('royalblue', 'red')) +
  facet_wrap(~ n.pop0, ncol = 1) +
  theme_bw()

bigboy %>%
  ggplot(aes(x = gen)) +
  geom_line(
    data = . %>% filter(ext.gen < 15),
    aes(
      y = v,
      group = trial,
      colour = ext.gen
    ),
    size = 0.1
  ) +
  geom_line(
    data = . %>% filter(ext.gen %in% 15),
    aes(
      y = v,
      group = trial
    ),
    colour = 'black',
    size = 0.1
  ) +
  geom_point(
    data = . %>% filter(ext.gen < 15),
    aes(
      y = v,
      colour = ext.gen
    ),
    position = position_jitter(width = 0.5),
    alpha = 0.1,
    size = 2
  ) +
  geom_point(
    data = . %>% filter(ext.gen %in% 15),
    aes(
      y = v,
    ),
    position = position_jitter(width = 0.5),
    colour = 'black',
    alpha = 0.025,
    size = 2
  ) +
  scale_color_gradient(low = 'lightpink', high = 'darkred')

bigboy %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = v,
      group = trial
    ),
    colour = 'black',
    size = 0.01
  ) +
  geom_point(
    aes(
      y = v,
      colour = extinct
    ),
    position = position_jitter(width = 0.5),
    alpha = 0.1,
    size = 2
  ) +
  scale_color_manual(values = c('royalblue', 'red')) +
  facet_wrap(~ n.pop0, ncol = 1)
