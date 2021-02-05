# Visualizing the extinction vortex in a non-evolving (stationary) population.
# Here, populations are at the phenotypic optimum but are initially small enough
# (and have low enough max. fitness) that they will go extinct.
# SN - February 4 2021

library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

rm(list = ls())

source('base_model_source/sim_functions.R')

pars = data.frame(
  init.row = 1e4,
  end.time = 15,
  n.pop0 = 40,
  n.loci = 25,
  w.max = 1,
  theta = 0,
  wfitn = sqrt(1 / 0.14 / 2),
  sig.e = sqrt(0.5),
  pos.p = 0.5,
  alpha = 0
)

n.trials = 400

liszt = vector('list', n.trials)

set.seed(992)

for (trial in 1:n.trials) {
  
  liszt[[trial]] = sim(params = pars)
  print(trial)
  
}

all.trials = liszt %>% unroller()

nrow(all.trials)

max.gen = pars$end.time

all.summ = merge(
  x = all.trials %>%
    group_by(trial, gen) %>%
    summarise(n = n(), wbar = mean(w_i), zbar = mean(z_i)),
  y = all.trials %>%
    select(-c(i, g_i, z_i, w_i, r_i, fem)) %>%
    gather(key = loc.copy, value = val, -c(gen, trial)) %>%
    mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
    group_by(trial, gen, locus) %>%
    summarise(p = mean(val > 0)) %>%
    group_by(trial, gen) %>%
    summarise(v = sum(2 * p * (1 - p)) / pars$n.loci[1]),
  by = c('trial', 'gen')
) %>%
  arrange(trial, gen)
  

all.summ.extended = all.summ %>%  
  group_by(trial) %>%
  arrange(desc(gen)) %>%
  distinct(trial, .keep_all = TRUE) %>%
  uncount(weight = max.gen - gen) %>%
  mutate(n = 0, gen = gen + (1:(max.gen-gen[1]))) %>%
  ungroup() %>% 
  rbind(all.summ) %>%
  group_by(trial) %>%
  mutate(extinct = any(!n)) %>% 
  group_by(gen, extinct) %>%
  summarise(n.bar = mean(n),
            n.var = var(n),
            v.bar = mean(v),
            v.var = var(v),
            w.bar = mean(wbar),
            w.var = var(wbar),
            z.bar = mean(zbar),
            z.var = var(zbar),
            nobs = n()
  ) %>%
  ungroup() %>%
  filter(!(extinct & gen %in% max.gen))

print(all.summ.extended, nrow = 40)

nv.plot = ggplot(all.summ.extended, aes(x = v.bar, y = n.bar)) +
  geom_path(
    aes(
      group = extinct,
      colour = extinct
    )
  ) +
  scale_y_log10() +
  scale_color_manual(values = c('black', 'red')) +
  labs(x = 'Additive variance', y = 'Population size')

nw.plot = ggplot(all.summ.extended, aes(x = w.bar, y = n.bar)) +
  geom_path(
    aes(
      group = extinct,
      colour = extinct
    )
  ) +
  scale_y_log10() +
  scale_color_manual(values = c('black', 'red')) +
  labs(x = 'Mean intrinsic fitness', y = 'Population size')
# what?

vw.plot = ggplot(all.summ.extended, aes(x = w.bar, y = v.bar)) +
  geom_path(
    aes(
      group = extinct,
      colour = extinct
    )
  ) +
  scale_color_manual(values = c('black', 'red')) +
  labs(x = 'Mean intrinsic fitness', y = 'Additive fitness')
# what...

wt.plot = ggplot(all.summ.extended, aes(x = gen, y = w.bar)) +
  geom_ribbon(
    aes(
      ymin = w.bar - 2 * sqrt(w.var / nobs),
      ymax = w.bar + 2 * sqrt(w.var / nobs),
      group = extinct
    ),
    alpha = 0.1
  ) +
  geom_path(
    aes(
      group = extinct,
      colour = extinct
    )
  ) +
  scale_color_manual(values = c('black', 'red')) +
  labs(x = 'Time step', y = 'Mean intrinsic fitness')

all.summ %>%
  filter(trial < 31) %>%
  ggplot(aes(x = gen, y = wbar)) +
  geom_line(aes(group = trial))

plot_grid(nv.plot, nw.plot, vw.plot, wt.plot,
          nrow = 2)

# Next: plot by tau.
