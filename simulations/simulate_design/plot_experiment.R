# Plot results

library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

d1 = read.csv('simulations/simulate_design/output/design1_summary.csv')
d2 = read.csv('simulations/simulate_design/output/design2_summary.csv')
d3 = read.csv('simulations/simulate_design/output/design3_summary.csv')
nb = read.csv('simulations/simulate_design/output/design_nb_summary.csv')

# Combine I guess
# (me reaping... me sowing)

# Design 1:

design1 = rbind(d1, nb)

design1 %>%
  mutate(bottleneck = factor(design > 0,  labels = c('hi var', 'lo var')),
         n.pop0     = factor(n.pop0 > 20, labels = c('small', 'large')),
         alpha      = factor(alpha  > 0,  labels = c('dens indep', 'dens dep'))) %>%
  group_by(trial, n.pop0, bottleneck, alpha, evolve) %>%
  mutate(extinct = max(gen) < 15) %>%
  ungroup() %>%
  ggplot(aes(x = gen, y = n)) +
  geom_line(
    aes(
      group    = interaction(trial, n.pop0, evolve, bottleneck),
      linetype = evolve,
      size     = evolve,
      colour   = extinct
    )
  ) +
  scale_linetype_manual(values = 2:1) +
  scale_size_manual(values = c(0.25, 0.75)) +
  scale_color_manual(values = c('black', 'red')) +
  scale_y_log10() +
  facet_wrap(alpha ~ paste(n.pop0, bottleneck), ncol = 4) +
  theme(legend.position = 'bottom')

# Design 2

design2 = rbind(d2 %>% select(names(nb)), nb)

design2 %>%
  mutate(bottleneck = factor(design > 0,  labels = c('hi var', 'lo var')),
         n.pop0     = factor(n.pop0 > 20, labels = c('small', 'large')),
         alpha      = factor(alpha  > 0,  labels = c('dens indep', 'dens dep'))) %>%
  group_by(trial, n.pop0, bottleneck, alpha, evolve) %>%
  mutate(extinct = max(gen) < 15) %>%
  ungroup() %>%
  ggplot(aes(x = gen, y = n)) +
  geom_line(
    aes(
      group    = interaction(trial, n.pop0, evolve, bottleneck),
      linetype = evolve,
      size     = evolve,
      colour   = extinct
    )
  ) +
  scale_linetype_manual(values = 2:1) +
  scale_size_manual(values = c(0.25, 0.75)) +
  scale_color_manual(values = c('black', 'red')) +
  scale_y_log10() +
  facet_wrap(alpha ~ paste(n.pop0, bottleneck), ncol = 4) +
  theme(legend.position = 'bottom')

### Design 3

design3 = rbind(d3 %>% select(names(nb)), nb)

design3 %>%
  mutate(bottleneck = factor(design > 0,  labels = c('hi var', 'lo var')),
         n.pop0     = factor(n.pop0 > 20, labels = c('small', 'large')),
         alpha      = factor(alpha  > 0,  labels = c('dens indep', 'dens dep'))) %>%
  group_by(trial, n.pop0, bottleneck, alpha, evolve) %>%
  mutate(extinct = max(gen) < 15) %>%
  ungroup() %>%
  ggplot(aes(x = gen, y = n)) +
  geom_line(
    aes(
      group    = interaction(trial, n.pop0, evolve, bottleneck),
      linetype = evolve,
      size     = evolve,
      colour   = extinct
    )
  ) +
  scale_linetype_manual(values = 2:1) +
  scale_size_manual(values = c(0.25, 0.75)) +
  scale_color_manual(values = c('black', 'red')) +
  scale_y_log10() +
  facet_wrap(alpha ~ paste(n.pop0, bottleneck), ncol = 4) +
  theme(legend.position = 'bottom')

### Combine all designs, bottleneck only

rbind(design1, design2, design3) %>%
  filter(design > 0) %>%
  mutate(n.pop0     = factor(n.pop0 > 20, labels = c('small', 'large')),
         alpha      = factor(alpha  > 0,  labels = c('dens indep', 'dens dep')),
         design     = paste0('design', design)) %>%
  group_by(trial, n.pop0, design, alpha, evolve) %>%
  mutate(extinct = max(gen) < 15) %>%
  ungroup() %>%
  ggplot(aes(x = gen, y = n)) +
  geom_line(
    aes(
      group    = interaction(trial, n.pop0, evolve, design),
      linetype = evolve,
      size     = evolve,
      colour   = extinct
    )
  ) +
  scale_linetype_manual(values = 2:1) +
  scale_size_manual(values = c(0.25, 0.75)) +
  scale_color_manual(values = c('black', 'red')) +
  scale_y_log10() +
  facet_wrap(n.pop0 ~ paste(alpha, design), nrow = 4) +
  theme(legend.position = 'bottom')

### Average population size

design.means = rbind(design1, design2, design3) %>%
  filter(design > 0) %>%
  group_by(trial, design, n.pop0, alpha, evolve) %>%
  arrange(desc(gen)) %>%
  distinct(trial, design, n.pop0, alpha, evolve, .keep_all = TRUE) %>%
  filter(gen < 15) %>%
  uncount(weight = 15 - gen) %>%
  mutate(n = 0, gen = gen + (1:(15-gen[1]))) %>%
  ungroup() %>% 
  rbind(rbind(design1, design2, design3) %>% filter(design > 0)) %>%
  group_by(n.pop0, design, n.pop0, alpha, evolve, gen) %>%
  summarise(nbar = mean(n),
            nvar = var(n),
            nobs = n()) %>%
  ungroup()
  
design.means %>%
  filter(nbar > 0) %>%
  mutate(n.pop0     = factor(n.pop0 > 20, labels = c('small', 'large')),
         alpha      = factor(alpha  > 0,  labels = c('dens indep', 'dens dep')),
         design     = paste0('design', design)) %>%
  ggplot(aes(x = gen, y = nbar)) +
  geom_line(
    aes(
      group    = interaction(n.pop0, evolve, design),
      linetype = evolve,
      size     = evolve,
    )
  ) +
  # geom_ribbon(
  #   aes(
  #     ymin = nbar - 2 * sqrt(nvar / nobs),
  #     ymax = nbar + 2 * sqrt(nvar / nobs),
  #     group = interaction(n.pop0, evolve, design)
  #   ),
  #   alpha = 0.1
  # ) +
  scale_linetype_manual(values = 2:1) +
  scale_size_manual(values = c(0.25, 0.75)) +
  scale_y_log10() +
  facet_wrap(~ paste(n.pop0, alpha, design, sep = ', '), nrow = 4) +
  theme(legend.position = 'bottom')
