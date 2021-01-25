library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

### Clear namespace

rm(list = ls())

### Global parameters

trials = 4000

pars = data.frame(
  end.time = 15,
  n.loci = 25,
  w.max = 2,
  theta = 2.75,
  wfitn = sqrt(1 / 0.14 / 2),
  sig.e = sqrt(0.5),
  pos.p = 0.5,
  alpha = 0
)

### Read in all data

all.data = 
  rbind(
    read.csv('simulations/outputs/alldata_n100_a000_hivar.csv'),
    read.csv('simulations/outputs/alldata_n100_a000_lowvar.csv'),
    read.csv('simulations/outputs/alldata_n20_a000_hivar.csv'),
    read.csv('simulations/outputs/alldata_n20_a000_lowvar.csv'),
    read.csv('simulations/outputs/alldata_n100_a035_hivar.csv'),
    read.csv('simulations/outputs/alldata_n100_a035_lowvar.csv'),
    read.csv('simulations/outputs/alldata_n20_a035_hivar.csv'),
    read.csv('simulations/outputs/alldata_n20_a035_lowvar.csv')
  )

### Inspect all data

nrow(all.data)

### Summarise

### Summary 1: population size

# ths requires inserting zeros where needed 

popsize.all = all.data %>%
  select(gen, trial, low.var, n.pop0, alpha, n) %>%
  # Add zeros if needed
  rbind(expand.grid(gen = 1:pars$end.time[1],
                    trial = 1:trials,
                    low.var = c(TRUE, FALSE),
                    n.pop0 = c(20, 100),
                    alpha = c(0, 0.0035),
                    n = 0)) %>%
  group_by(gen, trial, low.var, n.pop0, alpha) %>%
  summarise(n = sum(n)) %>%
  # Remove empty sims (never run for some reason - all 0)
  group_by(trial, low.var, n.pop0, alpha) %>%
  filter(any(n > 0)) %>%
  ungroup()

head(popsize.all)

# A summary of all population sizes.

popsize.all.summ = popsize.all %>%
  group_by(gen, low.var, n.pop0, alpha) %>%
  summarise(nbar = mean(n),
            nvar = var(n),
            n.trial = n())

# Quick and dirty plot of above:
popsize.all.summ %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = nbar,
      group = interaction(n.pop0, low.var, alpha),
      colour = factor(alpha)
    )
  ) +
  geom_ribbon(
    aes(
      ymin = nbar - 2 * sqrt(nvar / n.trial),
      ymax = nbar + 2 * sqrt(nvar / n.trial),
      group = interaction(n.pop0, low.var, alpha),
      fill = factor(alpha)
    ),
    alpha = 0.1
  ) +
  scale_colour_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  scale_y_log10() +
  facet_grid( ~ n.pop0)

# Related: time to extinction

ext.time.summ = popsize.all %>%
  group_by(trial, low.var, n.pop0, alpha) %>%
  mutate(ext.gen = max(gen[n>0])) %>%
  ungroup() %>%
  distinct(trial, ext.gen, n.pop0, low.var, alpha) %>%
  filter(ext.gen < 15) %>%
  group_by(n.pop0, low.var, alpha) %>%
  summarise(mean.ex.time = mean(ext.gen),
            var.ex.time = var(ext.gen),
            n = n())

ext.time.summ %>%
  ggplot() +
  geom_point(
    aes(
      x = interaction(alpha, n.pop0, low.var),
      y = mean.ex.time,
      colour = factor(alpha)
    )
  ) +
  geom_segment(
    aes(
      x = interaction(alpha, n.pop0, low.var),
      xend = interaction(alpha, n.pop0, low.var),
      y = mean.ex.time - 2 * sqrt(var.ex.time / n),
      yend = mean.ex.time + 2 * sqrt(var.ex.time / n)
    )
  )

### A statistical test

ext.times = popsize.all %>%
  group_by(trial, low.var, n.pop0, alpha) %>%
  mutate(ext.gen = max(gen[n>0])) %>%
  ungroup() %>%
  distinct(trial, ext.gen, n.pop0, low.var, alpha) %>%
  filter(ext.gen < 15) 

ext.times %>%
  aov(formula = ext.gen ~ factor(n.pop0) * low.var * factor(alpha)) %>%
  summary()

ext.times %>%
  aov(formula = ext.gen ~ factor(n.pop0) * low.var * factor(alpha)) %>%
  summary()

ext.times %>%
  aov(formula = ext.gen ~ factor(n.pop0) * low.var +
                          low.var * factor(alpha) +
                          factor(n.pop0) * factor(alpha)) %>%
  summary()

ext.times %>%
  aov(formula = ext.gen ~ low.var * factor(alpha) +
                          factor(n.pop0) * factor(alpha)) %>%
  summary()

ext.times %>%
  aov(formula = ext.gen ~ low.var * factor(alpha) +
                          factor(n.pop0) * factor(alpha)) %>%
  summary()

ext.times %>%
  lm(formula = ext.gen ~ factor(n.pop0) * low.var * factor(alpha)) %>%
  summary()

ext.times %>%
  lm(formula = log(ext.gen) ~ factor(n.pop0) * low.var * factor(alpha)) %>%
  summary()

ext.times %>%
  lm(formula = log(ext.gen) ~ factor(n.pop0) * low.var +
                              low.var * factor(alpha) +
                              factor(n.pop0) * factor(alpha)) %>%
  summary()

ext.times %>%
  lm(formula = log(ext.gen) ~ low.var * factor(alpha) +
                              factor(n.pop0) * factor(alpha)) %>%
  summary()

ext.times %>%
  lm(formula = log(ext.gen) ~ low.var * factor(alpha) +
                              factor(n.pop0) * factor(alpha)) %>%
  summary()

ext.times %>%
  lm(formula = log(ext.gen) ~ low.var + factor(n.pop0) * factor(alpha)) %>%
  summary()

ext.times %>%
  ggplot() +
  geom_boxplot(
    aes(
      x = interaction(alpha, n.pop0, low.var),
      y = ext.gen,
      colour = factor(alpha)
    ),
    notch = TRUE
  ) +
  scale_color_manual(
    values = c('black', 'purple'),
    labels = c("Dens. Ind.", "Dens. Dep.")
  ) +
  scale_x_discrete(labels = c(
    "Small, High\ndiversity", "", 
    "Large, High\ndiversity", "", 
    "Small, Low\ndiversity", "", 
    "Large, Low\ndiversity", "")) +
  labs(
    x = '', 
    y = 'Time to extinction',
    colour = ''
  ) +
  theme(
    axis.text.x = element_text(hjust = 0),
    legend.position = 'bottom'
  )


  scale_x_discrete(labels = c(
    "a", "b", "c", "d", "e", "f", "g", "h"
  ))

### Summary 2: 

######## Other junk (plots etc.)

# A summary of population size based on extinction status.

popsize.ext.summ = popsize.all %>%
  group_by(trial, low.var, n.pop0, alpha) %>%
  mutate(extinct = any(!n)) %>%
  group_by(gen, low.var, n.pop0, alpha, extinct) %>%
  summarise(nbar = mean(n),
            nvar = var(n),
            n.trial = n())

# think about whether zeros belong in the above plot

popsize.ext.summ.no.zero = popsize.all %>%
  group_by(trial, low.var, n.pop0, alpha) %>%
  mutate(extinct = any(!n)) %>%
  ungroup() %>%
  filter(n > 0) %>%
  group_by(gen, low.var, n.pop0, alpha, extinct) %>%
  summarise(nbar = mean(n),
            nvar = var(n),
            n.trial = n())

# QAD plot
popsize.ext.summ.no.zero %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = nbar,
      group = interaction(n.pop0, low.var, alpha, extinct),
      linetype = extinct,
      colour = factor(alpha)
    )
  ) +
  geom_ribbon(
    aes(
      ymin = nbar - 2 * sqrt(nvar / n.trial),
      ymax = nbar + 2 * sqrt(nvar / n.trial),
      group = interaction(n.pop0, low.var, alpha, extinct),
      linetype = extinct,
      fill = factor(alpha)
    ),
    alpha = 0.1
  ) +
  scale_colour_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  scale_y_log10() +
  facet_grid(n.pop0 ~ low.var) +
  ggsave(file = 'simulations/analysis_results/figure_drafts/test_figs/popsize_by_ext.pdf')

# Summary of population size based on generation of extinction

popsize.gen.summ.no.zero = popsize.all %>%
  group_by(trial, low.var, n.pop0, alpha) %>%
  mutate(ext.gen = max(gen[n>0])) %>%
  ungroup() %>%
  filter(n > 0) %>%
  group_by(gen, low.var, n.pop0, alpha, ext.gen) %>%
  summarise(nbar = mean(n),
            nvar = var(n),
            n.trial = n())

# QAD plot looking at population size by extinction
popsize.gen.summ %>%
  filter(nbar > 0) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    data = . %>% filter(ext.gen < 15),
    aes(
      y = nbar,
      group = interaction(n.pop0, low.var, alpha, ext.gen)
    ),
    colour = 'red'
  ) +
  geom_line(
    data = . %>% filter(ext.gen == 15),
    aes(
      y = nbar,
      group = interaction(n.pop0, low.var, alpha)
    ),
    colour = 'black'
  ) +
  scale_y_log10() +
  facet_wrap(alpha ~ interaction(n.pop0, low.var), ncol = 4) + 
  ggsave(file = 'simulations/analysis_results/figure_drafts/test_figs/popsize_by_ext_gen.pdf')

