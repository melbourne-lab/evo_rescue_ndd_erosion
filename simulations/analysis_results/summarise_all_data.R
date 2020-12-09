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

all.data = read.csv('simulations/outputs/final_results/alldata_combined.csv')

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
  facet_grid(n.pop0 ~ low.var)

# A summary of population size based on extinction status.

popsize.ext.summ = popsize.all %>%
  group_by(trial, low.var, n.pop0, alpha) %>%
  mutate(extinct = any(!n)) %>%
  group_by(gen, low.var, n.pop0, alpha, extinct) %>%
  summarise(nbar = mean(n),
            nvar = var(n),
            n.trial = n())

# think about whether zeros belong in the above plot


# QAD plot
popsize.ext.summ %>%
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
  facet_grid(n.pop0 ~ low.var)
