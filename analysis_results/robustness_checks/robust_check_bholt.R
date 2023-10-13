# Script to generate robustness plots with Beverton-Holt density dependence

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

rm(list = ls())

# Read in file
sim.out = read.csv('simulations/outputs/robust_checks_out/bholt_robust_check.csv')

# Check to see it looks correct
head(sim.out)
str(sim.out)

# Adding rows which treat extinctions as zeros
sim.out = sim.out %>%
  mutate(gen = paste0('gen', gen)) %>%
  pivot_wider(names_from = gen, values_from = n, values_fill = 0) %>%
  pivot_longer(cols = starts_with('gen'), names_to = 'gen', values_to = 'n') %>%
  mutate(gen = gsub('gen', '', gen) %>% as.numeric()) %>%
  arrange(trial, gen)

# Checking to see if any zeros were added inappropriately (size-censoring)
sim.out %>%
  group_by(trial) %>%
  filter(any(!n)) %>%
  mutate(nd = c(NA, diff(n))) %>%
  filter(nd < -1000)
# one populations (small?) that get censored

# How many rows remaining (will not be a neat even number)
nrow(sim.out) 

# Get a summary data frame
# outputs: mean pop size, variance in pop size, survival/extinciton rates, sample size
sim.sum = sim.out %>%
  group_by(gen, n.pop0, lowvar, alpha, delta, beta) %>%
  summarise(
    nbar = mean(n), nvar = var(n), 
    pext = mean(!n), nsurv = sum(n > 0), 
    nn = n()
  ) %>%
  ungroup() %>%
  mutate(
    n0 = factor(n.pop0, labels = c("Small", "Large")),
    lowvar = factor(lowvar, labels = c("High diversity", "Low diversity")),
    ndd = factor(ifelse(!alpha, 0, delta)),
    ndd = factor(ndd, levels = levels(ndd)[c(1, 4:2)]),
    beta = round(beta, 4)
  )

head(sim.sum)

sim.sum %>%
  ggplot(aes(x = gen)) +
  geom_segment(
    aes(
      x = 0, xend = 15,
      y = n.pop0, yend = n.pop0
    ),
    linetype = 3,
    linewidth = 0.5,
    colour = 'gray'
  ) +
  geom_line(
    aes(
      y = nbar,
      group = ndd,
      colour = ndd
    )
  ) +
  labs(x = 'Generation', y = 'Mean population size') +
  facet_grid( ~ paste0(n0, ', ', lowvar)) +
  scale_y_log10() +
  scale_colour_manual(
    values = c("black", "#458B00", "#66CD00", "#7FFF00"), 
    labels = c(
      "Density independent",
      expression(paste("Density dependent, ", delta == 1)),
      expression(paste("Density dependent, ", delta == 0.75)),
      expression(paste("Density dependent, ", delta == 0.5))
    ),
    ''
  ) +
  theme(
    axis.line = element_line(),
    legend.position = 'top',
    panel.background = element_blank()
  )

ggsave(
  'analysis_results/robustness_checks/figs/fig_robust_bholt_popsize.png',
  width = 8, height = 5
)

# A plot of extinction rates over time
sim.sum %>%
  # Get proper indexing for extinctions
  filter(gen > 0) %>% mutate(gen = gen - 1) %>%
  # Labels for arranging plot
  mutate(
    alpha = factor(alpha), 
    facet.x = paste0(n0, ', ', lowvar),
    facet.y = paste("delta ==", delta)
  ) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = pext,
      group = ndd,
      colour = ndd
    ),
    linewidth = 0.1
  ) +
  geom_ribbon(
    aes(ymin = 0, ymax = pext, group =  ndd, fill = ndd),
    alpha = 0.2
  ) +
  geom_line(
    data = sim.sum %>%
      arrange(desc(gen)) %>% 
      group_by(n.pop0, lowvar, alpha, delta) %>% 
      mutate(pins.ext = c(NA, diff(nsurv)) / nsurv) %>% 
      filter(nsurv > 10, !is.na(pins.ext)) %>%
      mutate(facet.x = paste0(n0, ', ', lowvar)),
    aes(y = pins.ext, group = ndd, colour = ndd)
  ) +
  labs(x = 'Generation', y = 'Probability of extinction') +
  facet_wrap(~ facet.x, nrow = 1) +
  scale_colour_manual(
    values = c("black", "#458B00", "#66CD00", "#7FFF00"), 
    labels = c(
      "Density independent",
      expression(paste("Density dependent, ", delta == 1)),
      expression(paste("Density dependent, ", delta == 0.75)),
      expression(paste("Density dependent, ", delta == 0.5))
    ),
    'Instantaneous:  '
  ) +
  scale_fill_manual(
    values = c("black", "#458B00", "#66CD00", "#7FFF00"), 
    labels = c(
      "Density independent",
      expression(paste("Density dependent, ", delta == 1)),
      expression(paste("Density dependent, ", delta == 0.75)),
      expression(paste("Density dependent, ", delta == 0.5))
    ),
    'Cumulative: '
  ) +
  guides(
    colour = guide_legend(direction = 'vertical', title.position = 'left'), 
    fill = guide_legend(direction = 'vertical', title.position = 'left')
  ) +
  theme(
    axis.line = element_line(),
    legend.position = 'top',
    legend.title = element_text(size = rel(0.8)),
    legend.text = element_text(hjust = 0),
    legend.title.align = 0,
    panel.background = element_blank(),
  )

ggsave(
  'analysis_results/robustness_checks/figs/fig_robust_bholt_extinct.png',
  width = 8, height = 5
)
