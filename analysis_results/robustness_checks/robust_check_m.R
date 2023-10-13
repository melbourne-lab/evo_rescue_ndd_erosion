# Script to generate robustness plots for number of loci (m)

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

rm(list = ls())

# Read in file
sim.out = read.csv('simulations/outputs/robust_checks_out/nloci_robust_check.csv')

# Check to see it looks correct
head(sim.out)

# Adding rows which treat extinctions as zeros
sim.out = sim.out %>%
  mutate(gen = paste0('gen', gen)) %>%
  pivot_wider(names_from = gen, values_from = n, values_fill = 0) %>%
  pivot_longer(cols = starts_with('gen'), names_to = 'gen', values_to = 'n') %>%
  mutate(gen = gsub('gen', '', gen) %>% as.numeric()) %>%
  arrange(trial, gen)
  
# Checking to see if any zeros were added inappropriately (size-censoring)
sim.out %>%
  group_by(trial, n.loci) %>%
  filter(any(!n)) %>%
  mutate(nd = c(NA, diff(n))) %>%
  filter(nd < -1000)
# great - nothing cut out prematurely!

# How many rows (should be a neat even number)
nrow(sim.out)

# Get a summary data frame
# outputs: mean pop size, variance in pop size, survival/extinciton rates, sample size
sim.sum = sim.out %>%
  group_by(gen, n.pop0, lowvar, alpha, n.loci) %>%
  summarise(
    nbar = mean(n), nvar = var(n), 
    pext = mean(!n), nsurv = sum(n > 0), 
    nn = n()
  ) %>%
  ungroup() %>%
  mutate(
    n0 = factor(n.pop0, labels = c("Small", "Large")),
    alpha = factor(alpha, labels = c("Density independent", "Density dependent")),
    lowvar = factor(lowvar, labels = c("High diversity", "Low diversity"))
  )

head(sim.sum)

# A basic plot of population size
sim.sum %>%
  mutate(alpha = factor(alpha)) %>%
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
      group = alpha,
      colour = alpha
    )
  ) +
  labs(x = 'Generation', y = 'Mean population size') +
  facet_grid(paste0('m = ', n.loci) ~ paste0(n0, ', ', lowvar)) +
  scale_y_log10() +
  scale_colour_manual(values = c('black', 'purple'), '') +
  scale_fill_manual(values = c('black', 'purple'), '') +
  theme(
    axis.line = element_line(),
    legend.position = 'top',
    panel.background = element_blank(),
    strip.background.y = element_blank()
  )

ggsave(
  'analysis_results/robustness_checks/figs/fig_robust_nloci_popsize.png',
  width = 8, height = 5
)

# A plot of extinction rates over time
sim.sum %>%
  # Get proper indexing for extinctions
  filter(gen > 0) %>% mutate(gen = gen - 1) %>%
  # Labels for arranging plot
  mutate(alpha = factor(alpha)) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = pext,
      group = alpha,
      colour = alpha
    ),
    linewidth = 0.1
  ) +
  geom_ribbon(
    aes(ymin = 0, ymax = pext, group = alpha, fill = alpha),
    alpha = 0.2
  ) +
  geom_line(
    data = sim.sum %>%
      arrange(desc(gen)) %>% 
      group_by(n.pop0, lowvar, alpha, n.loci) %>% 
      mutate(pins.ext = c(NA, diff(nsurv)) / nsurv) %>% 
      filter(nsurv > 10, !is.na(pins.ext)) %>%
      mutate(alpha = factor(alpha)),
    aes(y = pins.ext, group = alpha, colour = alpha)
  ) +
  labs(x = 'Generation', y = 'Probability of extinction') +
  facet_grid(paste0('m = ', n.loci) ~ paste0(n0, ', ', lowvar)) +
  scale_colour_manual(
    values = c('purple', 'black'), 
    breaks = c("Density dependent", "Density independent"),
    'Instantaneous: '
  ) +
  scale_fill_manual(
    values = c('purple' ,'black'), 
    breaks = c("Density dependent", "Density independent"),
    '  Cumulative: '
  ) +
  guides(
    colour = guide_legend(direction = 'vertical', title.position = 'left'), 
    fill = guide_legend(direction = 'vertical', title.position = 'left')
  ) +
  theme(
    legend.position = 'top',
    legend.title.align = 0.5,
    legend.title = element_text(size = rel(0.8)),
    axis.line = element_line(),
    panel.background = element_blank(),
    strip.text.y = element_text(size = rel(1)),
    strip.background.y = element_blank()
  )

ggsave(
  'analysis_results/robustness_checks/figs/fig_robust_nloci_extinct.png',
  width = 8, height = 5
)
