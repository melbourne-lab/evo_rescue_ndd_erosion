# Main text Figure 2: Population size over time
# and related summary statistics

# Clear namespace
rm(list = ls())

# Load packages
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

# Load in all simulation data
all.data = read.csv('simulations/outputs/alldata_combined.csv')

##### Figure 2

### Aggregate mean population size (and variance for standard error)
all.n = all.data %>%
  # Set generation time data to go 0-16
  mutate(gen = gen - 1) %>%
  # Get relevant columns only
  select(trial, gen, n, n.pop0, low.var, alpha) %>%
  # Add extinctions in (for mean population size)
  rbind(
    expand.grid(trial = 1:4000,
                gen = 1:15,
                n = 0,
                n.pop0 = c(20, 100),
                low.var = c(TRUE, FALSE),
                alpha = c(0, 0.0035))
  ) %>%
  group_by(trial, gen, n.pop0, low.var, alpha) %>%
  summarise(n = sum(n)) %>%
  # Get rid of trials where population was erroneously added (if any)
  group_by(trial, n.pop0, low.var, alpha) %>%
  filter(any(n > 0)) %>%
  # Get means and variances
  group_by(gen, n.pop0, low.var, alpha) %>%
  summarise(
    nbar = mean(n),
    nvar = var(n),
    n.trials = n()
  ) %>%
  ungroup() %>%
  mutate(n0 = factor(n.pop0, labels = c("Initially small", "Initially large")),
         alpha = factor(alpha, labels = c("Density independent", "Density dependent")),
         low.var = factor(low.var, labels = c("Low genetic diversity", "High genetic diversity")))

# Create figure

all.n %>%
  mutate(n0 = factor(n0, levels = levels(n0)[2:1])) %>%
  ggplot(aes(x = gen)) +
  geom_segment(
    aes(
      x = 0, xend = 15,
      y = n.pop0, yend = n.pop0
    ),
    linetype = 3,
    size = 0.5,
    colour = 'gray'
  ) +
  geom_line(
    aes(
      y = nbar,
      group = interaction(alpha, n.pop0, low.var),
      linetype = low.var,
      colour = alpha
    ),
    size = 1.25
  ) +
  geom_ribbon(
    aes(
      ymin = nbar - 2 * sqrt(nvar / n.trials),
      ymax = nbar + 2 * sqrt(nvar / n.trials),
      group = interaction(alpha, n.pop0, low.var),
      fill = alpha
    ),
    alpha = 0.2,
    size = 0.125
  ) +
  scale_color_manual(
    values = c('black', 'purple'),
    labels = c("Density\nindependent", "Density\ndependent"),
    name = ""
  ) +
  scale_fill_manual(
    values = c('black', 'purple'),
    labels = c("Density\nindependent", "Density\ndependent"),
    name = ""
  ) +
  scale_linetype(# _manual(
    # values = c(1, 5),
    name = "",
    labels = c("High diversity", "Low diversity")
  ) +
  labs(x = 'Generation', y = 'Mean population size') +
  scale_y_log10() +
  facet_wrap(~ n0) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    legend.background = element_rect(fill = NA),
    legend.direction = 'horizontal',
    legend.position = c(0.2, 0.85),
    legend.text = element_text(size = 12),
    strip.background = element_rect(colour = 'black'),
    strip.text = element_text(size = 12)
  ) #+

ggsave('analysis_results/figures/fig_pop_size.png',
       width = 8, height = 5)

### Summary statistics

# Minimum expected size, by treatment group:

all.n %>%
  group_by(n.pop0, low.var, alpha) %>%
  slice_min(order_by = nbar) %>%
  select(-c(nvar, n.trials)) %>%
  pivot_wider(names_from = alpha, values_from = c(gen, nbar)) %>%
  mutate(pop.pcts = `nbar_Density dependent` / `nbar_Density independent`)

# Final expected size, by treatment group:

all.n %>%
  filter(gen %in% 15) %>%
  group_by(n.pop0, low.var, alpha) %>%
  select(-c(nvar, n.trials, gen)) %>%
  pivot_wider(names_from = alpha, values_from = nbar) %>%
  mutate(pop.pcts = `Density dependent` / `Density independent`)
