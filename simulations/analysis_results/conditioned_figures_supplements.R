# Supplemental figures for various figures of:
#   - state variables conditioned on extinction/survival at end of simulation
#   - state variables conditioned on survival/generation of extinction

library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

# Load in all simulation data
all.data = read.csv('simulations/outputs/final_results/alldata_combined.csv')

### Fig C1: population size

# Summarise data by param combo, generation, extinction status
# n.b. need to add zeros for extinct populations
ext.n = all.data %>%
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
  # Get rid of trials where population was erroneously added
  group_by(trial, n.pop0, low.var, alpha) %>%
  filter(any(n>0)) %>%
  # Add extinction status to trial
  group_by(trial, n.pop0, low.var, alpha) %>%
  mutate(ext = any(!n)) %>%
  ungroup() %>%
  # Aggregate (mean and variance in each generation)
  group_by(gen, n.pop0, low.var, alpha, ext) %>%
  summarise(nbar = mean(n),
            nvar = var(n),
            n.trials = n()) %>%
  ungroup() %>%
  mutate(n0 = factor(n.pop0, labels = c("Initially small", "Initially large")),
         alpha = factor(alpha, labels = c("Density independent", "Density dependent")),
         low.var = factor(low.var, labels = c("Low genetic diversity", "High genetic diversity"))) %>%
  # For purposes of plotting on log scale, remove records where nbar is 0
  filter(nbar > 0)

ext.n %>%
  ggplot(aes(x = gen, y = nbar)) +
  geom_line(
    aes(
      colour = factor(alpha),
      group = interaction(ext, alpha),
      linetype = ext
      )
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000),
                labels = c(0.1, 1, 10, 100, 1000)) +
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 4) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_line(colour = 'gray88'),
        panel.grid.minor = element_blank(),
        legend.position = 'none') #,
        # strip.text = element_blank(),
        # strip.background = element_blank())
