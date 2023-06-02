# Re-create figure 5 (extinction vortex) and calculate some summary statistics

# Clear namespace
rm(list = ls())

# Load packages
library(cowplot)
library(ggplot2)
library(dplyr)
library(tidyr)

# Read in data
all.data = read.csv('simulations/outputs/alldata_combined.csv')

# Get data conditioned on extinction:
#   - r (log lambda) is log rate of population growth
#   - d is genotypic distance from the optimum
#   - k is proportion of genotypic load left per generation
#       (1 - k is rate of adaptation)
#   - varred is reduction in genetic variance
all.tau = all.data %>%
  # Reindex time to start at time = 0 (instead of time = 1)
  mutate(gen = gen - 1) %>%
  # Get data for each trial/population
  group_by(trial, n.pop0, low.var, alpha) %>%
  mutate(r = c(NA, diff(log(n))),
         d = abs(2.8 - gbar),
         k = c(NA, exp(diff(log(d)))),
         varred = c(NA, exp(diff(log(v))))) %>%
  # For each population in each time step, get tau, and extinction
  mutate(tau = max(gen) - gen + 1 - as.numeric(any(n < 2)),
         ext = max(gen) < 15 | any(n < 2)) %>%
  # Remove first generation (r = 0))
  filter(!is.na(r)) %>%
  # Remove generation of extinction (tau = 0) from df
  filter(tau > 0) %>%
  # For each "tau-step", get means and variances of 
  #   each variable of interest (r, k, nu)
  group_by(tau, ext, n.pop0, low.var, alpha) %>%
  summarise(rt = mean(r, na.rm = TRUE),
            kt = mean(k, na.rm = TRUE),
            nu = mean(varred, na.rm = TRUE),
            rtvar = var(r, na.rm = TRUE),
            ktvar = var(k, na.rm = TRUE),
            nuvar = var(varred, na.rm = TRUE),
            n  = n())

# Indices: tau = 0 is the first generation the population is extinct
# this *includes* a census with a population size of 1 
# (i.e. N = 1 corresponds to tau = 0).
# 
# A consequence of this: N = 1 does *not* get included in these estimated means
# (because we are filtering out tau = 0.)

### Individual panel plots

tau.r = all.tau %>%
  filter(ext) %>%
  ungroup() %>%
  mutate(n.pop0  = factor(n.pop0 > 20, labels = c('Small', 'Large')),
         low.var = factor(low.var, labels = c('High diversity', 'Low diversity'))) %>%
  ggplot(aes(x = tau, y = rt)) +
  geom_segment(
    aes(x = 1, xend = 14, y = 0, yend = 0),
    linetype = 3
  ) +
  geom_line(
    aes(
      colour = factor(alpha),
      linetype = ext
    ),
    linewidth = 1.25
  ) +
  geom_ribbon(
    aes(
      ymin = rt - 2 * sqrt(rtvar / n),
      ymax = rt + 2 * sqrt(rtvar / n),
      group = interaction(alpha, ext),
      fill = factor(alpha)
    ),
    alpha = 0.2
  ) +
  scale_x_reverse(breaks = (0:4)*3, labels = NULL) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  facet_wrap(~ paste(n.pop0, low.var, sep = ', '), ncol = 4) +
  labs(x = '', y = expression(atop("Population growth", " rate, " ~ group(langle, r[tau], rangle)))) +
  theme(legend.position = 'none',
        panel.grid.major = element_line(colour = 'gray88'),
        panel.background = element_rect(fill = 'white'),
        axis.ticks.x = element_blank(),
        plot.margin = margin(b = 0, r = 5, l = 5, unit = 'pt'),
        strip.background = element_rect(colour = 'black'))

tau.nu = all.tau %>%
  filter(ext) %>%
  ungroup() %>%
  mutate(n.pop0  = factor(n.pop0 > 20, labels = c('Small', 'Large')),
         low.var = factor(low.var, labels = c('High diversity', 'Low diversity'))) %>%
  ggplot(aes(x = tau, y = 1 - nu)) +
  geom_line(
    aes(
      colour = factor(alpha),
      linetype = ext
    ),
    linewidth = 1.25
  ) +
  geom_segment(
    aes(x = 1, xend = 14, y = 0, yend = 0),
    linetype = 3
  ) +
  geom_ribbon(
    aes(
      ymin = 1 - (nu - 2 * sqrt(nuvar / n)),
      ymax = 1 - (nu + 2 * sqrt(nuvar / n)),
      group = interaction(alpha, ext),
      fill = factor(alpha)
    ),
    alpha = 0.2
  ) +
  labs(x = '', y =  expression(atop("Loss of genetic", "variance, " ~ group(langle, nu[tau], rangle)))) +
  scale_x_reverse(breaks = (0:4)*3, labels = NULL) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  facet_wrap(~ paste(n.pop0, low.var, sep = ', '), ncol = 4) +
  theme(legend.position = 'none',
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.grid.major = element_line(colour = 'gray88'),
        panel.background = element_rect(fill = 'white'),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, b = 0, r = 5, l = 5, unit = 'pt'))

tau.k = all.tau %>%
  ungroup() %>%
  filter(ext) %>%
  mutate(n.pop0  = factor(n.pop0 > 20, labels = c('Small', 'Large')),
         low.var = factor(low.var, labels = c('High diversity', 'Low diversity'))) %>%
  ggplot(aes(x = tau, y = 1 - kt)) +
  geom_segment(
    aes(x = 1, xend = 14, y = 0, yend = 0),
    linetype = 3
  ) +
  geom_line(
    aes(
      colour = factor(alpha),
      linetype = ext
    ),
    linewidth = 1.25
  ) +
  geom_ribbon(
    aes(
      ymin = 1 - (kt - 2 * sqrt(ktvar / n)),
      ymax = 1 - (kt + 2 * sqrt(ktvar / n)),
      group = interaction(alpha, ext),
      fill = factor(alpha)
    ),
    alpha = 0.2
  ) +
  scale_x_reverse(breaks = (0:4)*3) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  facet_wrap(~ paste(n.pop0, low.var, sep = ', '), ncol = 4) +
  labs(x = expression('Time to extinction, ' ~ tau), 
       y = expression(atop('Rate of adaptation,', ~ 1-group(langle, k[tau], rangle)))) +
  theme(legend.position = 'none',
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.grid.major = element_line(colour = 'gray88'),
        panel.background = element_rect(fill = 'white'),
        plot.margin = margin(t = 0, r = 5, l = 5, unit = 'pt'))

plot_grid(
  tau.r,
  tau.nu,
  tau.k,
  labels = c('(A)', '(B)', '(C)'),
  label_x = c(-0.015, -0.015, -0.015),
  label_y = c(1, 1.12, 1.12),
  label_size = 12,
  nrow = 3
) %>%
  save_plot(filename = 'analysis_results/figures/fig_vortex.png',
            base_width = 8, base_height = 6)

### Summary statistics

# Proportion of genetic variation (nu) lost among extinct populations
# in the generation preceding extinction (tau = 1)
all.tau %>%
  filter(tau %in% 1, ext) %>%
  mutate(nu = 1 - nu) %>%
  select(n.pop0, low.var, alpha, nu)

# DMaximum difference in nu (within a timestep) for DI/NDD populations
all.tau %>%
  filter(ext) %>%
  mutate(nu = 1 - nu) %>%
  select(tau, n.pop0, low.var, alpha, nu) %>%
  mutate(alpha = ifelse(alpha > 0, 'ndd', 'di')) %>%
  pivot_wider(names_from = alpha, values_from = nu) %>%
  mutate(alpha.diff = ndd - di) %>%
  group_by(n.pop0, low.var) %>%
  slice_max(alpha.diff)
