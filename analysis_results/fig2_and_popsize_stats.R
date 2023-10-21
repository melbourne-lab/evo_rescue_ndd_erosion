# Main text Figure 2: Population size over time
# and related summary statistics

# Load packages
library(ggplot2)
library(dplyr)
library(tidyr)

# Clear namespace
rm(list = ls())

### Get necessary info for analytical solutions

# Get observed initial values of k for each treatment (from raw data)

# Read in data
all.data = read.csv('simulations/outputs/alldata_combined.csv')

# k values
kvals = all.data %>%
  # Get generations 1 and 2 (this is raw data before re-indexing generation)
  filter(gen < 3) %>%
  # Re-scale genotype to be distance from optimum (theta - zbar)
  mutate(zbar = 2.8 - zbar) %>%
  # For each trial in each treatment group:
  group_by(n.pop0, low.var, alpha, trial) %>%
  # summarise by getting the mean change in genotype, k, and variance in time step 1
  summarise(k = exp(diff(log(zbar))), v = mean(v[1])) %>%
  # Get the mean for each of the variance treatment groups
  # (more precise, and, is ~same across DI and size treatments)
  group_by(low.var) %>%
  summarise(kbar = mean(k), vbar = mean(v))

# Create a skeletal object with population size over time
analytics = expand.grid(
  n.pop0  = c(20, 100),
  low.var = c(TRUE, FALSE),
  alpha   = c(0, 0.0035),
  gen     = 0:15
) %>%
  merge(kvals) %>%
  mutate(
    w.max = 2,
    theta = 2.8,
    wfitn = sqrt(3.5),
    sig.e = sqrt(0.5)
  ) %>%
  mutate(
    zt = theta * kbar^gen,
    Nt = ifelse(!gen, n.pop0, NA),
  ) %>%
  arrange(n.pop0, low.var, alpha, gen)

# Use iterative procedure to get population size over time
for (i in 1:nrow(analytics)) {
  if (analytics$gen[i] > 0) {
    analytics$Nt[i] = with(analytics[i-1,], 
                           Nt * 
                             w.max * sqrt(wfitn^2 / (wfitn^2 + vbar + sig.e^2)) *
                             exp(-(zt^2) / (2 * (wfitn^2 + vbar + sig.e^2))) *
                             exp(-alpha * Nt)
    )
  }
}

# Plot of expected population size over time  
analytics %>%  
  ggplot(aes(x = gen, y = Nt)) +
  geom_line(aes(group = interaction(low.var, alpha), colour = alpha, linetype = low.var)) +
  facet_wrap(~ n.pop0) +
  scale_y_log10()

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
  )

# Add simulation and analytical estimates
alln = rbind(
  all.n %>% mutate(var.const = 0),
  analytics %>%
    rename(nbar = Nt) %>%
    select(gen, n.pop0, low.var, alpha, nbar) %>%
    mutate(nvar = 0, n.trials = 1, var.const = 1)
)

alln %>%  
  mutate(alpha = factor(alpha)) %>%
  ggplot(aes(x = gen, y = nbar)) +
  geom_line(aes(group = interaction(low.var, alpha, var.const), 
                colour = alpha, linetype = low.var, linewidth = var.const)) +
  scale_linewidth(range = c(1, 1/3), breaks = c(0, 1)) +
  scale_color_manual(values = c('black', 'purple')) +
  facet_wrap(~ n.pop0) +
  scale_y_log10()
# First visualization

alln %>%  
  mutate(alpha = factor(alpha)) %>%
  mutate(
    diver = ifelse(low.var, 'Low diversity', 'High diversity'),
    size0 = ifelse(n.pop0 > 20, 'Large', 'Small')
  ) %>%
  ggplot(aes(x = gen, y = nbar)) +
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
      group = interaction(alpha, var.const), 
      colour = alpha, 
      linewidth = var.const
    )
  ) +
  scale_color_manual(
    values = c('black', 'purple'),
    # labels = c("Density\nindependent", "Density\ndependent"),
    labels = c("Density independent", "Density dependent"),
    name = ""
  ) +
  geom_ribbon(
    aes(
      ymin = nbar - 2 * sqrt(nvar / n.trials),
      ymax = nbar + 2 * sqrt(nvar / n.trials),
      group = interaction(alpha, n.pop0, low.var, var.const),
      fill = alpha
    ),
    alpha = 0.2,
  ) +
  scale_fill_manual(
    values = c('black', 'purple'),
    # labels = c("Density\nindependent", "Density\ndependent"),
    labels = c("Density independent", "Density dependent"),
    name = ""
  ) +
  scale_linewidth(
    range = c(1, 1/3), 
    breaks = c(0, 1), 
    labels = c('Simulation', 'Analytical'), 
    name = ""
  ) +
  labs(x = 'Generation', y = 'Mean population size') +
  scale_y_log10() +
  # facet_wrap(~ paste0(size0, ', ', diver)) +
  facet_wrap(~ paste0(size0, ', ', diver), nrow = 1) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    legend.background = element_rect(fill = NA),
    # legend.direction = 'horizontal',
    # legend.position = c(0.8, 0.9),
    # legend.direction = 'vertical',
    # legend.position = c(0.9, 0.7),
    # legend.text = element_text(size = 8),
    legend.position = 'top',
    strip.background = element_rect(colour = 'black'),
    strip.text = element_text(size = 12)
  )

ggsave('analysis_results/figures/fig_pop_size_wide.png',
       width = 8, height = 3)


# ggsave('analysis_results/figures/fig_pop_size.png',
#        width = 8, height = 8)

all.data %>%
  mutate(gen = gen - 1) %>%
  group_by(gen, n.pop0, low.var, alpha) %>%
  filter(n() == 4000) %>%
  summarise(
    zbar = mean(2.8 - zbar),
    zvar = var(2.8 - zbar)
  ) %>%
  ggplot(aes(x = gen, y = zbar)) +
  geom_line(aes(group = alpha, colour = alpha)) +
  geom_line(
    data = analytics %>% filter(gen < 5),
    aes(x = gen, y = zt, group = alpha, colour = alpha),
    linetype = 2
  ) +
  facet_wrap(low.var ~ n.pop0)
