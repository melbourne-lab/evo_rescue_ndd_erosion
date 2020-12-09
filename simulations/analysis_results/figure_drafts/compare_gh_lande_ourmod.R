# Here: looking at genotypic and population change in our model and comparing to
# predctions from Lande 1976/G&H 1995.
# Our model allows genetic variation to change over time. Our analysis is also
# complicated by non-random extinctions affecting the averages.
# The Lande quantitative genetics framework assumes a constant genetic variance
# (and a genotype/phenotype which geometrically approaches the optimum).
# In this script I take the rate of change of the phenotype in the first
# generation (before all extinctions) and estimate the expected phenotypes over
# time according to Lande 1976.
# SN - 8 Dec 2020

### Setup

# Load packages
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

# Clear namespace
rm(list = ls())

# Load in all simulation data
all.data = read.csv('simulations/outputs/final_results/alldata_combined.csv')

### Start generating data.

# Get phenotypes (conditioned on population being extant)

all.z = all.data %>%
  # d is phenotypic distance from optimum (sensu G&H)
  mutate(d = 2.75 - zbar) %>%
  group_by(n.pop0, low.var, alpha, gen) %>%
  # Get mean and variance of phenotypic distances
  summarise(zbar = mean(d),
            zvar = var(d),
            n = n()) %>%
  ungroup() #%>%
  # factorize levels (not run here)
  mutate(n0 = factor(n.pop0, labels = c("Small", "Large")),
         alpha = factor(alpha, labels = c("Density independent", "Density dependent")),
         low.var = factor(low.var, labels = c("High diversity", "Low diversity")))

# Run the above, but conditioned on surviving to the end of the experiment
ext.z = all.data %>%
  mutate(d = 2.75 - zbar) %>%
  group_by(n.pop0, low.var, alpha, gen, extinct) %>%
  summarise(zbar = mean(d),
            zvar = var(d),
            n = n()) %>%
  ungroup() #%>%
  mutate(n0 = factor(n.pop0, labels = c("Small", "Large")),
         alpha = factor(alpha, labels = c("Density independent", "Density dependent")),
         low.var = factor(low.var, labels = c("High diversity", "Low diversity")))

# Here: get the mean rate of phenotypic change (+ uncertainty) for each param combo
all.k0 = all.data %>%
  # Get only frst two generations AND density ind. populatons
  # (density independent populations serving as reference)
  filter(gen %in% 1:2, !alpha) %>%
  # Get d (pheno. distnace from optimum)
  mutate(d = 2.75 - zbar) %>%
  # In each trial, get the rate of phenotypic approach to opt.
  # (from G&H call the factor k s/t
  # d_{t+1} = k d_t)
  group_by(n.pop0, low.var, trial) %>%
  summarise(k = d[2] / d[1]) %>%
  # Get mean and standard error across trials
  group_by(n.pop0, low.var) %>%
  summarise(kbar = mean(k),
            kse  = sqrt(var(k) / n()))

# Two things to note here: 
# First: for unknown reasons, this model doesn't quite fit in with the
# expectations from Lande, i.e., Lande's expression for k is off in our model.
# I've looked at this several times and can't figure out why. So, instead of
# using his formula I simply use an empirircal estimate of k.
# Second: it's possible that $k$ is influenced by the rate of population
# decline, which is affected by genetic variance, population size, and perhaps
# NDD (this may cause the issue above).

### Get first time step's genotypic change

# Use the observed k values above to get estimated z_t (`zest` for z est.)
# also use the bounds of a ~95% uncertainty interval
# Merge this with observed phenotypic data for comparison.
all.dk = merge(x = all.z, y = all.k0, by = c('n.pop0', 'low.var')) %>%
  mutate(zest = 2.75 * kbar^(gen - 1),
         zmin = 2.75 * (kbar - 2 * kse)^(gen-1),
         zmax = 2.75 * (kbar + 2 * kse)^(gen-1)) %>%
  mutate(n0 = factor(n.pop0, labels = c("Small", "Large")),
         alpha = factor(alpha, labels = c("Density independent", "Density dependent")),
         low.var = factor(low.var, labels = c("High diversity", "Low diversity")))

# Plot 
#   - NDD in purple
#   - density independent in solid black
#   - dashed line is expectated geometric increase to pheno. optimum

all.dk %>%
  ggplot(
    aes(
      x = gen, 
      group = interaction(n.pop0, low.var, alpha)
    )
  ) +
  geom_line(
    aes(
      y = zbar,
      colour = alpha
    )
  ) +
  geom_ribbon(
    aes(
      ymin = zbar - 2 * sqrt(zvar / n),
      ymax = zbar + 2 * sqrt(zvar / n),
      fill = alpha
    ),
    alpha = 0.1
  ) +
  geom_line(
    aes(
      y = zest
    ),
    linetype = 5
  ) +
  geom_ribbon(
    aes(
      ymin = zmin,
      ymax = zmax
    ),
    alpha = 0.2
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  facet_wrap(n0 ~ low.var)

### Now do the same but conditioned on popultions going extinct/surviving

ext.dk = merge(x = ext.z, y = all.k0, by = c('n.pop0', 'low.var')) %>%
  mutate(zest = 2.75 * kbar^(gen - 1),
         zmin = 2.75 * (kbar - 2 * kse)^(gen-1),
         zmax = 2.75 * (kbar + 2 * kse)^(gen-1)) %>%
  mutate(n0 = factor(n.pop0, labels = c("Small", "Large")),
         alpha = factor(alpha, labels = c("Density independent", "Density dependent")),
         low.var = factor(low.var, labels = c("High diversity", "Low diversity")))

# Plot

ext.dk %>%
  mutate_at(vars(zest, zmin, zmax, zbar), list(~ 2.75 - .)) %>%
  ggplot(
    aes(
      x = gen, 
      group = interaction(n.pop0, low.var, alpha, extinct)
    )
  ) +
  geom_line(
    aes(
      y = zbar,
      colour = alpha,
      linetype = extinct
    )
  ) +
  geom_ribbon(
    aes(
      ymin = zbar - 2 * sqrt(zvar / n),
      ymax = zbar + 2 * sqrt(zvar / n),
      fill = alpha
    ),
    alpha = 0.1
  ) +
  geom_line(
    aes(
      y = zest
    ),
    linetype = 5
  ) +
  geom_ribbon(
    data = . %>% mutate(!extinct),
    aes(
      ymin = zmin,
      ymax = zmax
    ),
    alpha = 0.1
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  scale_linetype_manual(values = 1:2) +
  facet_wrap(n0 ~ low.var)

### Could we use this to get population size? (*eyes emoji*)
# G&H 1995 provide (in piecemeal form) an expression for log pop size over time
# Let's plot that and compare with observed population sizes
# Note that G&H's model is continuous valued and handles extinctions by an
# extincton threshold, i.e., this expression should account for extinctions. Our
# model likewise includes extinctions in the mean population size.

# Get population size
all.n = all.data %>%
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
  # Get means and variances
  group_by(gen, n.pop0, low.var, alpha) %>%
  summarise(nbar = mean(n),
            nvar = var(n),
            n.trials = n()) %>%
  ungroup()

# Get global initial variables used in the G&H expression
all.0 = all.data %>%
  # Density independent trials only
  filter(!alpha) %>%
  # Get intial conditions for each trial
  filter(gen %in% 1) %>%
  # Get mean phenotypic variance (genetic var. v + env. var.)
  # I checked in sims and this is pretty similar to the true variances
  group_by(n.pop0, low.var) %>%
  summarise(sig2 = mean(v) + 0.5) %>%
  # Selection pressure
  mutate(w2 = 1 / 0.14 / 2) %>%
  # Get `k` (the rate of phenotypic change)
  merge(y = all.k0, by = c('n.pop0', 'low.var'))

# Now, get the density independent predictions from G&H
did.preds = expand.grid(n.pop0 = c(20, 100),
            low.var = c(TRUE, FALSE),
            gen = 1:15) %>%
  merge(all.0, by = c('n.pop0', 'low.var')) %>%
  mutate(t = gen - 1) %>%
  mutate(nt = log(n.pop0) + 
              t * (log(2) + (1/2)*log(w2 / (w2 + sig2))) -
              (2.75)^2 / (2 * (w2 + sig2)) * ((1 - kbar^(2*t)) / (1 - kbar^2))) %>%
  arrange(n.pop0, low.var, gen)
  
merge(all.n, did.preds %>% select(n.pop0, low.var, gen, nt)) %>%
  mutate(Nt = exp(nt)) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = nbar,
      colour = factor(alpha > 0),
      group = factor(alpha > 0)
    )
  ) +
  geom_ribbon(
    aes(
      ymin = nbar - 2 * sqrt(nvar / n.trials),
      ymax = nbar + 2 * sqrt(nvar / n.trials),
      fill = factor(alpha > 0),
      group = factor(alpha > 0)
    ),
    alpha = 0.1
  ) +
  geom_line(
    aes(
      y = Nt
    ),
    linetype = 5
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  scale_y_log10() +
  facet_wrap(n.pop0 ~ low.var)


ggplot(did.preds) +
  geom_line(aes(x = t, y = nt, group = interaction(n.pop0, low.var)))
