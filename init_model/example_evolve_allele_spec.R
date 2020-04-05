# An example of allele-specific model.
# See script "model_building/sim_evolve_allele_spec_brs.R'
# Commit 7de5cf5; commit and example from 2020 Mar 18
# SN

##### Setup

# Load in packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

# Clear your namespace
rm(list = ls())

# Source code needed to run simulations.
source('model_building/sim_evolve_allele_spec_brs.R')

##### Run one instance of simulation

# Set seed for reproducibility
set.seed(40904)

# Run one simulation
test.sim = sim(a = c(-1/2, 1/2),
               params = data.frame(end.time = 15,      # Length of simulatin
                                   init.row = 1e4,     # Number of rows to initialize data
                                   n.loci = 25,        # Number of loci
                                   n.pop0 = 20,        # Initial population size
                                   w.max = 2,          # Max fitness of an individual
                                   theta = 2.6,        # Distance from phenotypic optimum
                                   w.sqd = 1 / 0.14))  # Selection strength

# First rows of simulation
head(test.sim)
# Columns:
#   - i = individual number
#   - b_i = "genetic value" akin to breeding value of individual
#   - w_i = fitness of individual i
#   - r_i = number of offspring of individual i (Poisson draw)
#   - fem = boolean for if individual is female
#   - gen = generation number of individual
#   - a1, b1, etc. = effect of locus on genotype
#       (note: a1 and b1 are two copies of locus 1, a2 and b2 are
#       two copies of locus 2, etc.)
nrow(test.sim)
# One row for every individual in the simualtin.

### Analyze this simulation

# Get summary statistics for this trial run
# n = population size
# b.bar, b.var = mean and variance in genetic vals.
# w.bar, w.var = mean and variance in fitness
test.sim.summary = test.sim %>%
    group_by(gen) %>%
    summarise(n = n(),
              b.bar = mean(b_i),
              b.var = var(b_i),
              w.bar = mean(w_i),
              w.var = var(w_i)
              )

# Plot ech of these.

# Form ggplot base.
test.gg.base = ggplot(test.sim.summary, aes(x = gen))

# Plot of population size over time
test.gg.base +
    geom_line(aes(y = n)) +
    labs(x = 'Time', y = 'Population size') +
    ggtitle('Population size')

# Plot of mean genetic value of individual over time.
test.gg.base +
    geom_ribbon(aes(ymin = b.bar - 2 * sqrt(b.var),
                    ymax = b.bar + 2 * sqrt(b.var)),
                alpha = 0.2) +
    geom_line(aes(y = b.bar)) +
    geom_segment(aes(x = 0, xend = 15, y = 2.6, yend = 2.6)) +
    geom_segment(aes(x = 0, y = 2.6 - sqrt(log(2) / 0.14),
                     xend = 15, yend = 2.6 - sqrt(log(2) / 0.14)),
                 linetype = 2) +
    labs(x = 'Time', y = 'Mean genetic value') +
    ggtitle('Genetic value')

# Plot of mean fitness over time
test.gg.base +
    geom_segment(aes(x = 0, xend = 15, y = 1, yend = 1)) +
    geom_ribbon(aes(ymin = w.bar - 2 * sqrt(w.var),
                    ymax = w.bar + 2 * sqrt(w.var)),
                alpha = 0.2) +
    geom_line(aes(y = w.bar)) +
    labs(x = 'Time', y = 'Mean fitness') +
    ggtitle('Fitness')

# One way that genetic variation is lost is fixation
# Here is a plot of frequencies of the _positive_ allele
# at each locus over time

test.sim.alleles = test.sim %>%
    gather(key = chr.locus, val = val, -c(i, b_i, w_i, r_i, fem, gen)) %>%
    mutate(locus = gsub('^[ab]', 'locus', chr.locus)) %>%
    group_by(gen, locus) %>%
    summarise(p.pos = mean(val > 0)) 

# Rows = time step
# Columns = time series of allele frequencies for each locus
test.sim.alleles %>%
    spread(key = locus, val = p.pos)

# Plot
# Purple line is average
test.sim.alleles %>%
    ggplot(aes(x = gen, y = p.pos)) +
    stat_summary(aes(y = p.pos), fun = 'mean', geom = 'line',
                 colour = 'purple', size = 3) +
    geom_line(aes(group = locus)) + 
    labs(x = 'Time', y = 'Frequency of "positive" allele') +
    ggtitle('Allele frequency')

# Fixation of positive, negative alleles over time
test.sim.alleles %>%
    group_by(gen) %>%
    summarise(pos.fixation = sum(p.pos == 1),
              neg.fixation = sum(p.pos == 0))
# There is some fixation of negative alleles;
# this is likely happening when population size is small
# Over time there's fixation of positive alleles as well

# Genetic variance

test.sim.gvar = test.sim %>%
    gather(key = chr.locus, val = val, -c(i, b_i, w_i, r_i, fem, gen)) %>%
    mutate(locus = gsub('^[ab]', '', chr.locus)) %>%
    group_by(gen, locus) %>%
    summarise(p.pos = mean(val > 0),
              g.var = 2 * p.pos * (1 - p.pos) / 25) %>% # note: if you run w/ a different number of alleles, change '25'
    group_by(gen) %>%
    summarise(g.var = sum(g.var))

# Plot
test.sim.gvar %>%
    ggplot(aes(x = gen, y = g.var)) +
    geom_line() +
    labs(x = 'Time', y = 'Genetic variance') + 
    scale_y_continuous(limits = c(0, 0.5)) +
    ggtitle('Genetic variance')

##### Look at ensemble of simulations

# Model output for 1000 simulations of the above
# All parameters are different, except n.loci = 40
k.sim.summary = read.csv('model_building/eg_1k_n0_20.csv')

k.sim.summary
# Columns:
#   gen = generation
#   n.trials = number of extant populations
#   v.bar, v.var = mean and variance of genetic variance
#   n.bar, n.var = mean and variance in population size
#   r.bar, r.var = mean and variance in # offspring
#       (note: r_i is Poisson draw, so this includes
#       process variance)
#   p.fixp, p.fixn = proportion of loci with positive/negative 
#       alleles fixed

# Plot of mean population size over time (including extinctions)
k.sim.summary %>%
    ggplot(aes(x = gen)) +
    geom_line(aes(y = n.bar)) +
    geom_ribbon(aes(ymin = n.bar - 2 * sqrt(n.var / n.trials),
                    ymax = n.bar + 2 * sqrt(n.var / n.trials)),
                alpha = 0.2) +
    scale_y_log10() +
    labs(x = 'Time', y = 'Population size') +
    ggtitle('Population size')

# Plot of mean breeding value over time
k.sim.summary %>%    
    ggplot(aes(x = gen)) +
    geom_line(aes(y = b.bar)) +
    geom_ribbon(aes(ymin = b.bar - 2 * sqrt(b.var / n.trials),
                    ymax = b.bar + 2 * sqrt(b.var / n.trials)),
                alpha = 0.2) +
    labs(x = 'Time', y = 'Mean genotpyic value') +
    ggtitle('Genetic value')


# Plot of mean fitness over time
k.sim.summary %>%    
    ggplot(aes(x = gen)) +
    geom_line(aes(y = w.bar)) +
    geom_ribbon(aes(ymin = w.bar - 2 * sqrt(w.var / n.trials),
                    ymax = w.bar + 2 * sqrt(w.var / n.trials)),
                alpha = 0.2) +
    labs(x = 'Time', y = 'Mean fitness') +
    ggtitle('Fitness')

# Plot of fixation probabilities
k.sim.summary %>%
    ggplot(aes(x = gen)) +
    geom_line(aes(y = p.fixp), colour = 'blue') +
    geom_ribbon(aes(ymin = p.fixp - 2 * sqrt(40 * p.fixp * (1 - p.fixp) / n.trials),
                    ymax = p.fixp + 2 * sqrt(40 * p.fixp * (1 - p.fixp) / n.trials)),
                alpha = 0.2, fill = 'blue') +
    geom_line(aes(y = p.fixn), colour = 'red') +
    geom_ribbon(aes(ymin = p.fixn - 2 * sqrt(40 * p.fixn * (1 - p.fixn) / n.trials),
                    ymax = p.fixn + 2 * sqrt(40 * p.fixn * (1 - p.fixn) / n.trials)),
                alpha = 0.2, fill = 'red') +
    labs(x = 'Time', y = 'Proportion of loci at fixation') +
    ggtitle('Fixation probability')
# Note: this starts with a non-zero amount of fixation on average... 
# not sure why this would be...

# Plot of mean genetic variance over time
k.sim.summary %>%
    ggplot(aes(x = gen)) +
    geom_line(aes(y = v.bar)) +
    geom_ribbon(aes(ymin = v.bar - 2 * sqrt(v.var / n.trials),
                    ymax = v.bar + 2 * sqrt(v.var / n.trials)),
                alpha = 0.2) +
    labs(x = 'Time', y = 'Genetic variance') +
    ggtitle('Genetic variance')

# Above plot but with scaled axes
k.sim.summary %>%
    ggplot(aes(x = gen)) +
    geom_line(aes(y = v.bar)) +
    geom_ribbon(aes(ymin = v.bar - 2 * sqrt(v.var / n.trials),
                    ymax = v.bar + 2 * sqrt(v.var / n.trials)),
                alpha = 0.2) +
    scale_y_continuous(limits = c(0, 0.25)) +
    labs(x = 'Time', y = 'Genetic variance') +
    ggtitle('Genetic variance')

### Genetic variance in the 1000 trials

# Read in script
all.sigma.a = read.csv('model_building/allelic_all_sigma_a.csv')

# Plot variance in the first 50 trials
all.sigma.a %>%
    filter(trial %in% 1:50) %>%
    group_by(trial) %>%
    mutate(extinct = max(gen) < 15) %>%
    ggplot(aes(x = gen)) +
    geom_line(aes(y = sigma.a, group = trial, colour = extinct),
              size = 0.4) +
    scale_color_manual(values = c('black', 'red'))

# Plot all trials
all.sigma.a %>%
    group_by(gen) %>%
    summarise(sig.bar = mean(sigma.a),
              sig.var = var(sigma.a),
              n = n()) %>%
    ggplot() +
    geom_line(aes(x = gen, y = sig.bar)) +
    geom_point(aes(x = gen, y = sig.bar,
                   size = n / max(n))) +
    geom_ribbon(aes(x = gen,
                    ymin = sig.bar - 2 * sqrt(sig.var / n),
                    ymax = sig.bar + 2 * sqrt(sig.var / n)),
                alpha = 0.2)
# (points are proportional to the number of trials which have not gone extinct)
