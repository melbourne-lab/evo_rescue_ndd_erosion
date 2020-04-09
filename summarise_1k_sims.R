# Plot summaries of 1000 simulation runs
# Run on server, 8 April 2020, SN

# Parameters:
# end.time = 15,
# init.row = 1e4,
# n.loci = 25, 
# n.pop0 = 20,
# w.max = 2, 
# theta = 2.6, 
# wfitn = sqrt(1 / 0.14 / 2),
# sig.e = 0
# Note: model does not yet feature density dependence (alpha = 0)

library(ggplot2)
library(dplyr)
library(tidyr)

# Model output for 1000 simulations of the above
# All parameters are different, except n.loci = 25
k.sim.summary = read.csv('eg_1k_n0_20.csv')

k.sim.summary

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
  geom_line(aes(y = g.bar)) +
  geom_ribbon(aes(ymin = g.bar - 2 * sqrt(g.var / n.trials),
                  ymax = g.bar + 2 * sqrt(g.var / n.trials)),
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
  geom_ribbon(aes(ymin = p.fixp - 2 * sqrt(25 * p.fixp * (1 - p.fixp) / n.trials),
                  ymax = p.fixp + 2 * sqrt(25 * p.fixp * (1 - p.fixp) / n.trials)),
              alpha = 0.2, fill = 'blue') +
  geom_line(aes(y = p.fixn), colour = 'red') +
  geom_ribbon(aes(ymin = p.fixn - 2 * sqrt(25 * p.fixn * (1 - p.fixn) / n.trials),
                  ymax = p.fixn + 2 * sqrt(25 * p.fixn * (1 - p.fixn) / n.trials)),
              alpha = 0.2, fill = 'red') +
  labs(x = 'Time', y = 'Proportion of loci at fixation') +
  ggtitle('Fixation probability')
# Standard errors are negative...

# Plot allele frequencies
k.sim.summary %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = p.bar))

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
  scale_y_continuous(limits = c(0, 0.5)) +
  labs(x = 'Time', y = 'Genetic variance') +
  ggtitle('Genetic variance')
