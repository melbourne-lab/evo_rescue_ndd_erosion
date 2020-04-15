# Evaluating a script which iterates population dynamics for several genomes
# (i.e., several $m$ values)
# Script which generated this data: test_sim_n.R
# (this was run on melb2 on April 14)
# SN April 14 2020

library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

# Read in data

testm = read.csv('test_m.csv')

head(testm)
nrow(testm)

table(testm$m)
table(testm$gen)

# Make a ggplot object.

ggm  = ggplot(testm)
ggmt = ggplot(testm, aes(x = gen))

### First: plot population sizes

ggmt +
  geom_line(aes(y = n.bar, group = factor(m), colour = factor(m)),
            size = 1.5) +
  scale_color_brewer(palette = 'RdGy') +
  theme_bw()

ggmt +
  geom_ribbon(aes(ymin = n.bar + 2 * sqrt(n.var / n.trials),
                  ymax = n.bar - 2 * sqrt(n.var / n.trials),
                  fill = factor(m), group = m),
              alpha = 0.1) +
  geom_line(aes(y = n.bar, group = factor(m), colour = factor(m)),
            size = 1.5) +
  scale_color_brewer(palette = 'RdGy') +
  scale_fill_brewer(palette = 'RdGy') +
  scale_y_log10() +
  theme_bw()

# Large error bars, but, looks like for sufficiently large $m$ (>4) dynamics are
# the same.

### Plot extinction probabilities

ggmt +
  geom_ribbon(aes(ymin = 1 - (n.trials/250) + 2 * sqrt((n.trials/250)*(1- n.trials/250) / n.trials),
                  ymax = 1 - (n.trials/250) - 2 * sqrt((n.trials/250)*(1- n.trials/250) / n.trials),
                  fill = factor(m), group = m),
              alpha = 0.1) +
  geom_line(aes(y = 1 - (n.trials / 250), group = factor(m), colour = factor(m))) +
  scale_color_brewer(palette = 'RdGy') +
  scale_fill_brewer(palette = 'RdGy') +
  theme_bw()
  

### Now, plot genetic values

ggmt +
  geom_ribbon(aes(ymin = g.bar + 2 * sqrt(g.var / n.trials),
                  ymax = g.bar - 2 * sqrt(g.var / n.trials),
                  fill = factor(m), group = m),
              alpha = 0.1) +
  geom_line(aes(y = g.bar, group = factor(m), colour = factor(m)),
            size = 1.5) +
  scale_color_brewer(palette = 'RdGy') +
  scale_fill_brewer(palette = 'RdGy') +
  theme_bw()

# man I need to run more trials.
# Generally... evidence of increasing speed of adaptation?
# I.e., more loci means a slightly more rapid rate of genotypic increase.
# However, differences are quite small as $m$ grows.

### Plot mean fitnesses

ggmt +
  geom_ribbon(aes(ymin = w.bar + 2 * sqrt(w.var / n.trials),
                  ymax = w.bar - 2 * sqrt(w.var / n.trials),
                  fill = factor(m), group = m),
              alpha = 0.1) +
  geom_line(aes(y = w.bar, group = factor(m), colour = factor(m)),
            size = 1.5) +
  scale_color_brewer(palette = 'RdGy') +
  scale_fill_brewer(palette = 'RdGy') +
  theme_bw()

# Data is coarse but it looks like the rate at which these populations approach
# w.bar > 1 is about the same regardless of $m$, but there's a levelling off
# that happens for small $m$, likely as a population with few loci saturates its
# postive alleles (i.e., there's an upper limit to how positive the genotypic
# values can get).

### Plot allele frequencies.

ggmt +
  geom_ribbon(aes(ymin = p.bar + 2 * sqrt(p.bar*(1-p.bar) / n.trials),
                  ymax = p.bar - 2 * sqrt(p.bar*(1-p.bar) / n.trials),
                  fill = factor(m), group = m),
              alpha = 0.1) +
  geom_line(aes(y = p.bar, group = factor(m), colour = factor(m)),
            size = 1.5) +
  scale_color_brewer(palette = 'RdGy') +
  scale_fill_brewer(palette = 'RdGy') +
  theme_bw()

### Plot genetic variance

ggmt +
  geom_ribbon(aes(ymin = v.bar + 2 * sqrt(v.var / n.trials),
                  ymax = v.bar - 2 * sqrt(v.var / n.trials),
                  fill = factor(m), group = m),
              alpha = 0.1) +
  geom_line(aes(y = v.bar, group = factor(m), colour = factor(m)),
            size = 1.5) +
  scale_color_brewer(palette = 'RdGy') +
  scale_fill_brewer(palette = 'RdGy') +
  theme_bw()

# Synthesis:
# As shown below, a small number of loci (m) is very constraining for changing
# allele frequencies. For a fixed phenotypic optimum (theta), as the number of
# loci increases, the max genotypic number increases (in this model). Related to
# this, the mean genotype needed to reach this optimum decreases (i.e., slowly
# approaches p = 1/2) as the number of loci increases. Th consequences of this
# are two-fold: 1) change in allele freuqency over time is much slower (because
# it's approaching a lower optimum) and 2) there's less loss of genetic
# diversity, as allele frequencies are closer to parity and the probability of
# fixation occurring is smaller.
# An interesting pattern: for a very small number of loci, population growth is
# smaller than for a larger number of loci. But, curiously, the probability of
# extinction is smaller as $m$ grows. I don't have a great explanation for this;
# it could be because it's faster (easier) to reach the new phenotypic optimum;
# however, concurrently there's a larger loss in phenotypic diversity when the
# number of loci is small; this may inhibit population growth.

# Plot initial allele distribution (p = 0.5) as $m$ changes

rm(list = ls())

li = vector('list', 25)

for (m in 10 + 2*(1:25)) {
  li[[m]] = data.frame(g = ((0:(2*m)) - m) / sqrt(m) ,
                       p = dbinom(0:(2*m), p = 1/2, size = 2*m),
                       m = rep(m, 2*m + 1))
}

lim = do.call(li, what = rbind)

ggplot(lim) +
  geom_line(aes(x = g, y = p, group = m, colour = m)) +
  geom_segment(aes(x = 2.6, xend = 2.6, y = 0, yend = 0.2)) +
  scale_color_viridis_c() +
  theme(panel.background = element_blank())

ggplot(lim) +
  geom_line(aes(x = g, y = p, group = m)) +
  geom_segment(aes(x = 2.6, xend = 2.6, y = 0, yend = 0.2)) +
  theme(panel.background = element_blank(),
        legend.position = 'none') +
  facet_wrap(~ m)


# Plot allele distribution when the mean population genotypic value
# is theta (2.6)

litheta = vector('list', 25)

for (m in 10 + 2*(1:25)) {
  litheta[[m]] = data.frame(g = ((0:(2*m)) - m) / sqrt(m) ,
                            p = dbinom(0:(2*m), p = 1/2*(2.6/sqrt(m) + 1), size = 2*m),
                            m = rep(m, 2*m + 1))
}

limtheta = do.call(litheta, what = rbind)

ggplot(limtheta) +
  geom_line(aes(x = g, y = p, group = m, colour = m)) +
  geom_segment(aes(x = 2.6, xend = 2.6, y = 0, yend = 0.2)) +
  scale_color_viridis_c() +
  theme(panel.background = element_blank())

ggplot(limtheta) +
  geom_line(aes(x = g, y = p, group = m)) +
  geom_segment(aes(x = 2.6, xend = 2.6, y = 0, yend = 0.2)) +
  theme(panel.background = element_blank(),
        legend.position = 'none') +
  facet_wrap(~ m)
