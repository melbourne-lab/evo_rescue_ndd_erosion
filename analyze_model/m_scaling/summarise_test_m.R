# Evaluating a script which iterates population dynamics for several genomes
# (i.e., several $m$ values)
# Script which generated this data: test_sim_n.R
# (this was run on melb2 on April 14)
# SN April 14 2020

library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

setwd('analyze_model/m_scaling/')

##### New data (1000 trials per m)

# Read in data

testm = read.csv('test_m_k.csv')

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

testm %>%
  filter(m > 1) %>%
  filter(gen > 8) %>%
  ggplot() +
  geom_ribbon(aes(x = gen,
                  ymin = n.bar + 2 * sqrt(n.var / n.trials),
                  ymax = n.bar - 2 * sqrt(n.var / n.trials),
                  fill = factor(m), group = m),
              alpha = 0.2) +
  geom_line(aes(x = gen, y = n.bar, group = factor(m), colour = factor(m)),
            size = 1) +
  scale_color_brewer(palette = 'RdGy') +
  scale_fill_brewer(palette = 'RdGy') +
  scale_y_log10() +
  theme_bw()

# Hmm... still can't distinguish much between any m > 1

### Plot extinction probabilities

ggmt +
  geom_ribbon(aes(ymin = 1 - (n.trials/1000) + 2 * sqrt((n.trials/1000)*(1- n.trials/1000) / n.trials),
                  ymax = 1 - (n.trials/1000) - 2 * sqrt((n.trials/1000)*(1- n.trials/1000) / n.trials),
                  fill = factor(m), group = m),
              alpha = 0.1) +
  geom_line(aes(y = 1 - (n.trials / 1000), group = factor(m), colour = factor(m))) +
  scale_color_brewer(palette = 'RdGy') +
  scale_fill_brewer(palette = 'RdGy') +
  theme_bw()

ggm +
  geom_line(aes(x = m, y = n.trials, group = gen, colour = gen))
  
# Looking at the data, it looks like P(extant trials) is high for low and high m
# but lower for intermediate m.

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

ggm +
  geom_ribbon(aes(x = m,
                  ymin = g.bar + 2 * sqrt(g.var / n.trials),
                  ymax = g.bar - 2 * sqrt(g.var / n.trials),
                  group = gen, fill = gen),
              alpha = 0.2) +
  geom_line(aes(x = m, y = g.bar, group = gen, colour = gen)) +
  scale_color_gradient(low = 'black', high = 'lightblue') +
  scale_fill_gradient(low = 'black', high = 'lightblue') +
  theme_bw()

# Again, rate of adaptation is pretty similar past having one allele.

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

# Beyond m = 1, change in fitness is about the same (as predicted above)

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

# This makes sense. Allele frequencies change more quickly when there is less
# room for adaptation. It asymptotically approaches... something.
# It doesn't look like it's approaching staying at p = 0.5 though.
# So eve when there is a massive number of loci the allele frequencies still change.

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

# Similarly to above, stll a loss in genetic variance.
# Notably the loss in variance is still somewhat large even for large m.

### Plot fixation probabilities
# note these measures are censored; as populations go extinct they're no longer included in the summation.

ggmt +
  geom_line(aes(y = p.fixp, group = factor(m), colour = factor(m)),
            size = 1.5, linetype = 2) +
  geom_line(aes(y = p.fixn, group = factor(m), colour = factor(m)),
            size = 1.5, linetype = 5) +
  scale_color_brewer(palette = 'RdGy') +
  theme_bw()

# Damn well this is interesting.
# First, there is fixation happening in all trials, no matter how many loci.
# More positive fixation occurs for smaller number of loci, likely due to the
# changes in allele frequencies.
# Note that over time, the positive curve flattens with larger m, likely because 
# the allele frequencies are closer to 1/2.
# There are negative fixations occurring too. But, the opposite pattern is
# opposite in m here - as m increases, the curve gets higher. What is happening
# here? The decrease is due to extinction; once a population goes extinct, its
# fixations are no longer being included in this count. So, it looks like there
# is selection happening, where populations with fixations are more likely to go
# extinct if they have fewer loci. This is worth pondering in more detail.

### Synthesis:
# For reasonably large $m$ (m > 4 to be safe but even m > 1) the population
# dynamics, rate of genotypic change, and rate of change of fitness are the
# same. So, the genotypic change (g) is most important for determining
# population dynamics. Similarly, extinction probabilities are similar across
# locus-numbers.
# There are differences in the change n allele frquencies and, therefore, the
# change in genotypic variance. For small $m$, there are fewer steps required to
# reach a fixed $\theta$, so the change in allele frequencies gets larger. As
# $m$ increases, though, there are more steps required to reach the optimum;
# this means that the rate of allele frequncy change is smaller s $m$ increases.
# However, despite the fact that there is less change in allele frequencies as
# $m$ increases (and the new optimal allele frequency is close to 0.5), the
# bottleneck is stll evident in change in genetic variance. No matter the allele
# frequency, during the descending phase, there is still (presumably) a loss of
# genetic variance as fixation occurs.
# Interestingly, the population dynamics and the allelic-dynamics appear to
# operate on different scales of m. Population dynamics (and genotypic change)
# very quickly converge (m > 1), but the allele frequencies and fixation
# probabilities, which control the genetic variance, converge on a more slowly
# (m > 9). This suggests to me the genotypic change is somewhat decoupled from
# the genetic variance.

ggm +
  geom_line(aes(x = g.bar, y = v.bar, group = m, colour = factor(m))) +
  scale_color_brewer(palette = 'RdGy')

ggm +
  geom_line(aes(x = p.bar, y = v.bar, group = m, colour = factor(m))) +
  scale_color_brewer(palette = 'RdGy')

ggm +
  geom_point(aes(x = p.bar, y = g.bar, colour = factor(m))) +
  geom_line(aes(x = p.bar, y = g.bar, group = m, colour = factor(m))) +
  scale_color_brewer(palette = 'RdGy')

lm(g.bar ~ p.bar * factor(m), testm) %>% summary()

# I feel like ths is informative... 2*sqrt(m) is slope btwn pbar and gbar
# So as m increases, smaller changes in allele frequency are needed to shift the mean genotype...
# I feel like this is just the relationship between p and g
# Yes - that's exactly what this is - it's sqrt(m) ( 2p - 1 )
# At least this confirms what we thought before about geno change.
# Now it's just a matter of knowing how allele frequencies change.

# This also seems to explain the decoupling of pheno and geno changes with
# changes in m; this relationship scales with sqrt(m), so you'll see thse two
# scale slightly differently with m.

##### Old data (250 trials per m)

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
  geom_segment(aes(x = 2.6, xend = 2.6, y = 0, yend = max(p))) +
  scale_color_viridis_c() +
  theme(panel.background = element_blank())

ggplot(lim) +
  geom_line(aes(x = g, y = p, group = m)) +
  geom_segment(aes(x = 2.6, xend = 2.6, y = 0, yend = max(p))) +
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
  geom_segment(aes(x = 2.6, xend = 2.6, y = 0, yend = max(p))) +
  scale_color_viridis_c() +
  theme(panel.background = element_blank())

ggplot(limtheta) +
  geom_line(aes(x = g, y = p, group = m)) +
  geom_segment(aes(x = 2.6, xend = 2.6, y = 0, yend = max(p))) +
  theme(panel.background = element_blank(),
        legend.position = 'none') +
  facet_wrap(~ m)
