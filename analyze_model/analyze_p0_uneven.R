# Quick and dirty sim to satsfy my curiosity
#
# Look at dynamics when the populatio is maladapted AND starting at a point with
# lower genetic diversity due to the initial allele frequency (not fxation).
#
# Doing this by having the optimum phenotype be 0 (p = 0.5) but initializing the
# population at p = 0.25. With p approx 0.25, the genetic variance should be
# lower. But, genetic variance may increase over time as the population
# approches allele frequencies of 0.5.

library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

n.trials = 150

liszt = vector('list', n.trials)

set.seed(7756758)

for (trial in 1:n.trials) {
  liszt[[trial]] = sim(
    params = data.frame(
      end.time = 10,
      init.row = 1e3,
      n.loci = 25, 
      n.pop0 = 20,
      w.max = 2, 
      theta = 0, 
      wfitn = wfitn,
      sig.e = 0,
      pos.p = 0.25
    )
  )
  print(trial)
}

unrolled = unroller(liszt)

head(unrolled)

# Get mean phenotype (z_i) for each trial in each generation, see if the rate changes.
trials.z = unrolled %>%
  group_by(trial, gen) %>%
  summarize(n = n(),
            z = mean(z_i))

head(trials.z)

# Plot of phenotypic change over time
ggplot(trials.z) +
  geom_line(aes(x = gen, y = z, group = trial),
            size = 1/2)
# Can't actually pick anything up here.

# Plot of population size over time
ggplot(trials.z) +
  geom_line(aes(x = gen, y = n, group = trial),
            size = 1/2) +
  scale_y_log10()
# Interesting.


# Get mean frequency of the positive allele
# (note this is frequency across all genotypes)
trials.p = unrolled %>%
  select(-c(g_i, z_i, w_i, r_i, fem)) %>%
  gather(key = loc.copy, value = val, -c(i, gen, trial)) %>%
  group_by(trial, gen) %>%
  summarize(phat = mean(val > 0))

ggplot(trials.p) +
  geom_line(aes(x = gen, y = phat, group = trial),
            size = 1/2)
# Still not all the way there! Interesting.


# Now, look at the genetic variance over time.
trials.v = unrolled %>%
  select(-c(g_i, z_i, w_i, r_i, fem)) %>%
  gather(key = loc.copy, value = val, -c(i, gen, trial)) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(trial, gen, locus) %>%
  summarize(p = mean(val > 0)) %>%
  group_by(trial, gen) %>%
  summarise(v = sum(2 * p * (1-p) / 25))

# Now, plot dis.
trials.v %>%
  group_by(gen) %>%
  mutate(v.bar = mean(v)) %>%
  ggplot() +
  geom_line(aes(x = gen, y = v, group = trial),
            size = 1/2) +
  geom_line(aes(x = gen, y = v.bar),
            size = 2, colour = 'red')

# There IS a loss of diveristy from fixations happening
# Then there is a levelling out.
# This could be due to p getting larger or due to extinctions.
# Hard to know (should run this without extinctions to find out.)

