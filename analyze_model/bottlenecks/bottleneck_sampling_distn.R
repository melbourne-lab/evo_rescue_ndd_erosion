# Look at the distribution of allelic diversity left when bottlenecking a
# population at two individuals (male and female, let's say). Sources w/in the
# NSF proposal say that there should be a 25% reduction in heterozygosity for 29
# generations (e.g., Allendorf 1986). Here I'll test loss of heterozygosity
# simply from sampling. Note that this is something that should be easily
# expressed analytically.

# SN June 9 2020

##### Load packages and associated scripts. Clear namespace.
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

rm(list = ls())

source('base_model_source/sim_functions.R')

##### Sampling at a trait complexity (i.e., n. loci)

# I'm doing hierarchical partioning of variance.
# Variance in outcomes (loss of heterozygosity and/or fixation) can be
# attributed to two sources: the two individuals sampled, and the population
# they are sampled from. Here I'll partition variance as follows:
# For a fixed number of alleles (25), 
# - Draw 500 source populations
# - For each source population, subsample 100 pairs from the source population.
# Do this at both small (20) and large (100) population sizes.

n.popns = 20 # 500
n.pairs = 20 # 500
n.loci = 25
theta = 0
wfitn = sqrt(1 / 0.14 / 2)

ressies = expand.grid(n0 = c(20, 100),
                      popn = 1:n.popns,
                      pair = 1:n.pairs) %>%
  arrange(n0, popn, pair) %>%
  mutate(p.fixn = 0,
         p.fixp = 0,
         popn.H = 0)

ressies %>% head()
ressies %>% filter(n0 %in% 20 & popn %in% 1) %>% head()
ressies %>% filter(n0 %in% 20 & popn %in% 1) %>% nrow()


row.idx = 1
for (n0 in c(20, 100)) {
  for (popn.no in 1:n.popns) {
    popn.init = init.sim(
      a = c(-1/2, 1/2),
      params = data.frame(init.row = 1e2,
                          n.loci = n.loci, 
                          n.pop0 = n0,
                          w.max = 2, 
                          theta = theta, 
                          wfitn = wfitn,
                          sig.e = 0)
    )
    for (pair.no in 1:n.pairs) {
      pair.stats = popn.init[c(
        sample(which(popn.init$fem), 1),
        sample(which(!popn.init$fem), 1)
      ), ] %>%
        select(-c(g_i, z_i, w_i, r_i, fem, gen)) %>%
        gather(key = loc.copy, value = val, -i) %>%
        mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
        group_by(i, locus) %>%
        mutate(H.indiv = !sum(val)) %>%
        group_by(locus) %>%
        summarise(hetz = mean(H.indiv),
                  fixp = all(val > 0),
                  fixn = all(val < 0)) %>%
        ungroup() %>% 
        summarise(popn.H = mean(hetz == 0.5),
                  p.fixp = mean(fixp),
                  p.fixn = mean(fixn))
      ressies$popn.H[row.idx] = pair.stats$popn.H
      ressies$p.fixp[row.idx] = pair.stats$p.fixp
      ressies$p.fixn[row.idx] = pair.stats$p.fixn
      row.idx = row.idx + 1
    }
  }
}

ressies %>%
  ggplot() +
  geom_histogram(aes(p.fixn), binwidth = 0.04) +
  facet_wrap(~ n0)

ressies %>%
  ggplot() +
  geom_histogram(aes(p.fixp), binwidth = 0.04) +
  facet_wrap(~ n0)

ressies %>%
  ggplot() +
  geom_histogram(aes(popn.H), binwidth = 0.04) +
  facet_wrap(~ n0)

# Now, what is the expected popn.H?
# For a population with p = 0.5,
# Expected heterozygosity is 2 * p * (1-p) = 0.5

ressies %>%
  group_by(n0) %>%
  mutate(mean.H = mean(popn.H)) %>%
  ggplot() +
  geom_histogram(aes(popn.H), binwidth = 0.04,
                 fill = 'gray66') +
  geom_segment(aes(x = mean.H, xend = mean.H,
                   y = 0, yend = 100)) +
  geom_segment(aes(x = 0.5, xend = 0.5,
                   y = 0, yend = 100),
               colour = 'red') +
  facet_wrap(~ n0)

### Scratch work.

popn.init = init.sim(
  a = c(-1/2, 1/2),
  params = data.frame(init.row = 1e2,
                      n.loci = n.loci, 
                      n.pop0 = 100,
                      w.max = 2, 
                      theta = theta, 
                      wfitn = wfitn,
                      sig.e = 0)
)

popn.init[c(
  sample(which(popn.init$fem), 1),
  sample(which(!popn.init$fem), 1)
), ] %>%
  select(-c(g_i, z_i, w_i, r_i, fem, gen)) %>%
  gather(key = loc.copy, value = val, -i) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(i, locus) %>%
  mutate(H.indiv = !sum(val)) %>%
  group_by(locus) %>%
  summarise(hetz = mean(H.indiv),
            fixp = all(val > 0),
            fixn = all(val < 0)) %>%
  ungroup() %>% 
  summarise(popn.H = sum(hetz == 0.5),
            p.fixp = mean(fixp),
            p.fixn = mean(fixn))
  print(n = nrow(.))
