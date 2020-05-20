##### A script to look at the first time step (again)
##### Here, I want to model changes in allele freqs. over itme.
##### Do just the first time step to initialize w/ similar allele freqs., pop'n
##### size, avoid extinction.

library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

### Parameters used

params = data.frame(end.time = 2,
                    init.row = 1e4,
                    n.loci = 1, 
                    n.pop0 = 20,
                    w.max = 2, 
                    theta = 2, 
                    wfitn = sqrt(1 / 0.14 / 2),
                    sig.e = 0)

n.trials = 2000

# Note:
#   end.time = 2
#   n.loci = 1
#   theta = 2 (not 2.6, just makes the math slightly simpler)

### Source files
source('base_model_source/sim_functions.R')


### Run simulation 

liszt = vector(mode = 'list', length = n.trials)

set.seed(85200)

for (trial in 1:n.trials) {
  liszt[[trial]] = sim(
    a = c(-1/2, 1/2),
    params = params
  )  
  print(trial)
}

# Put all trials in one (large) object.
m1 = unroller(liszt)

### Analysis

m1.pvz = m1 %>%
  select(-c(g_i, w_i, r_i, fem)) %>%
  gather(key = all.copy, value = gval, -c(i, z_i, gen, trial)) %>%
  mutate(locus = gsub('^[ab]', '', all.copy)) %>%
  group_by(trial, gen, locus) %>%
  summarise(p = mean(gval > 0),
            z = mean(z_i)) %>%
  group_by(trial, gen) %>%
  mutate(v = sum(2 * p * (1-p) / sqrt(params$n.loci)))

# in each trial, for each locus, in each generation, records
#   p - mean frequency of positive allele at that locus
#   z - mean genotype
#   v - overall amount of genetic variance

head(m1.pvz)


# Now, what do we want to test?
# Want to use the expression from previously to see what the mean allele
# frequency should be in gen 2.

m1.preds = m1.pvz %>%
  mutate(m = params$n.loci,
         w = params$wfitn,
         d = params$theta - z,
         q = 1 - p,
         w_A = exp(-(d - 1/sqrt(m))^2 / (2*w^2)),
         w_0 = exp(-d^2 / (2*w^2)),
         w_a = exp(-(d + 1/sqrt(m))^2 / (2*w^2)),
         N_p = p * (w_A*p*(p+1) + w_0*q*(2*p+1) + w_a * q^2),
         N   = 2 * (p^2 * w_A + 2*p*q*w_0 + q^2 * w_a),
         pp2 = ifelse(p > 1, N_p / N, p))
# Has a lot of extra stuff, but:
# pp2 is predicted p at time 2
# (for time-step 1, just let it be pp2 = p.)

head(m1.preds)

gg.m1 = m1.preds %>% ggplot()

gg.m1 +
  geom_point(aes(x = gen, y = p),
             position = position_jitter(width = 1/3)) +
  stat_smooth(aes(x = gen, y = p), method = lm, 
              size = 2, colour = 'skyblue') +
  theme_bw()

gg.m1 +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               size = 0.5, colour = 'red') +
  geom_point(aes(x = p, y = pp2, alpha = factor(gen))) +
  scale_alpha_manual(values = c(0, 1)) +
  theme_bw()

# YES. On the one-to-one line. This means that the formula
# *exactly* predicts the change in allele frequencies.

### Next, try this with m = 2.