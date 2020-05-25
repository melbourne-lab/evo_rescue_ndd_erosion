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
  mutate(v = sum(2 * p * (1-p) / params$n.loci))

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
         pp2 = ifelse(gen > 1, N_p / N, p))
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

# No! Okay. Incorrect. But close?
# There is convergence for high $p$., but it's off for lower p.

gg.m1 +
  geom_point(aes(x = p, y = pp2 - p, alpha = factor(gen))) +
  scale_alpha_manual(values = c(0, 1)) +
  theme_bw()
# this looks a little... Gaussian.
# but the peak in in a weird place...

m1.preds %>%
  filter(gen > 1) %>%
  ungroup() %>%
  arrange(p) %>%
  ggplot() +
  geom_line(aes(x = p, y = pp2 - p)) +
  xlim(c(0, 1))
# yes that looks gaussian to me...
# but we know it has to be at zero on both ends...

gg.m1 +
  geom_point(aes(x = p, y = ((pp2 - p) / ((1-p))) , alpha = factor(gen))) +
  scale_alpha_manual(values = c(0, 1)) +
#  coord_fixed() +
  theme_bw()


### Next, try this with m = 2.

param2 = data.frame(end.time = 2,
                    init.row = 1e4,
                    n.loci = 2, 
                    n.pop0 = 20,
                    w.max = 2, 
                    theta = 2, 
                    wfitn = sqrt(1 / 0.14 / 2),
                    sig.e = 0)

n.trials = 2000

### Run simulation 

lisz2 = vector(mode = 'list', length = n.trials)

set.seed(5280)

for (trial in 1:n.trials) {
  lisz2[[trial]] = sim(
    a = c(-1/2, 1/2),
    params = param2
  )  
  print(trial)
}

# Put all trials in one (even larger) object.
m2 = unroller(lisz2)

m2.pvz = m2 %>%
  select(-c(g_i, w_i, r_i, fem)) %>%
  gather(key = all.copy, value = gval, -c(i, z_i, gen, trial)) %>%
  mutate(locus = gsub('^[ab]', '', all.copy)) %>%
  group_by(trial, gen, locus) %>%
  summarise(p = mean(gval > 0),
            z = mean(z_i)) %>%
  group_by(trial, gen) %>%
  mutate(v = sum(2 * p * (1-p) / param2$n.loci))

head(m2.pvz)

m2.preds = m2.pvz %>%
  mutate(m = param2$n.loci,
         w = param2$wfitn,
         d = param2$theta - z,
         q = 1 - p,
         w_A = exp(-(d - 1/sqrt(m))^2 / (2*w^2)),
         w_0 = exp(-d^2 / (2*w^2)),
         w_a = exp(-(d + 1/sqrt(m))^2 / (2*w^2)),
         N_p = p * (w_A*p*(p+1) + w_0*q*(2*p+1) + w_a * q^2),
         N   = 2 * (p^2 * w_A + 2*p*q*w_0 + q^2 * w_a),
         pp2 = ifelse(gen > 1, N_p / N, p))

head(m2.preds)

gg.m2 = m2.preds %>% ggplot()

gg.m2 +
  geom_point(aes(x = gen, y = p),
             position = position_jitter(width = 1/3)) +
  stat_smooth(aes(x = gen, y = p), method = lm, 
              size = 2, colour = 'skyblue') +
  theme_bw()

gg.m2 +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               size = 0.5, colour = 'red') +
  geom_point(aes(x = p, y = pp2, alpha = factor(gen))) +
  scale_alpha_manual(values = c(0, 1)) +
  theme_bw()

# This is noisier, but still close. Following a similar pattern...

gg.m2 +
  geom_point(aes(x = p, y = pp2 - p, alpha = factor(gen))) +
  scale_alpha_manual(values = c(0, 1)) +
  theme_bw()
# lso looks Gaussian, although the variance depends on p
# (or something like p*(1-p)... but not quite?)
# this function maxes out at p approx 1/3

gg.m2 +
  geom_point(aes(x = p, y = (pp2 - p) / ((1 - p)), alpha = factor(gen))) +
  scale_alpha_manual(values = c(0, 1)) +
  theme_bw()

### Well... what happens with m = 3?
# Expect variance to increase but will the peak residual also increase?

lisz3 = vector(mode = 'list', length = n.trials)

set.seed(111997)

for (trial in 1:n.trials) {
  lisz3[[trial]] = sim(
    a = c(-1/2, 1/2),
    params = param3
  )  
  print(trial)
}

# Put all trials in one (even largerer!) object.
m3 = unroller(lisz3)

m3.pvz = m3 %>%
  select(-c(g_i, w_i, r_i, fem)) %>%
  gather(key = all.copy, value = gval, -c(i, z_i, gen, trial)) %>%
  mutate(locus = gsub('^[ab]', '', all.copy)) %>%
  group_by(trial, gen, locus) %>%
  summarise(p = mean(gval > 0),
            z = mean(z_i)) %>%
  group_by(trial, gen) %>%
  mutate(v = sum(2 * p * (1-p) / param3$n.loci))

head(m3.pvz)

m3.preds = m3.pvz %>%
  mutate(m = param3$n.loci,
         w = param3$wfitn,
         d = param3$theta - z,
         q = 1 - p,
         w_A = exp(-(d - 1/sqrt(m))^2 / (2*w^2)),
         w_0 = exp(-d^2 / (2*w^2)),
         w_a = exp(-(d + 1/sqrt(m))^2 / (2*w^2)),
         N_p = p * (w_A*p*(p+1) + w_0*q*(2*p+1) + w_a * q^2),
         N   = 2 * (p^2 * w_A + 2*p*q*w_0 + q^2 * w_a),
         pp2 = ifelse(gen > 1, N_p / N, p))

gg.m3 = m3.preds %>% ggplot()

gg.m3 +
  geom_point(aes(x = p, y = pp2 - p, alpha = factor(gen))) +
  scale_alpha_manual(values = c(0, 1)) +
  theme_bw()
# It does.
# I bet this asymptotically approaches p = 1/2.
# But still... what the heck.

gg.m3 +
  geom_point(aes(x = p, y = (pp2 - p)/((1-p)), alpha = factor(gen))) +
  scale_alpha_manual(values = c(0, 1)) +
  theme_bw()

### Try wth a larger number (m = 100)
parami = data.frame(end.time = 2,
                    init.row = 1e4,
                    n.loci = 100, 
                    n.pop0 = 20,
                    w.max = 2, 
                    theta = 2, 
                    wfitn = sqrt(1 / 0.14 / 2),
                    sig.e = 0)

n.trials = 2000

### Run simulation 

liszi = vector(mode = 'list', length = n.trials)

set.seed(811990)

for (trial in 1:n.trials) {
  liszi[[trial]] = sim(
    a = c(-1/2, 1/2),
    params = parami
  )  
  print(trial)
}

# Put all trials in one (even largerer!) object.
mi = unroller(liszi)

mi.pvz = mi %>%
  select(-c(g_i, w_i, r_i, fem)) %>%
  gather(key = all.copy, value = gval, -c(i, z_i, gen, trial)) %>%
  mutate(locus = gsub('^[ab]', '', all.copy)) %>%
  group_by(trial, gen, locus) %>%
  summarise(p = mean(gval > 0),
            z = mean(z_i)) %>%
  group_by(trial, gen) %>%
  mutate(v = sum(2 * p * (1-p) / parami$n.loci))

head(mi.pvz)

mi.preds = mi.pvz %>%
  mutate(m = parami$n.loci,
         w = parami$wfitn,
         d = parami$theta - z,
         q = 1 - p,
         w_A = exp(-(d - 1/sqrt(m))^2 / (2*w^2)),
         w_0 = exp(-d^2 / (2*w^2)),
         w_a = exp(-(d + 1/sqrt(m))^2 / (2*w^2)),
         N_p = p * (w_A*p*(p+1) + w_0*q*(2*p+1) + w_a * q^2),
         N   = 2 * (p^2 * w_A + 2*p*q*w_0 + q^2 * w_a),
         pp2 = ifelse(gen > 1, N_p / N, p))

write.csv(mi.preds %>% filter(as.numeric(locus) < 21),
          row.names = FALSE,
          file = 'm_100_locus_1-20_gen2.csv')

gg.mi = mi.preds %>% filter(as.numeric(locus) < 25) %>% ggplot()

gg.mi +
  geom_point(aes(x = p, y = pp2 - p, alpha = factor(gen), colour = v)) +
  scale_alpha_manual(values = c(0, 1)) +
  theme_bw()

gg.mi +
  geom_point(aes(x = p, y = (pp2 - p)/(p*(1-p)), alpha = factor(gen))) +
  scale_alpha_manual(values = c(0, 1)) +
  theme_bw()
