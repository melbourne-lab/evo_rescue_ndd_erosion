# Trial run at calculating and plotting 3D trajectories of populations

library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)
library(rgl)

source('base_model_source/sim_functions.R')

pars = data.frame(
  end.time = 15,
  init.row = 1e4,
  n.pop0 = 30,
  n.loci = 25,
  w.max = 2,
  theta = 3,
  wfitn = sqrt(1 / 0.14 / 2),
  sig.e = sqrt(0.5),
  pos.p = 0.5,
  alpha = 0
)

htrials = vector('list', 100)

set.seed(101)

for (trial in 1:100) {
  htrials[[trial]] = sim(params = pars)
  print(trial)
}

htrials = unroller(htrials)

# Three variables: N, sigma, W?

h1 = htrials %>% filter(trial %in% 1)

# h1 demo
h1.demo = h1 %>%
  group_by(gen) %>%
  summarise(n = n(),
            w = mean(w_i))

h1.gene = h1 %>%
  select(-c(i, g_i, z_i, w_i, r_i, fem)) %>%
  gather(key = loc.copy, value = val, -c(gen)) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(gen, locus) %>%
  summarise(p = mean(val > 0)) %>%
  group_by(gen) %>%
  summarise(p.fix.pos = mean(p == 1),
            p.fix.neg = mean(p == 0),
            v = sum(2 * p * (1 - p)) / pars$n.loci[1])

h1.gene

h1.all = merge(h1.demo, h1.gene)

h1.plot = h1.all %>% 
  rbind(., ., .) %>% 
  arrange(gen) %>%
  mutate(v = c(v[3:(nrow(.))], NA, NA),
         w = c(w[2:(nrow(.))], NA)) %>%
  filter(!is.na(n))

h1.plot

with(h1.plot, plot3d(n, v, w, type = 'l'))

# h1 demo
h.demo = htrials %>%
  group_by(trial, gen) %>%
  summarise(n = n(),
            w = mean(w_i))

h.gene = htrials %>%
  select(-c(i, g_i, z_i, w_i, r_i, fem)) %>%
  gather(key = loc.copy, value = val, -c(gen, trial)) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(gen, locus, trial) %>%
  summarise(p = mean(val > 0)) %>%
  group_by(gen, trial) %>%
  summarise(p.fix.pos = mean(p == 1),
            p.fix.neg = mean(p == 0),
            v = sum(2 * p * (1 - p)) / pars$n.loci[1])

h.all = merge(h.demo, h.gene) %>%
  group_by(trial) %>%
  mutate(extinct = max(gen) < 15) %>%
  ungroup()

head(h.all)

h.plot = h.all %>%
  rbind(., ., .) %>%
  arrange(trial, gen) %>%
  group_by(trial) %>%
  mutate(v = c(v[3:(n())], NA, NA),
         w = c(w[2:(n())], NA)) %>%
  filter(!is.na(v))

h.plot

with(h.plot[h.plot$trial %in% 10:20,], plot3d(n, w, v, xlab = c(0, 60), type = 'n'))

rbo = rainbow(1:pars$end.time, start = 0.5, end = 1)

for (k in 1:4) {
  with(h.plot[h.plot$trial %in% k,], 
       plot3d(n, w, v, type = 'l', add = TRUE, lwd = 2,
              col = c('black', 'red')[1 + as.numeric(extinct)]))
  with(h.plot[h.plot$trial %in% k,],
       spheres3d(n, w, v, radius = 0.1, cols = rbo[gen]))
}

# Not incredibly helpful.
# Try log lambda?

