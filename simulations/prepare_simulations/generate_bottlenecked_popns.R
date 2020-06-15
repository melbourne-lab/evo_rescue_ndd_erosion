library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

rm(list = ls())

source('base_model_source/sim_functions.R')

##### Parameters

n.trials = 500
n.loci = 25
theta = 0
wfitn = sqrt(1 / 0.14 / 2)

##### Run simulations

liszt = vector('list', n.trials)

set.seed(87795)

for (trial in 1:n.trials) {
  liszt[[trial]] = sim(
    a = c(-1/2, 1/2),
    params = data.frame(end.time = 8,
                        init.row = 1e2,
                        n.loci = n.loci, 
                        n.pop0 = 2,
                        w.max = 2, 
                        theta = theta, 
                        wfitn = wfitn,
                        sig.e = 0)
  )
  print(trial)
}

all.trials = unroller(liszt)

head(all.trials)

init.pops = all.trials %>%
  group_by(trial) %>%
  filter(mean(fem[gen %in% 1]) == 0.5)

length(unique(init.pops$trial))

pop20s = init.pops %>%
  group_by(trial, gen) %>%
  mutate(n = n()) %>%
  group_by(trial) %>%
  filter(any(n > 20)) %>%
  ungroup() %>%
  select(-n)

pop100s = init.pops %>%
  group_by(trial, gen) %>%
  mutate(n = n()) %>%
  group_by(trial) %>%
  filter(any(n > 100)) %>%
  ungroup() %>%
  select(-n)

nrow(pop20s)
length(unique(pop20s$trial))
nrow(pop100s)
length(unique(pop100s$trial))

pop100s %>%
  group_by(trial, gen) %>%
  summarise(n = n()) %>%
  group_by(trial) %>%
  filter(n == max(n)) %>%
  ungroup() %>%
  select(n) %>% unlist() %>% hist()

# One question that just occurred to me... does it matter when a population is selected? Should I sample from the first time step with >100 individuals, or all at the same timestep?

var.t.20s = pop20s %>%
  select(-c(g_i, z_i, w_i, r_i, fem)) %>%
  gather(key = loc.copy, value = val, -c(trial, i, gen)) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(trial, gen, locus) %>%
  summarise(p = mean(val > 0),
            n = length(unique(i))) %>%
  group_by(trial, gen) %>%
  summarise(var.t = sum(2 * p * (1 - p) / 25),
            n.t = n[1])

var.t.20s %>%
  group_by(gen) %>%
  mutate(vt.bar = mean(var.t)) %>%
  ggplot() +
  geom_segment(aes(x = 1, xend = 8, y = 0.48, yend = 0.48)) +
  geom_line(aes(x = gen, y = var.t, group = trial)) +
  geom_point(aes(x = gen, y = var.t, fill = log(n.t / 20)),
             shape = 21, colour = 'black', stroke = 1) +
  geom_line(aes(x = gen, y = vt.bar), linetype = 3,
            colour = 'red', size = 3) +
  scale_fill_gradient2(low = 'darkblue', mid = 'white', high = 'red') +
  theme(panel.background = element_rect('white'))

# Do the same but for the populations of size 100

var.t.100s = pop100s %>%
  select(-c(g_i, z_i, w_i, r_i, fem)) %>%
  gather(key = loc.copy, value = val, -c(trial, i, gen)) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(trial, gen, locus) %>%
  summarise(p = mean(val > 0),
            n = length(unique(i))) %>%
  group_by(trial, gen) %>%
  summarise(var.t = sum(2 * p * (1 - p) / 25),
            n.t = n[1])

var.t.100s %>%
  group_by(gen) %>%
  mutate(vt.bar = mean(var.t)) %>%
  ggplot() +
  geom_segment(aes(x = 1, xend = 8, y = 0.48, yend = 0.48)) +
  geom_line(aes(x = gen, y = var.t, group = trial)) +
  geom_point(aes(x = gen, y = var.t, fill = log(n.t / 100)),
             shape = 21, colour = 'black', stroke = 1) +
  geom_line(aes(x = gen, y = vt.bar), linetype = 3,
            colour = 'red', size = 3) +
  scale_fill_gradient2(low = 'darkblue', mid = 'white', high = 'red') +
  theme(panel.background = element_rect('white'))

nrow(pop20s)

nrow(pop100s)

## Looking at the proposal, it specifies having inbred population lines bred out to N = 1000 individuals. Try doing this here, adding an extra time step (will slow things down).

##### Parameters

n.trials = 500
n.loci = 25
theta = 0
wfitn = sqrt(1 / 0.14 / 2)

##### Run simulations

liszt = vector('list', n.trials)

set.seed(87795)

for (trial in 1:n.trials) {
  liszt[[trial]] = sim(
    a = c(-1/2, 1/2),
    params = data.frame(end.time = 9,
                        init.row = 1e2,
                        n.loci = n.loci, 
                        n.pop0 = 2,
                        w.max = 2, 
                        theta = theta, 
                        wfitn = wfitn,
                        sig.e = 0)
  )
  print(trial)
}

all.trials = unroller(liszt)

pop.k.s = all.trials %>%
  group_by(trial, gen) %>%
  mutate(n = n()) %>%
  group_by(trial) %>%
  filter(any(n > 1000)) %>%
  ungroup() %>%
  select(-n)

length(unique(pop.k.s$trial))
nrow(pop.k.s)

pop.k.s %>%
  group_by(trial, gen) %>%
  summarise(n = n())

var.k = pop.k.s %>%
  select(-c(g_i, z_i, w_i, r_i, fem)) %>%
  gather(key = loc.copy, value = val, -c(trial, i, gen)) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(trial, gen, locus) %>%
  summarise(p = mean(val > 0),
            n = length(unique(i))) %>%
  group_by(trial, gen) %>%
  summarise(var.t = sum(2 * p * (1 - p) / 25),
            n.t = n[1])

var.k %>%
  group_by(gen) %>%
  mutate(vt.bar = mean(var.t)) %>%
  ggplot() +
  geom_segment(aes(x = 1, xend = 8, y = 0.48, yend = 0.48)) +
  geom_line(aes(x = gen, y = var.t, group = trial)) +
  geom_point(aes(x = gen, y = var.t, fill = n.t > 1000),
             shape = 21, colour = 'black', stroke = 1) +
  geom_line(aes(x = gen, y = vt.bar), linetype = 3,
            colour = 'red', size = 3) +
  #scale_fill_gradient2(low = 'darkblue', mid = 'white', high = 'red') +
  scale_fill_manual(values = c('black', 'lightblue')) +
  theme(panel.background = element_rect('white'))

# With the seed provided, I have eight of these populations.


# Well - for now just use these I suppose.

write.csv(pop.k.s, na = '', row.names = FALSE,
          file = 'simulations/prepare_simulations/outputs/eight_bottlenecked_popns.R')

# Check: Linkage disequilibrium.

