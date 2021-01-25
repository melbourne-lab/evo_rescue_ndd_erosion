# Trial run at calculating and plotting F_{IS} over tim.

library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

source('base_model_source/sim_functions.R')

pars = data.frame(
  end.time = 15,
  init.row = 1e4,
  n.pop0 = 50,
  n.loci = 25,
  w.max = 2,
  theta = 2.75,
  wfitn = sqrt(1 / 0.14 / 2),
  sig.e = sqrt(0.5),
  pos.p = 0.5,
  alpha = 0
)

set.seed(9900)

sim.test = sim(params = pars, evolve = TRUE)

table(sim.test$gen)

sim.test %>%
  select(-c(i, g_i, z_i, w_i, r_i, fem)) %>%
  gather(key = loc.copy, value = val, -gen) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(gen, locus) %>%
  summarise(p = mean(val > 0)) %>%
  group_by(gen) %>%
  summarise(p.fix.pos = mean(p == 1),
            p.fix.neg = mean(p == 0),
            v = sum(2 * p * (1 - p)) / pars$n.loci[1])

sim.test %>%
  select(-c(g_i, z_i, w_i, r_i, fem)) %>%
  gather(key = loc.copy, value = val, -c(i, gen)) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(gen, i, locus) %>%
  mutate(heter = !sum(val)) %>%
  group_by(gen, locus) %>%
  summarise(p = mean(val > 0),
            h = mean(heter)) %>%
  group_by(gen) %>%
  summarise(p.fix.pos = mean(p == 1),
            p.fix.neg = mean(p == 0),
            h = 1 - (mean(h) / mean(2 * p * (1-p))),
            v = sum(2 * p * (1 - p)) / pars$n.loci[1])

uhm = sim.test %>%
  select(-c(g_i, z_i, w_i, r_i, fem)) %>%
  gather(key = loc.copy, value = val, -c(i, gen)) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(gen, i, locus) %>%
  mutate(heter = !sum(val)) %>%
  group_by(gen, locus) %>%
  summarise(p = mean(val > 0),
            h = mean(heter),
            pq2 = 2 * p * (1-p))

# uhm...

uhm %>%
  ggplot(aes(x = pq2, y = h)) +
  geom_line(aes(group = locus)) +
  geom_point(aes(colour = gen))

uhm %>%
  mutate(fis = 1 - h/pq2) %>%
  ggplot() +
  geom_point(aes(x = gen, y = fis)) +
  geom_line(aes(x = gen, y = fis, group = locus))

sim.test %>%
  filter(gen %in% 1) %>%
  filter((a1 < 0 & b1 > 0) | (a1 > 0 & b1 < 0))

### Look at sampling distribution of h in t = 0

init.trials = vector('list', 100)

for (trial in 1:100) {
  init.trials[[trial]] = init.sim(params = pars) %>%
    select(-c(gen, g_i, z_i, w_i, r_i, fem)) %>%
    gather(key = loc.copy, value = val, -c(i)) %>%
    mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
    group_by(i, locus) %>%
    mutate(heter = !sum(val)) %>%
    group_by(locus) %>%
    summarise(p = mean(val > 0),
              h = mean(heter),
              pq2 = 2 * p * (1-p))
}

init.arr = do.call(rbind, init.trials)

init.arr %>% print(n = 30)

init.arr %>% mutate(idiff = c(0, diff(as.numeric(locus)))) %>% mutate(trial = 1 + cumsum(idiff == -8)) %>% print(n = 30)

init.arr = init.arr %>%
  mutate(idiff = c(0, diff(as.numeric(locus)))) %>% 
  mutate(trial = 1 + cumsum(idiff == -8)) %>% 
  select(-idiff)

table(init.arr$trial)

init.arr %>%
  group_by(trial) %>%
  summarise(hbar = mean(h),
            pq2b = mean(pq2)) %>%
  mutate(fis = 1 - hbar/pq2b) %>%
  summarise(mean(fis), sd(fis), mean(fis) + 2 * sqrt(sd(fis) / 50))

### Hmm... okay try several of these.

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

hn = htrials %>% group_by(gen, trial) %>% summarise(n = n())

hf = htrials %>%
  select(-c(g_i, z_i, w_i, r_i, fem)) %>%
  gather(key = loc.copy, value = val, -c(i, gen, trial)) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(i, gen, locus, trial) %>%
  mutate(heter = !sum(val)) %>%
  group_by(gen, locus, trial) %>%
  summarise(p = mean(val > 0),
            h = mean(heter),
            pq2 = 2 * p * (1-p)) %>%
  group_by(gen, trial) %>%
  summarise(hbar = mean(h),
            hvar = var(h),
            pbar = mean(pq2),
            pvar = var(pq2),
            nobs = n())

hall = merge(hf, hn, by = c('trial', 'gen'))

head(hall)

hall %>%
  ggplot() +
  geom_line(aes(x = gen, y = n, group = trial)) +
  scale_y_log10()

hall %>%
  ggplot() +
  geom_line(aes(x = gen, y = hbar, group = trial))

hall %>%
  ggplot() +
  geom_line(aes(x = gen, y = 1 - hbar/pbar, group = trial))

hall = hall %>% group_by(trial) %>% mutate(extinct = max(gen) < 15)

hall %>%
  ggplot() +
  geom_line(aes(x = gen, y = 1 - hbar/pbar, group = trial, colour = extinct)) +
  scale_color_manual(values = c('black', 'red'))

hall %>%
  group_by(trial) %>%
  mutate(tau = max(gen) - gen + 1 + as.numeric(extinct & all(n > 1))) %>%
  ggplot() +
  geom_line(aes(x = tau, y = 1 - hbar/pbar, group = trial, colour = extinct), size = 0.5) +
  scale_color_manual(values = c('black', 'red')) +
  scale_x_reverse()

hall %>%
  group_by(trial) %>%
  mutate(tau = max(gen) - gen + 1 + as.numeric(extinct & all(n > 1))) %>%
  group_by(extinct, tau) %>%
  summarise(fisbar = mean(1 - (hbar / pbar)),
            fisvar = var(1 - (hbar / pbar)),
            n = n()) %>%
  ggplot() +
  geom_ribbon(aes(tau, 
                  ymin = fisbar - 2*sqrt(fisvar/n),
                  ymax = fisbar + 2*sqrt(fisvar/n),
                  group = extinct),
              alpha = 0.2) +
  geom_line(aes(tau, fisbar, group = extinct, linetype = extinct)) +
  scale_x_reverse()

