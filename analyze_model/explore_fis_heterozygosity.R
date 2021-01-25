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

hn = htrials %>% group_by(gen, trial) %>% summarise(n = n(), g = mean(g_i))

hf = htrials %>%
  select(-c(g_i, z_i, w_i, r_i, fem)) %>%
  gather(key = loc.copy, value = val, -c(i, gen, trial)) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(i, gen, locus, trial) %>%
  mutate(heter = !sum(val)) %>%
  group_by(gen, locus, trial) %>%
  summarise(p = mean(val > 0),
            f = 1 - (mean(heter) / (2 * p * (1-p)))) %>%
  group_by(gen, trial) %>%
  summarise(f = mean(f, na.rm = TRUE))

hall = merge(hf, hn, by = c('trial', 'gen'))

head(hall)

hall %>%
  ggplot() +
  geom_line(aes(x = gen, y = n, group = trial)) +
  scale_y_log10()

hall %>%
  ggplot() +
  geom_line(aes(x = gen, y = f, group = trial))

hall = hall %>% group_by(trial) %>% mutate(extinct = max(gen) < 15) %>% ungroup()

hall %>%
  ggplot() +
  geom_line(aes(x = gen, y = f, group = trial, colour = extinct)) +
  scale_color_manual(values = c('black', 'red'))

hall %>%
  group_by(trial) %>%
  mutate(tau = max(gen) - gen + 1 + as.numeric(extinct & all(n > 1))) %>%
  ggplot() +
  geom_line(aes(x = tau, y = f, group = trial, colour = extinct), size = 0.5) +
  scale_color_manual(values = c('black', 'red')) +
  scale_x_reverse()

hall %>%
  group_by(trial) %>%
  mutate(tau = max(gen) - gen + 1 + as.numeric(extinct & all(n > 1))) %>%
  group_by(extinct, tau) %>%
  summarise(fbar = mean(f),
            fvar = var(f),
            n = n()) %>%
  ggplot() +
  geom_ribbon(aes(tau, 
                  ymin = fbar - 2*sqrt(fvar/n),
                  ymax = fbar + 2*sqrt(fvar/n),
                  group = extinct),
              alpha = 0.2) +
  geom_line(aes(tau, fbar, group = extinct, linetype = extinct)) +
  scale_x_reverse()

# lookng great
# too great...?

hall %>% filter(extinct) %>% 
  arrange(trial, desc(gen)) %>% 
  distinct(trial, .keep_all = TRUE) %>%
  summarise(fbar = mean(f)) # uh??? is this wrong?

# ohhh this will have all n = 1 at f = -1 yes?
# wat wtf is happening with negative f values..

htrials %>% filter(trial %in% 1 & gen %in% 7)
# uhh????? there's a ton of heterzygosity here...

htrials %>% 
  filter(trial %in% 1 & gen %in% 6) %>%
  select(-c(gen, g_i, z_i, w_i, r_i, fem, trial)) %>%
  gather(key = loc.copy, value = val, -i) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(i, locus) %>%
  mutate(heter = !sum(val)) %>%
  group_by(locus) %>%
  summarise(p = mean(val > 0),
            f = 1 - (mean(heter) / (2 * p * (1-p))))

 htrials %>% 
    filter(trial %in% 1) %>%
    select(-c(g_i, z_i, w_i, r_i, fem, trial)) %>%
    gather(key = loc.copy, value = val, -c(i, gen)) %>%
    mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
    group_by(i, locus, gen) %>%
    mutate(heter = !sum(val)) %>%
    group_by(locus, gen) %>%
    summarise(p = mean(val > 0),
              h = mean(heter),
              f = 1 - (mean(heter) / (2 * p * (1-p)))) %>%
   ggplot() +
   geom_line(aes(x = gen, y = h), colour = 'blue') +
   geom_line(aes(x = gen, y = p), colour = 'red') +
   facet_wrap(~ locus)
 
 htrials %>% 
   filter(trial %in% 1) %>%
   select(-c(g_i, z_i, w_i, r_i, fem, trial)) %>%
   gather(key = loc.copy, value = val, -c(i, gen)) %>%
   mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
   group_by(i, locus, gen) %>%
   mutate(heter = !sum(val)) %>%
   group_by(locus, gen) %>%
   summarise(p = mean(val > 0),
             h = mean(heter),
             f = 1 - (mean(heter) / (2 * p * (1-p)))) %>%
   ggplot() +
   geom_segment(aes(x = 0, xend = 7, y = 0, yend = 0), 
                linetype = 2) +
   geom_line(aes(x = gen, y = 1 - h/(2*p*(1-p))), colour = 'blue') +
   geom_line(aes(x = gen, y = p), colour = 'red') +
   facet_wrap(~ locus)

 htrials %>%
   filter(trial %in% 1) %>%
   select(i, gen, a14, b14, a19, b19, a20, b20) %>%
   gather(key = loc.copy, value = val, -c(i, gen)) %>%
   mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
   group_by(i, locus, gen) %>%
   mutate(heter = !sum(val)) %>%
   group_by(locus, gen) %>%
   summarise(p = mean(val > 0),
             h = mean(heter),
             f = 1 - (mean(heter) / (2 * p * (1-p))))
   
hf2 = htrials %>%
   select(-c(g_i, z_i, w_i, r_i, fem)) %>%
   gather(key = loc.copy, value = val, -c(i, gen, trial)) %>%
   mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
   group_by(i, gen, locus, trial) %>%
   mutate(heter = !sum(val)) %>%
   group_by(gen, locus, trial) %>%
   summarise(p = mean(val > 0),
             f = 1 - (mean(heter) / (2 * p * (1-p))),
             f = ifelse(is.nan(f), 1, f)) %>%
   group_by(gen, trial) %>%
   summarise(f = mean(f, na.rm = TRUE))

hf2 %>% 
  group_by(trial) %>%
  filter(max(gen) < 15) %>% 
  ungroup() %>%
  arrange(trial, desc(gen)) %>% 
  distinct(trial, .keep_all = TRUE) %>%
  summarise(fbar = mean(f))

hall2 = merge(hn, hf2, by = c('trial', 'gen'))

hall2 %>%
  group_by(trial) %>%
  mutate(extinct = max(gen) < 15,
         tau = max(gen) - gen + 1 + as.numeric(extinct & all(n > 1))) %>%
  group_by(extinct, tau) %>%
  summarise(fbar = mean(f),
            fvar = var(f),
            n = n()) %>%
  ggplot() +
  geom_ribbon(aes(tau, 
                  ymin = fbar - 2*sqrt(fvar/n),
                  ymax = fbar + 2*sqrt(fvar/n),
                  group = extinct),
              alpha = 0.2) +
  geom_line(aes(tau, fbar, group = extinct, linetype = extinct)) +
  scale_x_reverse()

hall2 %>%
  group_by(trial) %>%
  mutate(extinct = max(gen) < 15) %>%
  group_by(extinct, gen) %>%
  summarise(fbar = mean(f),
            fvar = var(f),
            n = n()) %>%
  ggplot() +
  geom_ribbon(aes(gen, 
                  ymin = fbar - 2*sqrt(fvar/n),
                  ymax = fbar + 2*sqrt(fvar/n),
                  group = extinct),
              alpha = 0.2) +
  geom_line(aes(gen, fbar, group = extinct, linetype = extinct))

# These latter plots are interesting but may just be weighted averages of F and 1...

# Above we may have detected patterns with mean(h) / mean(2pq)

hf3 = htrials %>%
  select(-c(g_i, z_i, w_i, r_i, fem)) %>%
  gather(key = loc.copy, value = val, -c(i, gen, trial)) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(i, gen, locus, trial) %>%
  mutate(heter = !sum(val)) %>%
  group_by(gen, locus, trial) %>%
  summarise(p = mean(val > 0),
            h = mean(heter),
            exph = 2 * mean(p * (1-p))) %>%
  group_by(gen, trial) %>%
  summarise(f = 1 - (mean(h) / mean(exph)))

# This also seems bad... averaging over proportions.
# (although we're also sorta averaging ovr proportions by taking mean of h across loci)

hall3 = merge(hn, hf3, by = c('trial', 'gen'))

hall3 %>%
  group_by(trial) %>%
  mutate(extinct = max(gen) < 15,
         tau = max(gen) - gen + 1 + as.numeric(extinct & all(n > 1))) %>%
  group_by(extinct, tau) %>%
  summarise(fbar = mean(f),
            fvar = var(f),
            n = n()) %>%
  ggplot() +
  geom_ribbon(aes(tau, 
                  ymin = fbar - 2*sqrt(fvar/n),
                  ymax = fbar + 2*sqrt(fvar/n),
                  group = extinct),
              alpha = 0.2) +
  geom_line(aes(tau, fbar, group = extinct, linetype = extinct)) +
  scale_x_reverse()

# seems to be same issue as above...

### Synthesis:
# Okay a couple of potntially bad problems here. First, for some reason, there's
# consistently more heterozygosity than expected here. WTAF. Over the summer
# there was an issue with inheritance. I don't know if that's what's happening
# here. This may require doing some sort of extra checking of sorts... But even
# still, F_{IS} may not be good to use here. There are two different ways to
# estimate mean heterozygosity here... One is to get mean h over all loci and
# mean 2pq over all loci, then tak that ratio. That may be missing some info,
# though, I feel. There's also the issue of averaging over ratios which seems
# problematic (would bias results somehow). If we instead calculated mean F_{IS}
# at each locus, we also may have the ratio problem... but more importantly,
# F_{IS} is undefined for allles at fixation because 2*p*(1-q) is zero. Taking
# these alleles out may be biasing, especially at low populaton size (I think).
# It may be better instead to use fixation and/or genetic variation.

