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

#####
##### 3D plots
#####

h.plot

with(h.plot[h.plot$trial %in% 1:10,], plot3d(n, w, v, xlim = c(0, 60), type = 'n'))

rbo = rainbow(1:pars$end.time, start = 0.5, end = 1)

for (k in 1:10) {
  with(h.plot[h.plot$trial %in% k,], 
       plot3d(n, w, v, type = 'l', add = TRUE, lwd = 2,
              col = c('black', 'red')[1 + as.numeric(extinct)]))
  with(h.plot[h.plot$trial %in% k,],
       spheres3d(n, w, v, radius = 0.1, cols = rbo[gen]))
}

# Not incredibly helpful.
# Try log lambda?

#####
##### What about a 2d plot?
#####

# h.plot2d = h.plot %>% distinct(gen, n, v, .keep_all = TRUE)

h.plot2d = rbind(h.all %>% mutate(win.n = 1), h.all %>% mutate(win.n = 2)) %>%
  arrange(trial, gen) %>%
  group_by(trial) %>%
  mutate(v = c(NA, v[1:(n()-1)])) %>%
  filter(!is.na(v))

# h.plot2d %>%
#   filter(trial %in% 1:20) %>%
#   ggplot(aes(x = n, y = v)) +
#   geom_line(aes(group = trial, colour = extinct)) +
#   xlim(c(0, 60))

# Plot with slanted lines

h.plot2d %>%
  ggplot(aes(x = n, y = v)) +
  geom_line(aes(group = trial, colour = extinct),
            size = 0.2) +
  scale_x_log10() +
  scale_color_manual(values = c('black', 'red'))


# Plot with points (and failed-stair-step)
h.plot2d %>%
  ggplot(aes(x = n, y = v)) +
  geom_line(aes(group = trial, colour = extinct),
            size = 0.05) +
  geom_point(aes(colour = extinct),
             size = 0.3) +
  scale_x_log10() +
  scale_color_manual(values = c('black', 'red'))

# Plot of just one population over time 
# (to show what the problem is with stair-steps)
h.plot2d %>%
  filter(trial %in% 1) %>%
  ggplot(aes(x = v, y = n)) +
  geom_line(aes(group = trial, colour = extinct),
            size = 0.5) +
  geom_point(aes(colour = extinct),
             size = 0.3) +
  scale_y_log10() +
  scale_color_manual(values = c('black', 'red'))

# Plotting just the first 30 trials, diagonal lines
h.all %>%
  filter(trial %in% 1:30) %>%
  ggplot(aes(x = v, y = n)) +
  geom_line(aes(group = trial, colour = extinct),
            size = 0.3) +
  geom_point(aes(colour = extinct),
             size = 0.4) +
  scale_y_log10() +
  scale_color_manual(values = c('black', 'red'))

### Sim some populations for Brett to explore.

rm(list = ls())

# Source
source('base_model_source/sim_functions.R')

# Pars, density independnt 
pars = data.frame(
  end.time = 10,
  init.row = 1e4,
  n.pop0 = 30,
  n.loci = 25,
  w.max = 2,
  theta = 2.75,
  wfitn = sqrt(1 / 0.14 / 2),
  sig.e = sqrt(0.5),
  pos.p = 0.5,
  alpha = 0
)

htrials.di = vector('list', 100)

set.seed(101)

for (trial in 1:100) {
  htrials.di[[trial]] = sim(params = pars)
  print(trial)
}

htrials.di = unroller(htrials.di)
# Now run density independent trials
pars = pars %>% mutate(alpha = 0.0035)

htrials.dd = vector('list', 100)

set.seed(7017)

for (trial in 1:100) {
  htrials.dd[[trial]] = sim(params = pars)
  print(trial)
}

htrials.dd = unroller(htrials.dd)

h.trials = rbind(
  htrials.di %>% mutate(alpha = 0), 
  htrials.dd %>% mutate(alpha = 0.0035)
)

h.demo = h.trials %>%
  group_by(trial, alpha, gen) %>%
  summarise(n = n(),
            g = mean(g_i),
            z = mean(z_i),
            w = mean(w_i))

h.gene = h.trials %>%
  select(-c(i, g_i, z_i, w_i, r_i, fem)) %>%
  gather(key = loc.copy, value = val, -c(gen, trial, alpha)) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(gen, locus, trial, alpha) %>%
  summarise(p = mean(val > 0)) %>%
  group_by(gen, trial, alpha) %>%
  summarise(p.fix.pos = mean(p == 1),
            p.fix.neg = mean(p == 0),
            v = sum(2 * p * (1 - p)) / pars$n.loci[1])

h.all = merge(h.demo, h.gene) %>% 
  arrange(trial, alpha, gen) %>%
  group_by(trial, alpha) %>%
  mutate(extinct = max(gen) < 10 | any(n < 2)) %>%
  ungroup()

head(h.all)

with(h.all[h.all$trial %in% 1:10,], plot3d(n, w, v, xlim = c(0, 60), type = 'n'))

rbo = rainbow(1:pars$end.time, start = 0.5, end = 1)

for (k in 1:10) {
  with(h.all[h.all$trial %in% k & !h.all$alpha,], 
       plot3d(n, w, v, type = 'l', add = TRUE, lwd = 2, xlim = c(0, 60),
              col = c('black', 'red')[1 + as.numeric(extinct)]))
  # with(h.plot[h.plot$trial %in% k,],
       # spheres3d(n, w, v, radius = 0.1, cols = rbo[gen]))
}

# For stairstep plot

h.plot = h.all %>%
  rbind(., ., .) %>%
  arrange(trial, alpha, gen) %>%
  group_by(trial, alpha) %>%
  mutate(n  = c(n[2:(n())], NA),
         v  = c(w[3:(n())], NA, NA),
         fn = c(p.fix.neg[3:(n())], NA, NA)) %>%
  filter(!is.na(w)) %>%
  arrange(trial, alpha, gen)

# Le stairstep plot
with(h.plot[h.plot$trial %in% 1:15,], plot3d(n, w, v, xlim = c(0, 60), add = FALSE, type = 'n'))

for (k in 1:15) {
  with(h.plot[h.plot$trial %in% k & !h.plot$alpha,], 
       plot3d(n, w, v, type = 'l', add = TRUE, lwd = 1,
              col = c('black', 'red')[1 + as.numeric(extinct)],
              xlim = c(0, 60)))
  # with(h.plot[h.plot$trial %in% k,],
  # spheres3d(n, w, v, radius = 0.1, cols = rbo[gen]))
}

with(h.plot[h.plot$trial %in% 1:15,], plot3d(n, w, v, xlim = c(0, 60), add = FALSE, type = 'n'))

for (k in 1:15) {
  with(h.plot[h.plot$trial %in% k & !!h.plot$alpha,], 
       plot3d(n, w, v, type = 'l', add = TRUE, lwd = 2,
              col = c('black', 'red')[1 + as.numeric(extinct)],
              xlim = c(0, 60), size = 0.5))
  # with(h.plot[h.plot$trial %in% k,],
  # spheres3d(n, w, v, radius = 0.1, cols = rbo[gen]))
}

### Make summaries (for extinct vs. surviving)

h.summ = h.all %>%
  group_by(trial, alpha) %>%
  mutate(tau = max(gen) - gen + 1 - as.numeric(any(n < 2))) %>%
  group_by(alpha, extinct, tau) %>%
  summarise(nbar = mean(n), nvar = var(n),
            wbar = mean(w), wvar = var(w),
            vbar = mean(v), vvar = var(v),
            pnbar = mean(p.fix.neg),
            n.trials = n())

# Note: doesn't carry variances over

h.summ.plot = h.summ %>%
  rbind(., ., .) %>%
  arrange(alpha, extinct, tau) %>%
  group_by(alpha, extinct) %>%
  mutate(nbar  = c(nbar[2:(n())], NA),
         vbar  = c(vbar[3:(n())], NA, NA),
         pnbar = c(pnbar[3:(n())], NA, NA)) %>%
  filter(!is.na(vbar)) %>%
  arrange(alpha, extinct, tau)

with(h.summ.plot, plot3d(nbar, wbar, vbar, xlim = c(0, 60), type = 'n'))

with(h.summ.plot[!h.summ.plot$alpha & h.summ.plot$extinct,], 
     plot3d(nbar, wbar, vbar, 
            xlim = c(0, 60), 
            type = 'l',
            col = 'red',
            lwd = 2)
     )
         
with(h.summ.plot[!h.summ.plot$alpha & !h.summ.plot$extinct,], 
     plot3d(nbar, wbar, vbar, 
            xlim = c(0, 60), 
            type = 'l',
            lwd = 2)
)

with(h.summ.plot[h.summ.plot$alpha & h.summ.plot$extinct,], 
     plot3d(nbar, wbar, vbar, 
            xlim = c(0, 60), 
            type = 'l',
            col = 'blue',
            lwd = 2)
)

with(h.summ.plot[h.summ.plot$alpha & !h.summ.plot$extinct,], 
     plot3d(nbar, wbar, vbar, 
            xlim = c(0, 60), 
            type = 'l',
            lwd = 2,
            col = 'purple')
)
