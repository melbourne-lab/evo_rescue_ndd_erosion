### Setup

# Load packages
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

# Clear namespace
rm(list = ls())

# Load in all simulation data
all.data = read.csv('simulations/outputs/final_results/alldata_combined.csv')

### Filter out extinctions

ext.data = all.data %>% filter(extinct)

nrow(ext.data)

# Add `tau` (time until extinction)
# extinction is population of size 0 or 1
# some populations go extinct without having a record of 0/1 so I will add that in manually
# (add in 1 because to makes log pop size easier to calculate)
ext.mod = ext.data %>%
  arrange(n.pop0, low.var, alpha, trial, desc(gen)) %>%
  distinct(n.pop0, low.var, alpha, trial, .keep_all = TRUE) %>%
  filter(n > 1) %>%
  mutate(gen = gen + 1,
         ext.gen = gen) %>%
  mutate(n = 1) %>%
  mutate_at(vars(gbar, zbar, wbar, p.fix.pos, p.fix.neg, v),
            list(~ NA)) %>%
  rbind(ext.data) %>%
  arrange(n.pop0, low.var, alpha, trial, gen)

head(ext.mod, 20)

ext.tau = ext.mod %>%
  group_by(n.pop0, alpha, low.var, trial) %>%
  mutate(log.lam = c(diff(log(n)), NA)) %>%
  mutate(tau = max(gen) - gen)%>%
  ungroup() %>%
  filter(!is.na(log.lam))

# potential confounders:
#   generation/time of extinction

ext.tau %>%
  ggplot(aes(x = tau, y = log.lam)) +
  geom_point(
    alpha = 0.2,
    position = position_jitter(width = 0.2)
  ) +
  stat_smooth() +
  facet_wrap(paste(n.pop0, low.var, sep = ', ') ~ alpha, ncol = 2)

ext.tau.summ = ext.tau %>%
  group_by(n.pop0, low.var, alpha, tau) %>%
  summarise(ll.bar = mean(log.lam),
            ll.var = var(log.lam),
            n = n())

# Plot mean log lambda against tau over time
# note: this low-balls standard error b/c of temporal autocorrelation
ext.tau.summ %>%
  ggplot(aes(x = tau)) +
  geom_line(
    aes(
      y = ll.bar,
      colour = factor(alpha)
    )
  ) +
  geom_ribbon(
    aes(
      ymin = ll.bar - 2 * sqrt(ll.var / n),
      ymax = ll.bar + 2 * sqrt(ll.var / n),
      fill = factor(alpha)
    ),
    alpha = 0.2
  ) +
  facet_wrap(n.pop0 ~ low.var) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple'))

# well... this seems to show that population growth rate doesn't vary by alpha group

# Try plotting by extincton generation...

ext.tau.summ.gen = ext.tau %>%
  group_by(n.pop0, low.var, alpha, trial) %>%
  mutate(ext.gen = max(gen)) %>%
  group_by(n.pop0, low.var, alpha, tau, ext.gen) %>%
  summarise(ll.bar = mean(log.lam),
            ll.var = var(log.lam),
            n = n())

ext.tau.summ.gen %>%
  ggplot(aes(x = tau)) +
  geom_line(
    aes(
      y = ll.bar,
      colour = ext.gen,
      group = ext.gen
    )
  ) +
  geom_ribbon(
    aes(
      ymin = ll.bar - 2 * sqrt(ll.var / n),
      ymax = ll.bar + 2 * sqrt(ll.var / n),
      group = ext.gen,
      fill = ext.gen
    ),
    alpha = 0.2
  ) +
  facet_wrap(paste(n.pop0, low.var, sep = ', ') ~ alpha, nrow = 4) +
  scale_color_gradient(low = 'pink', high = 'firebrick') +
  scale_fill_gradient(low = 'pink', high = 'firebrick')

### Hmm... okay now look at genotypic change

ext.tau.d = ext.mod %>%
  mutate(dbar = 2.75 - gbar) %>%
  group_by(n.pop0, alpha, low.var, trial) %>%
  mutate(deltad = c(exp(diff(log(dbar))), NA)) %>%
  mutate(tau = max(gen) - gen)%>%
  ungroup() %>%
  filter(!is.na(deltad))

ext.tau.d %>%
  ggplot(aes(x = tau, y = deltad)) +
  geom_point(
    alpha = 0.2,
    position = position_jitter(width = 0.2)
  ) +
  stat_smooth() +
  facet_wrap(paste(n.pop0, low.var, sep = ', ') ~ alpha, ncol = 2)

ext.tau.d.summ = ext.tau.d %>%
  group_by(n.pop0, low.var, alpha, tau) %>%
  summarise(dd.bar = mean(deltad),
            dd.var = var(deltad),
            n = n())

ext.tau.d.summ %>%
  ggplot(aes(x = tau)) +
  geom_line(
    aes(
      y = dd.bar,
      colour = factor(alpha)
    )
  ) +
  geom_ribbon(
    aes(
      ymin = dd.bar - 2 * sqrt(dd.var / n),
      ymax = dd.bar + 2 * sqrt(dd.var / n),
      fill = factor(alpha)
    ),
    alpha = 0.2
  ) +
  facet_wrap(n.pop0 ~ low.var) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple'))

# ah... as you get closer to extinction (tau -> 0) rate of adaptation slows down

# Separate everything out by extinction type.

all.deltas = all.data %>%
  filter(extinct) %>%
  arrange(n.pop0, low.var, alpha, trial, desc(gen)) %>%
  distinct(n.pop0, low.var, alpha, trial, .keep_all = TRUE) %>%
  filter(n > 1) %>%
  mutate(gen = gen + 1) %>%
  mutate(n = 1) %>%
  mutate_at(vars(gbar, zbar, wbar, p.fix.pos, p.fix.neg, v),
            list(~ NA)) %>%
  rbind(all.data) %>%
  arrange(n.pop0, low.var, alpha, trial, gen) %>%
  mutate(dbar = 2.75 - gbar) %>%
  group_by(n.pop0, alpha, low.var, trial) %>%
  mutate(log.lam = c(diff(log(n)), NA),
         delta.d = c(exp(diff(log(dbar))), NA)) %>%
  mutate(tau = max(gen) - gen)%>%
  ungroup() %>%
  filter(!is.na(log.lam))

delta.summs = all.deltas %>%
  group_by(n.pop0, low.var, alpha, extinct, tau) %>%
  summarise(ll.bar = mean(log.lam),
            ll.var = var(log.lam),
            dd.bar = mean(delta.d, na.rm = TRUE),
            dd.var = var(delta.d, na.rm = TRUE),
            n = n()) %>%
  ungroup() %>%
  mutate(n.pop0  = factor(n.pop0,  labels = c('small', 'large')),
         low.var = factor(low.var, labels = c('high div.', 'low div.')),
         alpha   = factor(alpha,   labels = c('dens. indp', 'dens. dep')),
         extinct = factor(extinct, labels = c('surviving', 'extinct')))

delta.summs %>%
  ggplot(aes(x = tau)) +
  geom_segment(
    aes(
      x = 1, xend = 14,
      y = 1, yend = 1
    ),
    linetype = 3,
    color = 'gray44'
  ) +
  geom_line(
    aes(
      y = dd.bar,
      colour = factor(alpha),
      linetype = extinct
    )
  ) +
  geom_ribbon(
    aes(
      ymin = dd.bar - 2 * sqrt(dd.var / n),
      ymax = dd.bar + 2 * sqrt(dd.var / n),
      fill = factor(alpha),
      group = interaction(alpha, extinct)
    ),
    alpha = 0.1
  ) +
  facet_wrap(paste(n.pop0, low.var, sep = ', ') ~ .) +
  scale_x_reverse() +
  labs(x = 'Time until extinction (tau)',
       y = 'Proportion of gen. load reduced') +
  scale_color_manual(values = c('black', 'purple'), '') +
  scale_fill_manual(values = c('black', 'purple'), '') +
  scale_linetype_manual(values = 1:2, '') +
  theme(legend.position = 'bottom') +
  ggtitle('Genotypic change (with standard error)') +
  ggsave('~/Dropbox/rescue_ndd_paper_2020/figures/geno_change_tau.pdf')

delta.summs %>%
  ggplot(aes(x = tau)) +
  geom_segment(
    aes(
      x = 1, xend = 14,
      y = 1, yend = 1
    ),
    linetype = 3,
    color = 'gray44'
  ) +
  geom_line(
    aes(
      y = dd.bar,
      colour = factor(alpha),
      linetype = extinct
    )
  ) +
  facet_wrap(paste(n.pop0, low.var, sep = ', ') ~ .) +
  scale_x_reverse() +
  labs(x = 'Time until extinction (tau)',
       y = 'Proportion of gen. load reduced') +
  scale_color_manual(values = c('black', 'purple'), '') +
  scale_linetype_manual(values = 1:2, '') +
  theme(legend.position = 'bottom') +
  ggtitle('Genotypic change (no standard error)') +
  ggsave('~/Dropbox/rescue_ndd_paper_2020/figures/geno_change_tau_no_se.pdf')

lm(delta.d ~ tau, all.deltas %>% filter(n.pop0 %in% 100 & low.var)) %>% summary()
lm(delta.d ~ tau, all.deltas %>% filter(n.pop0 %in% 100 & !low.var)) %>% summary()
lm(delta.d ~ tau, all.deltas %>% filter(n.pop0 %in% 20 & low.var)) %>% summary()
lm(delta.d ~ tau, all.deltas %>% filter(n.pop0 %in% 20 & !low.var)) %>% summary()

lm(delta.d ~ tau * extinct, all.deltas %>% filter(n.pop0 %in% 100 & low.var)) %>% summary()
lm(delta.d ~ tau * extinct, all.deltas %>% filter(n.pop0 %in% 100 & !low.var)) %>% summary()
lm(delta.d ~ tau * extinct, all.deltas %>% filter(n.pop0 %in% 20 & low.var)) %>% summary()
lm(delta.d ~ tau * extinct, all.deltas %>% filter(n.pop0 %in% 20 & !low.var)) %>% summary()

lm(delta.d ~ tau * extinct * factor(n.pop0) * low.var, all.deltas) %>%
  summary()

lme4::lmer(delta.d ~ (1 | trial) + tau * extinct * factor(n.pop0) * low.var, 
           all.deltas) %>%
  summary()
# error for singularity... hmm... random effect here is saying there's a mean effect for each trial
# feel like that shouldn't necessarily be true

# lme4::lmer(delta.d ~ (tau * extinct * factor(n.pop0) * low.var | trial), 
#            all.deltas) %>%
#   summary()

delta.summs %>%
  ggplot(aes(x = tau)) +
  geom_segment(
    aes(
      x = 1, xend = 14,
      y = 0, yend = 0
    ),
    linetype = 3,
    color = 'gray44'
  ) +
  geom_line(
    aes(
      y = ll.bar,
      colour = factor(alpha),
      linetype = extinct
    )
  ) +
  geom_ribbon(
    aes(
      ymin = ll.bar - 2 * sqrt(ll.var / n),
      ymax = ll.bar + 2 * sqrt(ll.var / n),
      fill = factor(alpha),
      group = interaction(alpha, extinct)
    ),
    alpha = 0.1
  ) +
  facet_wrap(paste(n.pop0, low.var, sep = ', ') ~ .) +
  scale_x_reverse() +
  labs(x = 'Time until extinction (tau)',
       y = 'Log population growth rate') +
  scale_color_manual(values = c('black', 'purple'), '') +
  scale_fill_manual(values = c('black', 'purple'), '') +
  scale_linetype_manual(values = 1:2, '') +
  theme(legend.position = 'bottom') +
  ggtitle('Population size change') +
  ggsave('~/Dropbox/rescue_ndd_paper_2020/figures/log_lambda_tau.pdf')

