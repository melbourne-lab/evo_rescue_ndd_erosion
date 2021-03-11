library(ggplot2)
library(dplyr)
library(tidyr)

all.data = 
  rbind(
    read.csv('simulations/outputs/longsims/longsims_n100_a000_hivar.csv'),
    read.csv('simulations/outputs/longsims/longsims_n100_a000_lowvar.csv'),
    read.csv('simulations/outputs/longsims/longsims_n20_a000_hivar.csv'),
    read.csv('simulations/outputs/longsims/longsims_n20_a000_lowvar.csv'),
    read.csv('simulations/outputs/longsims/longsims_n100_a035_hivar.csv'),
    read.csv('simulations/outputs/longsims/longsims_n100_a035_lowvar.csv'),
    read.csv('simulations/outputs/longsims/longsims_n20_a035_hivar.csv'),
    read.csv('simulations/outputs/longsims/longsims_n20_a035_lowvar.csv')
  )

# write.csv(all.data, row.names = FALSE,
#           'simulations/outputs/final_results/longsims/longsims_combined.csv')

nrow(all.data)

hmm = all.data %>%
  arrange(trial, n.pop0, low.var, alpha, desc(gen)) %>%
  distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE) %>%
  filter(gen < 50)

head(hmm)

hmm %>% ggplot(aes(x = n)) + geom_histogram() + facet_wrap(n.pop0 + low.var ~ alpha, nrow = 2)

all.data = all.data %>%
  arrange(trial, n.pop0, low.var, alpha, desc(gen)) %>%
  distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE) %>%
  group_by(n.pop0, low.var, alpha, trial) %>%
  summarise(extinct = gen < 50 & n < 10000) %>%
  merge(x = all.data, y = ., by = c('trial', 'n.pop0', 'low.var', 'alpha'))

# Now do some stuff

# Extinction probability:

all.data %>%
  distinct(n.pop0, low.var, alpha, trial, .keep_all = TRUE) %>%
  group_by(n.pop0, low.var, alpha) %>%
  summarise(n = mean(extinct))

all.ext = all.data %>%
  group_by(n.pop0, low.var, alpha, trial) %>%
  summarise(ext.gen = ifelse(extinct[1], max(gen), 50))

all.ext %>%
  filter(ext.gen < 50) %>%
  ggplot() +
  geom_histogram(aes(x = ext.gen, fill = factor(alpha)), 
                 position = 'identity',
                 binwidth = 1, alpha = 0.55) +
  scale_fill_manual(values = c('black', 'purple')) +
  facet_wrap(n.pop0 ~ low.var)

all.ext %>%
  filter(ext.gen < 50) %>%
  group_by(n.pop0, low.var, alpha) %>%
  summarise(tau = mean(ext.gen),
            tau.se = sd(ext.gen / n()))

all.ext %>%
  filter(ext.gen < 50) %>%
  group_by(n.pop0, low.var, alpha) %>%
  summarise(tau = mean(ext.gen)) %>%
  spread(key = alpha, value = tau)

all.ext %>% filter(ext.gen < 50) %>%
  glm(ext.gen ~ factor(n.pop0) * low.var * factor(alpha), data = ., family = 'poisson') %>%
  summary()

all.ext %>% filter(ext.gen < 50) %>%
  glm(ext.gen ~ factor(n.pop0) + low.var + 
                factor(n.pop0) * factor(alpha) + 
                low.var * factor(alpha), 
      data = ., family = 'poisson') %>%
  summary()

all.ext %>% filter(ext.gen < 50) %>%
  glm(ext.gen ~ factor(n.pop0) + low.var + factor(n.pop0) * factor(alpha), 
      data = ., family = 'poisson') %>%
  summary()
  
# Feel like a censored approach here would be helpful!

### Plot accumulation of extinctions
all.extinctions.long = all.data %>%
  select(gen, trial, n.pop0, low.var, alpha, extinct) %>%
  group_by(trial, n.pop0, low.var, alpha) %>%
  mutate(ext.this.gen = (gen %in% max(gen)) & extinct) %>%
  group_by(n.pop0, low.var, alpha, gen) %>%
  summarise(n.extinctions = sum(ext.this.gen)) %>%
  ungroup() %>%
  arrange(n.pop0, low.var, alpha, gen)

all.extinctions.long %>%
  group_by(n.pop0, low.var, alpha) %>%
  mutate(p.extinct = cumsum(n.extinctions / 1000)) %>%
  ggplot(aes(x = gen)) +
  geom_ribbon(
    aes(
      xmin = 0, xmax = 50,
      ymin = 0, ymax = p.extinct,
      group = interaction(n.pop0, low.var, alpha),
      fill = factor(alpha)
    ),
    alpha = 0.25
  ) +
  scale_fill_manual(values = c('black', 'purple')) +
  labs(x = 'Generation', y = 'Probability of extinction') +
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 1) +
  theme(legend.position = 'none',
        panel.background = element_blank(),
        strip.background = element_rect(colour = 'black'))

all.extinctions.long %>%
  group_by(n.pop0, low.var, alpha) %>%
  mutate(n.extant = 1000 - c(0, cumsum(n.extinctions)[-50]),
         p.instant.extinct = n.extinctions / n.extant) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = p.instant.extinct,
      group = factor(alpha),
      colour = factor(alpha)
    ),
    size = 3
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 1) +
  theme(legend.position = 'none',
        panel.background = element_blank(),
        strip.background = element_rect(colour = 'black'))

all.data %>%
  distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE) %>%
  glm(extinct ~ factor(n.pop0) * low.var * factor(alpha) * gbar,
      data = ., family = 'binomial') %>%
  summary()

all.data %>%
  distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE) %>%
  glm(extinct ~ factor(n.pop0) * low.var * factor(alpha) +
                factor(n.pop0) * low.var * gbar +
                factor(n.pop0) * gbar * factor(alpha) +
                gbar * low.var * factor(alpha),
      data = ., family = 'binomial') %>%
  summary()

all.data %>%
  distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE) %>%
  glm(extinct ~ factor(n.pop0) * low.var * factor(alpha) +
                factor(n.pop0) * low.var * gbar +
                factor(n.pop0) * gbar * factor(alpha),
      data = ., family = 'binomial') %>%
  summary()

all.data %>%
  distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE) %>%
  glm(extinct ~ factor(n.pop0) * low.var * factor(alpha) +
                factor(n.pop0) * gbar * factor(alpha) +
                low.var * gbar,
      data = ., family = 'binomial') %>%
  summary()

all.data %>%
  distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE) %>%
  glm(extinct ~ factor(n.pop0) * low.var * factor(alpha) +
                factor(n.pop0) * gbar +
                gbar * factor(alpha) +
                low.var * gbar,
      data = ., family = 'binomial') %>%
  summary()

all.data %>%
  distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE) %>%
  glm(extinct ~ factor(n.pop0) * low.var * factor(alpha) +
                gbar * factor(alpha) +
                low.var * gbar,
      data = ., family = 'binomial') %>%
  summary()

all.data %>%
  distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE) %>%
  glm(extinct ~ factor(n.pop0) * low.var * factor(alpha) +
                low.var * gbar,
      data = ., family = 'binomial') %>%
  summary()

# kinda puny model!

ext.mod.final =  glm(
  extinct ~ factor(n.pop0) * low.var * factor(alpha) +
            low.var * gbar,
  data = all.data %>% distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE), 
  family = 'binomial'
  ) 

expand.grid(gbar = (-10:10)/10,
            alpha = c(0, 0.0035),
            n.pop0 = factor(c(20, 100)),
            low.var = c(TRUE, FALSE)) %>%
  (function(x) cbind(x, p = predict(ext.mod.final, x, type = 'response'))) %>%
  ggplot() +
  geom_line(
    aes(
      x = gbar, y = p, group = factor(alpha),
      colour = factor(alpha)
    )
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 1)


all.data %>%
  ggplot() +
  geom_line(aes(x = gen, y = n, colour = extinct, group = trial),
            size = 0.2) +
  scale_y_log10() +
  scale_color_manual(values = c('black', 'red')) +
  facet_wrap(n.pop0 + low.var ~ alpha, ncol = 4)
  
### Old stuff - will almost surely delete

  # Get number of extant trials in each generation for each treatment
  group_by(n.pop0, low.var, alpha, gen) %>%
  summarise(n.trials = n()) %>%
  # Summarise extinctions
  mutate(# Number of trials going extinct
    d.trials = -c(0, diff(n.trials)),
    # Proportion of extant trials going extinct
    prop.extant = c(d.trials[-1], 0) / n.trials,
    # Proportion of all trials going extinct
    prop.totaln = c(d.trials[-1], 0) / max(n.trials),
    # Proportion of extinct trials extinct in this generation
    prop.condit = c(d.trials[-1], 0) / sum(d.trials),
    # Proportion of trials surviving
    prop.surviv = n.trials / max(n.trials)) %>%
  ungroup() %>%
  mutate(n.pop0 = factor(n.pop0, labels = c("Small", "Large")),
         alpha = factor(alpha, labels = c("Density independent", "Density dependent")),
         low.var = factor(low.var, labels = c("Low diversity", "High diversity")))


cuml.probs = ggplot(all.ext %>% filter(gen < 50), aes(x = gen)) +
  geom_ribbon(
    aes(
      xmin = 0, xmax = 14,
      ymin = 0, ymax = 1 - prop.surviv,
      group = interaction(n.pop0, low.var, alpha),
      fill = factor(alpha)
    ),
    alpha = 0.25
  ) +
  scale_fill_manual(values = c('black', 'purple')) +
  labs(x = 'Generation', y = 'Probability of extinction') +
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 1) +
  theme(legend.position = 'none',
        panel.background = element_blank(),
        strip.background = element_rect(colour = 'black'))

cuml.probs
