library(ggplot2)
library(dplyr)
library(tidyr)

all.data = read.csv('simulations/outputs/final_results/alldata_combined.csv')

# tau.n.v = all.data %>%
#   group_by(trial, n.pop0, low.var, alpha) %>%
#   mutate(tau = max(gen) - gen + 1 - as.numeric(any(n < 2))) %>%
#   group_by(tau, n.pop0, low.var, alpha, extinct) %>%
#   summarise(n = mean(n),
#             v = mean(v))

ext.n.v = all.data %>%
  group_by(trial, n.pop0, low.var, alpha) %>%
  mutate(ext.gen = max(gen)) %>%
  group_by(n.pop0, low.var, alpha, gen, ext.gen) %>%
  summarise(n = mean(n),
            v = mean(v),
            extinct = ext.gen[1] < 15)

# Plotting Pr(x | ext) in just in 2d
ext.n.v %>%
  filter(extinct) %>%
  ggplot() +
  geom_line(
    aes(
      x = v,
      y = n,
      group = interaction(alpha, ext.gen),
      colour = factor(alpha)
    )
  ) +
  geom_point(
    aes(
      x = v,
      y = n,
      colour = factor(alpha)
    )
  ) +
  scale_y_log10() +
  facet_wrap(n.pop0 ~ low.var)

ext.n.v %>%
  rbind(., .) %>%
  arrange(gen, n.pop0, low.var, alpha)

# Points - extincton for population size, variance

all.data %>%
  filter(trial < 1000) %>%
  ggplot(aes(x = v, y = n, colour = extinct)) +
  geom_point(alpha = 0.2, position = position_jitter(height = 0.5, width = 0.01)) +
  scale_color_manual(values = c('black', 'red')) +
  scale_y_log10() +
  facet_wrap(paste0(n.pop0, low.var) ~ alpha, nrow = 4)

# Try the binning exercise...

all.bins = all.data %>%
  mutate(nr = ifelse(n > 100, 101, n),
         vr = round(v * 50) / 50) %>%
  group_by(n.pop0, low.var, alpha, nr, vr) %>%
  summarise(p.extinct = mean(extinct))

# Viridis, small boxes
all.bins %>%
  ggplot() +
  geom_tile(
    aes(
      x = vr, y = nr,
      fill = p.extinct
    )
  ) +
  scale_fill_viridis_b(option = 'B') +
  facet_wrap(low.var ~ paste0(n.pop0, alpha) , nrow = 2) +
  theme_bw()

# Pr(extinct) > 0.5
all.bins %>%
  ggplot() +
  geom_tile(
    aes(
      x = vr, y = nr,
      fill = p.extinct > 0.50
    )
  ) +
  scale_fill_manual(values = c('black', 'red')) +
  facet_wrap(low.var ~ paste0(n.pop0, alpha, sep = ', ') , nrow = 2) +
  theme_bw()

# Try to go out wide
bins.spread = all.bins %>%
  ungroup() %>%
  mutate(alpha = ifelse(alpha > 0, 'dep', 'indep')) %>%
  spread(key = alpha, value = p.extinct)

fill.leg = c(
  # neither is at extinction risk
  FALSEFALSE = 'black',
  # indep. not at risk, dep. is
  FALSETRUE  = 'orange',
  # indep. at risk, dep. is not
  TRUEFALSE  = 'pink',
  # both are at risk
  TRUETRUE   = 'firebrick',
  # missing data
  TRUENA     = 'salmon',
  NATRUE     = 'salmon',
  FALSENA    = 'gray55',
  NAFALSE    = 'gray55'
)

bins.spread %>%
  mutate(fill.tab = paste0(indep > 0.5, dep > 0.5)) %>%
  ggplot() +
  geom_tile(
    aes(
      x = vr, y = nr,
      fill = fill.tab
    )
  ) +
  #scale_fill_manual(values = c('black', 'red')) +
  scale_fill_manual(values = fill.leg) +
  facet_wrap(low.var ~ n.pop0 , nrow = 2) +
  theme_bw()

### Bins 2:

all.bins2 = all.data %>%
  mutate(nr = ifelse(n > 100, 101, round(n / 2) * 2),
         vr = round(v * 50) / 50) %>%
  group_by(n.pop0, low.var, alpha, nr, vr) %>%
  summarise(p.extinct = mean(extinct),
            n.obs = n()) %>%
  uncount(weights = 2) %>%
  group_by(nr, vr, n.pop0, low.var, alpha) %>%
  mutate(n = nr + 0:1)

all.bins2 %>%
  ggplot() +
  geom_tile(
    aes(
      x = vr, y = n,
      fill = p.extinct > 0.50
    )
  ) +
  scale_fill_manual(values = c('black', 'red')) +
  facet_wrap(low.var ~ paste0(n.pop0, alpha) , nrow = 2) +
  theme_bw()

bins.spread2 = all.bins2 %>%
  ungroup() %>%
  filter(n.obs > 4) %>%
  select(-n.obs) %>%
  mutate(alpha = ifelse(alpha > 0, 'dep', 'indep')) %>%
  spread(key = alpha, value = p.extinct)

nrow(bins.spread2)

bins.spread2 %>%
  filter(!(is.na(dep) | is.na(indep))) %>%
  mutate(fill.tab = paste0(indep > 0.5, dep > 0.5)) %>%
  ggplot() +
  geom_tile(
    aes(
      x = vr, y = n,
      fill = fill.tab
    )
  ) +
  #scale_fill_manual(values = c('black', 'red')) +
  scale_fill_manual(values = fill.leg) +
  facet_wrap(low.var ~ n.pop0 , nrow = 2) +
  theme_bw()

bins.spread2 %>%
  filter(!(is.na(dep) | is.na(indep))) %>%
  mutate(fill.ratio = ifelse(indep > 0, dep / indep,
                             ifelse(dep > 0, NA, 1))) %>%
  ggplot() +
  geom_tile(
    aes(
      x = vr, y = n,
      fill = log(fill.ratio)
    )
  ) +
  #scale_fill_manual(values = c('black', 'red')) +
  #scale_fill_manual(values = fill.leg) +
  scale_fill_gradient2(low = 'white', mid = 'black', high = 'firebrick',
                       midpoint = 0) +
  facet_wrap(low.var ~ n.pop0 , nrow = 2) +
  theme_bw()
