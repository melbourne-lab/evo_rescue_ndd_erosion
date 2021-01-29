library(ggplot2)
library(dplyr)
library(tidyr)

yes = data.frame(
  x = rep(0:9, each = 2) + rep(1:2, times = 10),
  y = rep(0:4, each = 4) + rep(1:2, times = 10)
)

plot(yes, type = 'l')

yes %>%
  ggplot(aes(x = x, y = y)) +
  geom_line()

# Try for rizzle

all.data = read.csv('simulations/outputs/final_results/alldata_combined.csv')

# Want to plot two populatons.

test.two = all.data %>% 
  filter(trial %in% 1:2, !alpha, n.pop0 < 50, low.var) %>%
  rbind(., .) %>%
  arrange(trial, gen) %>%
  group_by(trial) %>%
  mutate(v2 = c(NA, v[1:(n()-1)])) %>%
  ungroup() %>%
  filter(!is.na(v2))
  #ungroup()

head(test.two)

test.two %>%
  ggplot() +
  geom_path(aes(x = n, y = v2, group = trial)) +
  geom_point(aes(x = n, y = v2, colour = factor(gen), shape = factor(trial)), size = 3) +
  scale_color_brewer('BuPu')

###

steps.all = all.data %>% 
  rbind(., .) %>%
  arrange(trial, alpha, low.var, n.pop0, gen) %>%
  group_by(trial, alpha, low.var, n.pop0) %>%
  mutate(v2 = c(NA, v[1:(n()-1)])) %>%
  ungroup() %>%
  filter(!is.na(v2))

steps.all %>% filter(trial %in% 1:2, !alpha, n.pop0 < 50, low.var)

steps.all %>%
  filter(trial %in% 1:100) %>%
  ggplot(aes(x = n, y = v2)) +
  geom_path(aes(colour = extinct, group = trial), alpha = 0.5) +
  scale_x_log10() +
  scale_color_manual(values = c('black', 'red')) +
  facet_wrap(n.pop0 + low.var ~ alpha, nrow = 2)

### Try with summary

all.summ = all.data %>%
  group_by(trial, alpha, low.var, n.pop0) %>%
  mutate(tau = max(gen) - gen + 1 - as.numeric(any(n < 2) & extinct)) %>%
  group_by(alpha, low.var, n.pop0, extinct, tau) %>%
  summarise(n.bar = mean(n),
            n.var = var(n),
            v.bar = mean(v),
            v.var = var(v),
            w.bar = mean(wbar),
            w.var = var(wbar),
            nobs = n()
  ) %>%
  ungroup()

all.summ %>%
  ggplot(aes(x = n.bar, y = v.bar)) +
  geom_path(
    aes(
      group = interaction(alpha, extinct), 
      linetype = extinct, 
      colour = factor(alpha)
      )
    ) +
  geom_point(
    aes(
      colour = factor(alpha)
    )
  ) +
  # geom_ribbon(
  #   aes(
  #     group = interaction(alpha, extinct), 
  #     linetype = extinct, 
  #     colour = factor(alpha)
  #   )
  # ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_x_log10() +
  facet_wrap(n.pop0 ~ low.var, ncol = 4)

all.summ %>%
  ggplot(aes(x = n.bar, y = w.bar)) +
  geom_path(
    aes(
      group = interaction(alpha, extinct), 
      linetype = extinct, 
      colour = factor(alpha)
    )
  ) +
  geom_point(
    aes(
      colour = factor(alpha)
    )
  ) +
  # geom_ribbon(
  #   aes(
  #     group = interaction(alpha, extinct), 
  #     linetype = extinct, 
  #     colour = factor(alpha)
  #   )
  # ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_x_log10() +
  facet_wrap(n.pop0 ~ low.var, ncol = 4)

all.summ %>%
  ggplot(aes(x = v.bar, y = w.bar)) +
  geom_path(
    aes(
      group = interaction(alpha, extinct), 
      linetype = extinct, 
      colour = factor(alpha)
    )
  ) +
  geom_point(
    aes(
      colour = factor(alpha)
    )
  ) +
  # geom_ribbon(
  #   aes(
  #     group = interaction(alpha, extinct), 
  #     linetype = extinct, 
  #     colour = factor(alpha)
  #   )
  # ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_x_log10() +
  facet_wrap(n.pop0 ~ low.var, ncol = 4)

### Trying one plot with confints...

all.summ %>%
  filter(tau > 0) %>%
  ggplot() +
  geom_path(
    aes(
      x = v.bar,
      y = w.bar,
      group = interaction(alpha, extinct), 
      linetype = extinct, 
      colour = factor(alpha)
    )
  ) +
  geom_path(
    aes(
      x = v.bar - 2 * sqrt(v.var / sqrt(nobs)),
      y = w.bar - 2 * sqrt(w.var / sqrt(nobs)),
      group = interaction(alpha, extinct),
      colour = factor(alpha)
    ),
    linetype = 3,
    size = 0.4
  ) +
  geom_path(
    aes(
      x = v.bar + 2 * sqrt(v.var / sqrt(nobs)),
      y = w.bar + 2 * sqrt(w.var / sqrt(nobs)),
      group = interaction(alpha, extinct),
      colour = factor(alpha)
    ),
    linetype = 3,
    size = 0.4
  ) +
  geom_point( 
    aes(
        x = v.bar,
        y = w.bar,
        colour = factor(alpha)
    )
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  facet_wrap(n.pop0 ~ low.var, ncol = 4)

### Try putting everything together.

library(cowplot)

plot.nv = all.summ %>%
  ggplot(aes(x = n.bar, y = v.bar)) +
  geom_path(
    aes(
      group = interaction(alpha, extinct), 
      linetype = extinct, 
      colour = factor(alpha)
    )
  ) +
  geom_point(
    aes(
      colour = factor(alpha)
    )
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_x_log10() +
  facet_wrap(low.var ~ n.pop0, ncol = 4) +
  theme(legend.position = 'none')

plot.nw = all.summ %>%
  ggplot(aes(x = n.bar, y = w.bar)) +
  geom_path(
    aes(
      group = interaction(alpha, extinct), 
      linetype = extinct, 
      colour = factor(alpha)
    )
  ) +
  geom_point(
    aes(
      colour = factor(alpha)
    )
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_x_log10() +
  facet_wrap(low.var ~ n.pop0, ncol = 4) +
  theme(legend.position = 'none')

plot.vw = all.summ %>%
  ggplot(aes(x = v.bar, y = w.bar)) +
  geom_path(
    aes(
      group = interaction(alpha, extinct), 
      linetype = extinct, 
      colour = factor(alpha)
    )
  ) +
  geom_point(
    aes(
      colour = factor(alpha)
    )
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  facet_wrap(low.var ~ n.pop0, ncol = 4) +
  theme(legend.position = 'none')

plot_grid(plot.nv, plot.nw, plot.vw, 
          nrow = 3, rel_heights = c(1, 1, 1))

