# Plot results from 24 source populations

library(ggplot2)
library(dplyr)
library(tidyr)

expmt.sum = read.csv('experiment/expmt_bottleneck24/out/sim_results.csv')

# Plot with treatments in panels
experiment.summ %>%
  group_by(trial) %>%
  mutate(extinct = max(gen) < 15) %>%
  ggplot(aes(x = gen, y = n)) +
  geom_line(
    aes(
      group = trial,
      linetype = evolve,
      size = evolve,
      colour = extinct
    )
  ) +
  scale_color_manual(values = c('black', 'red')) +
  scale_linetype_manual(values = 2:1) +
  scale_size_manual(values = c(0.25, 0.75)) +
  scale_y_log10() +
  facet_wrap(n.pop0 ~ alpha)

# Plot with source population in panel
experiment.summ %>%
  ggplot(aes(x = gen, y = n)) +
  geom_line(
    aes(
      group = trial,
      colour = factor(alpha > 0),
      linetype = evolve
    )
  ) +
  geom_point(
    aes(
      colour = factor(alpha > 0),
      fill = factor(alpha > 0),
      shape = factor(n.pop0)
    )
  ) +
  scale_fill_manual(values = c('black', 'purple')) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_shape_manual(values = c(1, 21)) +
  scale_y_log10() +
  facet_wrap(~ sourcepops, nrow = 4) +
  theme(legend.position = 'none')
