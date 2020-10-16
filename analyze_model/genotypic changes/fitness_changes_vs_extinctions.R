# extinctons

extns = all.data %>%
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
  filter(gen < 15)

# changes in fitness

wbars = all.data %>%
  # Get number of extant trials in each generation for each treatment
  group_by(n.pop0, low.var, alpha, trial) %>%
  # get changes in w
  mutate(dw = c(diff(wbar), 0)) %>%
  # now get mean delta w
  group_by(n.pop0, low.var, alpha, gen) %>%
  summarise(dwbar = mean(dw),
            dwvar = var(dw),
            n = n()) %>%
  filter(gen < 15)

# For each generation (w/in var x size treatments),
# which dwbar is larger?
comp.dw = wbars %>%
  select(-c(n, dwvar)) %>%
  ungroup() %>%
  mutate(alpha = factor(alpha, labels = c('di', 'ndd'))) %>%
  spread(key = alpha, value = dwbar) %>%
  mutate(ndd.gt = ndd > di)

all.wbars = merge(extns, wbars, 
                  by = c('n.pop0', 'low.var', 'alpha', 'gen'), 
                  all.x = TRUE, all.y = TRUE) %>%
  # Now, merge with comparison data structure above
  merge(y = comp.dw, by = c('n.pop0', 'low.var', 'gen')) %>%
  arrange(n.pop0, low.var, alpha, gen) #%>%
  # Change factor labels
  mutate(n0 = factor(n.pop0, labels = c("Small", "Large")),
         alpha = factor(alpha, labels = c("Density dependent", "Density independent")),
         low.var = factor(low.var, labels = c("High variance", "Low variance")))

nrow(all.wbars)
head(all.wbars)

all.wbars %>%
  ggplot() +
  geom_line(
    aes(
      x = prop.extant,
      y = dwbar,
      group = interaction(n.pop0, low.var, gen),
      colour = ndd.gt
    ),
    size = 2
  ) +
  geom_point(
    aes(
      x = prop.extant,
      y = dwbar,
      fill = factor(alpha)
    ),
    shape = 21,
    size = 2
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple'))

all.wbars %>%
  ungroup() %>%
  summarise(mean(ndd.gt))

all.wbars %>%
  group_by(n.pop0, low.var) %>%
  summarise(ndd.gt = mean(ndd.gt),
            n = n())

all.wbars %>%
  ggplot() +
  geom_line(
    aes(
      x = prop.extant,
      y = dwbar,
      group = interaction(n.pop0, low.var, gen),
      colour = ndd.gt
    ),
    size = 1
  ) +
  geom_point(
    aes(
      x = prop.extant,
      y = dwbar,
      fill = factor(alpha)
    ),
    shape = 21,
    size = 3
  ) +
  scale_color_manual(values = c('black', 'plum')) +
  scale_fill_manual(values = c('black', 'plum'), c('DI', 'NDD')) +
  guides(color = 'none') +
  facet_wrap(n0 ~ low.var) +
  labs(x = "Prob. of extinction", 
       y = expression(paste('Change in mean intrinsic fitness, ', Delta, W[t]))) +
  theme(panel.background = element_blank(),
        legend.position = 'bottom')


all.wbars %>%
  ggplot() +
  geom_line(
    aes(
      x = prop.extant,
      y = dwbar,
      group = interaction(n.pop0, low.var, gen),
      colour = ndd.gt
    ),
    size = 1
  ) +
  geom_point(
    aes(
      x = prop.extant,
      y = dwbar,
      fill = gen,
      shape = factor(alpha)
    ),
    size = 3
  ) +
  scale_fill_continuous(high = 'yellow', low = 'forestgreen') +
  scale_colour_manual(values = c('red', 'blue')) +
  scale_shape_manual(values = c(21, 24)) +
  theme(panel.background = element_blank())

diff.bars = all.wbars %>%
  group_by(n.pop0, low.var, gen) %>%
  summarise(d_dw = dwbar[alpha == 0] - dwbar[alpha > 0],
            p_di = prop.extant[alpha == 0],
            p_dd = prop.extant[alpha >  0])

diff.bars %>%
  ggplot() +
  geom_segment(
    aes(
      x = p_di, xend = p_dd,
      y = d_dw, yend = d_dw
    )
  ) +
  geom_point(
    aes(
      x = p_di, 
      y = d_dw
    ),
    colour = 'black',
    size = 2
  ) +
  geom_point(
    aes(
      x = p_dd,
      y = d_dw
    ),
    colour = 'plum',
    size = 2
  ) +
  facet_wrap(n.pop0 ~ low.var)

