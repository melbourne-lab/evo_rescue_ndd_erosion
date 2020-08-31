library(ggplot2)
library(dplyr)
library(tidyr)

all.d = rbind(
  read.csv('simulations/outputs/thrice_bottlenecked_results/a000_hivar_type_all.csv')  %>%
    mutate(gse = sqrt(gvar / n),
           zse = sqrt(zvar / n),
           wse = sqrt(wvar / n)) %>%
    mutate(bottleneck = FALSE, ndd = FALSE),
  read.csv('simulations/outputs/thrice_bottlenecked_results/a035_hivar_type_all.csv') %>%
    mutate(gse = sqrt(gvar / n),
           zse = sqrt(zvar / n),
           wse = sqrt(wvar / n)) %>%
    mutate(bottleneck = FALSE, ndd = TRUE),
  read.csv('simulations/outputs/thrice_bottlenecked_results/a000_lovar_type_all.csv') %>%
    mutate(gvar = gvar_win + gvar_amg,
           zvar = zvar_win + zvar_amg,
           wvar = wvar_win + wvar_amg) %>%
    mutate(gse = sqrt(gvar_win / n + gvar_amg / 27),
           zse = sqrt(zvar_win / n + zvar_amg / 27),
           wse = sqrt(wvar_win / n + wvar_amg / 27)) %>%
    select(n.pop0, gen, gbar, gvar, zbar, zvar, wbar, wvar, n, gse, zse, wse) %>%
    mutate(bottleneck = TRUE, ndd = FALSE),
  read.csv('simulations/outputs/thrice_bottlenecked_results/a035_lovar_type_all.csv') %>%
    mutate(gvar = gvar_win + gvar_amg,
           zvar = zvar_win + zvar_amg,
           wvar = wvar_win + wvar_amg) %>%
    mutate(gse = sqrt(gvar_win / n + gvar_amg / 27),
           zse = sqrt(zvar_win / n + zvar_amg / 27),
           wse = sqrt(wvar_win / n + wvar_amg / 27)) %>%
    select(n.pop0, gen, gbar, gvar, zbar, zvar, wbar, wvar, n, gse, zse, wse) %>%
    mutate(bottleneck = TRUE, ndd = TRUE)
)

unique(all.d$n.pop0) # phew!

all.d %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = wbar,
      group = interaction(n.pop0, bottleneck, ndd),
      colour = ndd
    )
  ) +
  geom_ribbon(
    aes(
      ymin = wbar - 2 * wse,
      ymax = wbar + 2 * wse,
      group = interaction(n.pop0, bottleneck, ndd),
      fill = ndd
    ),
    alpha = 0.1
  ) +
  facet_wrap(n.pop0 ~ bottleneck) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  theme(panel.background = element_blank())

### Old plots.

all.d %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = wbar,
      group = interaction(n.pop0, bottleneck, ndd),
      colour = ndd
    )
  ) +
  geom_ribbon(
    aes(
      ymin = wbar - 2 * sqrt(wvar / n),
      ymax = wbar + 2 * sqrt(wvar / n),
      group = interaction(n.pop0, bottleneck, ndd),
      fill = ndd
    ),
    alpha = 0.1
  ) +
  facet_wrap(n.pop0 ~ bottleneck) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  theme(panel.background = element_blank())

###

ext.d = rbind(
  read.csv('simulations/outputs/thrice_bottlenecked_results/a000_hivar_type_ext.csv')  %>%
    mutate(gse = sqrt(gvar / n),
           zse = sqrt(zvar / n),
           wse = sqrt(wvar / n)) %>%
    mutate(bottleneck = FALSE, ndd = FALSE),
  read.csv('simulations/outputs/thrice_bottlenecked_results/a035_hivar_type_ext.csv') %>%
    mutate(gse = sqrt(gvar / n),
           zse = sqrt(zvar / n),
           wse = sqrt(wvar / n)) %>%
    mutate(bottleneck = FALSE, ndd = TRUE),
  read.csv('simulations/outputs/thrice_bottlenecked_results/a000_lovar_type_ext.csv') %>%
    mutate(gvar = gvar_win + gvar_amg,
           zvar = zvar_win + zvar_amg,
           wvar = wvar_win + wvar_amg) %>%
    mutate(gse = sqrt(gvar_win / n + gvar_amg / 27),
           zse = sqrt(zvar_win / n + zvar_amg / 27),
           wse = sqrt(wvar_win / n + wvar_amg / 27)) %>%
    select(n.pop0, extinct, gen, gbar, gvar, zbar, zvar, wbar, wvar, n, gse, zse, wse) %>%
    mutate(bottleneck = TRUE, ndd = FALSE),
  read.csv('simulations/outputs/thrice_bottlenecked_results/a035_lovar_type_ext.csv') %>%
    mutate(gvar = gvar_win + gvar_amg,
           zvar = zvar_win + zvar_amg,
           wvar = wvar_win + wvar_amg) %>%
    mutate(gse = sqrt(gvar_win / n + gvar_amg / 27),
           zse = sqrt(zvar_win / n + zvar_amg / 27),
           wse = sqrt(wvar_win / n + wvar_amg / 27)) %>%
    select(n.pop0, extinct, gen, gbar, gvar, zbar, zvar, wbar, wvar, n, gse, zse, wse) %>%
    mutate(bottleneck = TRUE, ndd = TRUE)
)

ext.d %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = wbar,
      group = interaction(n.pop0, bottleneck, ndd, extinct),
      colour = ndd
    )
  ) +
  geom_ribbon(
    aes(
      ymin = wbar - 2 * wse,
      ymax = wbar + 2 * wse,
      group = interaction(n.pop0, bottleneck, ndd, extinct),
      fill = ndd
    ),
    alpha = 0.1
  ) +
  facet_wrap(ndd ~ interaction(n.pop0, bottleneck), ncol = 4) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  theme(panel.background = element_blank())

ext.d %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = wbar,
      group = interaction(n.pop0, bottleneck, ndd, extinct),
      linetype = extinct,
      colour = ndd
    )
  ) +
  geom_ribbon(
    aes(
      ymin = wbar - 2*wse,
      ymax = wbar + 2*wse,
      group = interaction(n.pop0, bottleneck, ndd, extinct),
      fill = ndd
    ),
    alpha = 0.1
  ) +
  facet_wrap(n.pop0 ~ bottleneck) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  theme(panel.background = element_blank())

# old versions

ext.d %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = wbar,
      group = interaction(n.pop0, bottleneck, ndd, extinct),
      colour = ndd
    )
  ) +
  geom_ribbon(
    aes(
      ymin = wbar - 2 * sqrt(wvar / n),
      ymax = wbar + 2 * sqrt(wvar / n),
      group = interaction(n.pop0, bottleneck, ndd, extinct),
      fill = ndd
    ),
    alpha = 0.1
  ) +
  facet_wrap(ndd ~ interaction(n.pop0, bottleneck), ncol = 4) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  theme(panel.background = element_blank())

ext.d %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = wbar,
      group = interaction(n.pop0, bottleneck, ndd, extinct),
      linetype = extinct,
      colour = ndd
    )
  ) +
  geom_ribbon(
    aes(
      ymin = wbar - 2 * sqrt(wvar / n),
      ymax = wbar + 2 * sqrt(wvar / n),
      group = interaction(n.pop0, bottleneck, ndd, extinct),
      fill = ndd
    ),
    alpha = 0.1
  ) +
  facet_wrap(n.pop0 ~ bottleneck) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  theme(panel.background = element_blank())

###

gen.d = rbind(
  read.csv('simulations/outputs/thrice_bottlenecked_results/a000_hivar_type_gen.csv') %>%
    mutate(bottleneck = FALSE, ndd = FALSE),
  read.csv('simulations/outputs/thrice_bottlenecked_results/a035_hivar_type_gen.csv') %>%
    mutate(bottleneck = FALSE, ndd = TRUE),
  read.csv('simulations/outputs/thrice_bottlenecked_results/a000_lovar_type_gen.csv') %>%
    mutate(gvar = gvar_win + gvar_amg,
           zvar = zvar_win + zvar_amg,
           wvar = wvar_win + wvar_amg) %>%
    select(n.pop0, ext.gen, gen, gbar, gvar, zbar, zvar, wbar, wvar, n) %>%
    mutate(bottleneck = TRUE, ndd = FALSE),
  read.csv('simulations/outputs/thrice_bottlenecked_results/a035_lovar_type_gen.csv') %>%
    mutate(gvar = gvar_win + gvar_amg,
           zvar = zvar_win + zvar_amg,
           wvar = wvar_win + wvar_amg) %>%
    select(n.pop0, ext.gen, gen, gbar, gvar, zbar, zvar, wbar, wvar, n) %>%
    mutate(bottleneck = TRUE, ndd = TRUE)
)

gen.d %>%
  ggplot(aes(x = gen)) +
  geom_line(
    data = . %>% filter(ext.gen < 15),
    aes(
      y = wbar,
      group = interaction(n.pop0, bottleneck, ndd, ext.gen),
      colour = ext.gen
    )
  ) +
  geom_ribbon(
    data = . %>% filter(ext.gen < 15),
    aes(
      ymin = wbar - 2 * sqrt(wvar / n),
      ymax = wbar + 2 * sqrt(wvar / n),
      group = interaction(n.pop0, bottleneck, ndd, ext.gen),
      fill = ext.gen
    ),
    alpha = 0.1
  ) +
  geom_line(
    data = . %>% filter(ext.gen %in% 15),
    aes(
      y = wbar,
      group = interaction(n.pop0, bottleneck, ndd)
    ),
    colour = 'black',
  )  +
  geom_ribbon(
    data = . %>% filter(ext.gen %in% 15),
    aes(
      ymin = wbar - 2 * sqrt(wvar / n),
      ymax = wbar + 2 * sqrt(wvar / n),
      group = interaction(n.pop0, bottleneck, ndd)
    ),
    fill = 'black',
    alpha = 0.1
  ) +
  scale_colour_gradient(low = 'lightpink1', high = 'darkred') +
  scale_fill_gradient(low = 'lightpink1', high = 'darkred') +
  facet_wrap(ndd ~ interaction(n.pop0, bottleneck), ncol = 4) +
  theme_bw()

gen.d %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = wbar,
      group = interaction(n.pop0, ext.gen),
      colour = ext.gen # factor(ext.gen < 15)
    )
  ) +
  scale_color_gradient(high = 'black', low = 'red') +
  facet_wrap(~ n.pop0)
