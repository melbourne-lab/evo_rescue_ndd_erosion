library(ggplot2)
library(dplyr)
library(tidyr)

all.d = rbind(
  read.csv('simulations/outputs/thrice_bottlenecked_results/a000_hivar_extinctions.csv') %>%
    mutate(bottleneck = FALSE, ndd = FALSE),
  read.csv('simulations/outputs/thrice_bottlenecked_results/a000_lovar_extinctions.csv') %>%
    mutate(bottleneck = TRUE, ndd = FALSE),
  read.csv('simulations/outputs/thrice_bottlenecked_results/a035_hivar_extinctions.csv') %>%
    mutate(bottleneck = FALSE, ndd = TRUE),
  read.csv('simulations/outputs/thrice_bottlenecked_results/a035_lovar_extinctions.csv') %>%
    mutate(bottleneck = TRUE, ndd = TRUE)
) %>%
  group_by(bottleneck, ndd, n.pop0) %>%
  mutate(n.trials = pext * 24000,
         d.trials = -c(0, diff(n.trials)),
         p.ins1 = c(d.trials[-1], 0) / n.trials,
         p.ins2 = c(d.trials[-1], 0) / 24000,
         p.ins3 = c(d.trials[-1], 0) / sum(d.trials)) %>% # (24000 - min(n.trials))) %>%
  filter(gen < 15) %>%
  ungroup() %>%
  mutate(n.pop0 = factor(n.pop0, labels = c('Initially small', 'Initially large')),
         bottleneck = factor(bottleneck, labels = c('High variation', 'Low variation')),
         ndd = factor(ndd, labels = c('Density independent', 'Density dependent')))

head(all.d)

##### Plots

### Conditional (instantaneous) probability
# Prob. of going extinct in time step t given
# being extant in time t-1

# With ribbon for overall extinction risk

all.d %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes( 
      y = p.ins1, 
      group = interaction(n.pop0, bottleneck, ndd),
      colour = ndd
    )
  ) +
  geom_ribbon(
    aes(
      ymin = p.ins1 - 2 * sqrt(p.ins1 * (1-p.ins1) / n.trials),
      ymax = p.ins1 + 2 * sqrt(p.ins1 * (1-p.ins1) / n.trials),
      group = interaction(n.pop0, bottleneck, ndd),
      fill = ndd
    ),
    alpha = 0.2) +
  geom_ribbon(
    aes(
      xmin = 0, xmax = 14,
      ymin = 0, ymax = 1 - pext,
      group = interaction(n.pop0, bottleneck, ndd),
      fill = ndd
    ),
    alpha = 0.1
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  facet_wrap(n.pop0 ~ bottleneck)

# Without the ribbon
all.d %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes( 
      y = p.ins1, 
      group = interaction(n.pop0, bottleneck, ndd),
      colour = ndd
    )
  ) +
  geom_ribbon(
    aes(
      ymin = p.ins1 - 2 * sqrt(p.ins1 * (1-p.ins1) / n.trials),
      ymax = p.ins1 + 2 * sqrt(p.ins1 * (1-p.ins1) / n.trials),
      group = interaction(n.pop0, bottleneck, ndd),
      fill = ndd
    ),
    alpha = 0.2
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  facet_wrap(n.pop0 ~ bottleneck)

##

all.d %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes( 
      y = p.ins2, 
      group = interaction(n.pop0, bottleneck, ndd),
      colour = ndd
    )
  ) +
  geom_ribbon(
    aes(
      ymin = p.ins2 - 2 * sqrt(p.ins2 * (1-p.ins2) / n.trials),
      ymax = p.ins2 + 2 * sqrt(p.ins2 * (1-p.ins2) / n.trials),
      group = interaction(n.pop0, bottleneck, ndd),
      fill = ndd
    ),
    alpha = 0.2) +
  geom_ribbon(
    aes(
      xmin = 0, xmax = 14,
      ymin = 0, ymax = 1 - pext,
      group = interaction(n.pop0, bottleneck, ndd),
      fill = ndd
    ),
    alpha = 0.1
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  facet_wrap(n.pop0 ~ bottleneck)

# Extinction risk conditioned on being extinct
# (when does extinction happen if it happens at all)

all.d %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes( 
      y = p.ins3, 
      group = interaction(n.pop0, bottleneck, ndd),
      colour = ndd
    )
  ) +
  # geom_ribbon(
  #   aes(
  #     ymin = p.ins3 - 2 * sqrt(p.ins3 * (1-p.ins3) / n.trials),
  #     ymax = p.ins3 + 2 * sqrt(p.ins3 * (1-p.ins3) / n.trials),
  #     group = interaction(n.pop0, bottleneck, ndd),
  #     fill = ndd
  #   ),
  #   alpha = 0.2) +
  # geom_ribbon(
  #   aes(
  #     xmin = 0, xmax = 14,
  #     ymin = 0, ymax = 1 - pext,
  #     group = interaction(n.pop0, bottleneck, ndd),
  #     fill = ndd
  #   ),
  #   alpha = 0.1
  # ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  facet_wrap(n.pop0 ~ bottleneck)


###

# Population size

all.d %>%
  mutate(init.size = paste0('Initially ',
                            ifelse(n.pop0 > 50, 'large', 'small'))) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = nbar,
      group = interaction(n.pop0, bottleneck, ndd),
      colour = ndd
    )
  ) +
  geom_ribbon(
    aes(
      ymin = nbar - 2 * sqrt(nvar / n.trials),
      ymax = nbar + 2 * sqrt(nvar / n.trials),
      group = interaction(n.pop0, bottleneck, ndd),
      fill = ndd
    ),
    alpha = 0.2
  ) +
  scale_color_manual(values = c('black', 'slateblue1'),
                     labels = c("Density\nindependent",
                                "Density\ndependent"),
                     name = "Growth form") +
  scale_fill_manual(values = c('black', 'slateblue1'),
                    labels = c("Density\nindependent",
                               "Density\ndependent"),
                    name = "Growth form") +
  scale_linetype(name = "Genetic variance",
                 labels = c("High\n(no bottleneck)",
                            "Low\n(bottlenecked)")) +
  labs(x = 'Generation',
       y = 'Mean population size') +
  scale_y_log10() +
  facet_wrap(~ init.size) +
  theme(panel.background = element_blank(),
        legend.position = 'bottom') 

### Extinction figure (with cowplot)

library(cowplot)

ins.probs = ggplot(all.d, aes(x = gen)) +
  geom_line(
    aes( 
      y = p.ins1, 
      group = interaction(n.pop0, bottleneck, ndd),
      colour = ndd
    )
  ) +
  geom_ribbon(
    aes(
      ymin = p.ins1 - 2 * sqrt(p.ins1 * (1-p.ins1) / n.trials),
      ymax = p.ins1 + 2 * sqrt(p.ins1 * (1-p.ins1) / n.trials),
      group = interaction(n.pop0, bottleneck, ndd),
      fill = ndd
    ),
    alpha = 0.2
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  facet_wrap( ~ paste(n.pop0, bottleneck, sep = ', '), ncol = 1) +
  labs(x = 'Generation', y = 'Probability of extinction') +
  theme(legend.position = 'none',
        panel.background = element_blank())

cum.probs = ggplot(all.d, aes(x = gen)) +
  geom_ribbon(
    aes(
      xmin = 0, xmax = 14,
      ymin = 0, ymax = 1 - pext,
      group = interaction(n.pop0, bottleneck, ndd),
      fill = ndd
    ),
    alpha = 0.1
  ) +
  scale_fill_manual(values = c('black', 'purple')) +
  labs(x = 'Generation') +
  facet_wrap( ~ paste(n.pop0, bottleneck, sep = ', '), ncol = 1) +
  theme(legend.position = 'none',
        panel.background = element_blank())

extinct.legend = get_legend(
  ins.probs + 
    guides(fill = guide_legend(''), colour = guide_legend('')) +
    theme(legend.position = 'bottom',
          legend.box.margin = margin(1, 0, 0, 0))
  )

data.plots = plot_grid(ins.probs, cum.probs)
plot_grid(data.plots, extinct.legend, ncol = 1, rel_heights = c(1, .1))

# plot_grid(ins.probs, cum.probs, extinct.legend)
