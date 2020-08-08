library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

setwd('simulations/outputs/thrice_bottlenecked_results/')

### Population size data

a0n = read.csv('a000_evotrue_all_n.csv') %>%
  mutate(ndd = FALSE) %>%
  merge(y = read.csv('a000_evotrue_all_demo.csv') %>%
          select(n.pop0, bottleneck, gen, n.trials))
a3n = read.csv('a035_evotrue_all_n.csv') %>%
  mutate(ndd = TRUE) %>%
  merge(y = read.csv('a035_evotrue_all_demo.csv') %>%
          select(n.pop0, bottleneck, gen, n.trials))

alln = rbind(a0n, a3n)

## All treatment figure
ggplot(alln, aes(x = gen)) +
  geom_segment(
    aes(x = 1, xend = max(gen), y = n.pop0, yend = n.pop0),
    linetype = 2,
    colour = 'gray'
  ) +
  geom_line(
    aes(y = nbar,
        group = interaction(ndd, n.pop0, bottleneck),
        colour = factor(ndd),
        linetype = bottleneck)
  ) +
  geom_ribbon(
    aes(ymin = nbar - 2 * sqrt(nvar / n.trials),
        ymax = nbar + 2 * sqrt(nvar / n.trials),
        group = interaction(ndd, n.pop0, bottleneck),
        colour = factor(ndd),
        fill = factor(ndd)),
    alpha = 0.2,
    size = 0.125
  ) +
  scale_color_manual(values = c('black', 'purple'),
                     labels = c("Density\nindependent",
                                "Density\ndependent"),
                     name = "Growth form") +
  scale_fill_manual(values = c('black', 'purple'),
                    labels = c("Density\nindependent",
                               "Density\ndependent"),
                    name = "Growth form") +
  scale_linetype(name = "Genetic variance",
                 labels = c("High\n(no bottleneck)",
                            "Low\n(bottlenecked)")) +
  labs(x = 'Generation',
       y = 'Mean population size') +
  scale_y_log10(limits = c(9, 1180)) +
  theme(panel.background = element_blank(),
        legend.position = 'bottom')


##

ext.plot = alln %>%
  group_by(n.pop0, bottleneck, ndd) %>%
  arrange(n.pop0, bottleneck, ndd, gen) %>%
  mutate(n.extinct = c(0, -1 * diff(n.trials)),
         ins.p.ext1 = n.extinct / n.trials,
         ins.p.ext2 = n.extinct / max(n.trials),
         cum.p.ext  = n.trials / max(n.trials),
         con.p.ext  = n.extinct / max(n.extinct)) %>%
  mutate(varn = ifelse(bottleneck, 'High variation', 'Low variation'),
         pop0 = ifelse(n.pop0 < 50, 'Initially small', 'Initially large')) %>%
  ggplot(aes(x = gen))

ext.plot +
  geom_ribbon(
    aes(ymin = 0,
        ymax = 1 - cum.p.ext,
        group = interaction(ndd, n.pop0, bottleneck),
        fill = ndd),
    alpha = 0.15
  ) +
  geom_line(
    aes(y = ins.p.ext1,
        group = interaction(ndd, n.pop0, bottleneck),
        colour = ndd)
  ) +
  geom_ribbon(
    aes(ymin = ins.p.ext1 - 2 * sqrt(ins.p.ext1 * (1 - ins.p.ext1) / n.trials),
        ymax = ins.p.ext1 + 2 * sqrt(ins.p.ext1 * (1 - ins.p.ext1) / n.trials),
        group = interaction(ndd, n.pop0, bottleneck),
        fill = ndd),
    alpha = 0.3
  ) +
  scale_color_manual(values = c('black', 'purple'),
                     labels = c("Density\nindependent",
                                "Density\ndependent"),
                     name = "Growth form") +
  scale_fill_manual(values = c('black', 'purple'),
                    labels = c("Density\nindependent",
                               "Density\ndependent"),
                    name = "Growth form") +
  facet_grid(pop0 ~ varn) +
  theme(panel.background = element_blank(),
        legend.position = 'bottom')

### Demographic data

a0demo = read.csv('a000_evotrue_all_demo.csv') %>% 
  mutate(ndd = FALSE) %>%
  merge(y = read.csv('a000_evotrue_all_n.csv') %>%
          select(n.pop0, bottleneck, gen, pext))
a3demo = read.csv('a035_evotrue_all_demo.csv') %>% 
  mutate(ndd = TRUE) %>%
  merge(y = read.csv('a035_evotrue_all_n.csv') %>%
          select(n.pop0, bottleneck, gen, pext))

alldemo = rbind(a0demo, a3demo)

ggplot(alldemo, aes(x = gen)) +
  geom_segment(
    aes(x = 1, xend = 15, y = 1, yend = 1),
    linetype = 2, colour = 'gray'
  ) +
  geom_line(
    aes(y = zbar,
        group = interaction(n.pop0, ndd, bottleneck),
        colour = ndd)
  ) +
  geom_ribbon(
    aes(ymin = zbar - 2 * sqrt(zvar / n.trials),
        ymax = zbar + 2 * sqrt(zvar / n.trials),
        group = interaction(n.pop0, ndd, bottleneck),
        fill = ndd),
    alpha = 0.2
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  facet_wrap(n.pop0 ~ bottleneck) +
  theme(panel.background = element_blank(),
        legend.position = 'bottom')

#

ggplot(alldemo, aes(x = gen)) +
  geom_segment(
    aes(x = 1, xend = 15, y = 1, yend = 1),
    linetype = 2, colour = 'gray'
  ) +
  geom_line(
    aes(y = wbar,
        group = interaction(n.pop0, ndd, bottleneck),
        colour = ndd)
  ) +
  geom_ribbon(
    aes(ymin = wbar - 2 * sqrt(wvar / n.trials),
        ymax = wbar + 2 * sqrt(wvar / n.trials),
        group = interaction(n.pop0, ndd, bottleneck),
        fill = ndd),
    alpha = 0.2
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  facet_wrap(n.pop0 ~ bottleneck) +
  theme(panel.background = element_blank(),
        legend.position = 'bottom')

### Genetic data

a0gene = read.csv('a000_evotrue_all_gene.csv') %>%
  mutate(ndd = FALSE) %>%
  merge(y = read.csv('a000_evotrue_all_demo.csv') %>%
          select(n.pop0, bottleneck, gen, n.trials))
a3gene = read.csv('a035_evotrue_all_gene.csv') %>%
  mutate(ndd = TRUE) %>%
  merge(y = read.csv('a035_evotrue_all_demo.csv') %>%
          select(n.pop0, bottleneck, gen, n.trials))

allgene = rbind(a0gene, a3gene)

## All treatment plot

ggplot(allgene, aes(x = gen)) +
  geom_line(
    aes(y = vbar,
        group = interaction(ndd, n.pop0, bottleneck),
        colour = interaction(ndd, n.pop0))
  ) +
  geom_ribbon(
    aes(ymin = vbar - 2 * sqrt(vvar / n.trials),
        ymax = vbar + 2 * sqrt(vvar / n.trials),
        group = interaction(ndd, n.pop0, bottleneck),
        fill = interaction(ndd, n.pop0)),
    alpha = 0.2
  ) +
  scale_colour_manual(values = c('gray11', 'plum', 'black', 'purple'),
                      labels = c("Density independent, initially small",
                                 "Density dependent, initially small",
                                 "Density independent, initially large",
                                 "Density dependent, initially small"),
                      name = "Treatment") +
  scale_fill_manual(values = c('gray11', 'plum', 'black', 'purple'),
                    labels = c("Density independent, initially small",
                               "Density dependent, initially small",
                               "Density independent, initially large",
                               "Density dependent, initially small"),
                    name = "Treatment") +
  guides(color = guide_legend(ncol = 2, byrow = TRUE)) +
  labs(x = "Generation", y = "Genetic variance") +
  theme(panel.background = element_blank(),
        legend.position = 'bottom') +
  ggsave('~/Documents/Research/boulder/scott_rescue/esa_2020/genvar_fig.pdf',
         width = 10, height = 5)

## Only high-variation plot
ggplot(allgene %>% filter(!bottleneck), aes(x = gen)) +
  geom_line(
    aes(y = vbar,
        group = interaction(ndd, n.pop0, bottleneck),
        colour = interaction(ndd, n.pop0))
  ) +
  geom_ribbon(
    aes(ymin = vbar - 2 * sqrt(vvar / n.trials),
        ymax = vbar + 2 * sqrt(vvar / n.trials),
        group = interaction(ndd, n.pop0, bottleneck),
        fill = interaction(ndd, n.pop0)),
    alpha = 0.2
  ) +
  scale_colour_manual(values = c('gray11', 'plum1', 'black', 'purple'),
                      labels = c("Density independent, initially small",
                                 "Density dependent, initially small",
                                 "Density independent, initially large",
                                 "Density dependent, initially small"),
                      name = "Treatment") +
  scale_fill_manual(values = c('gray11', 'plum1', 'black', 'purple'),
                    labels = c("Density independent, initially small",
                               "Density dependent, initially small",
                               "Density independent, initially large",
                               "Density dependent, initially small"),
                    name = "Treatment") +
  guides(color = guide_legend(ncol = 2, byrow = TRUE)) +
  labs(x = "Generation", y = "Genetic variance") +
  lims(y = c(0.25, 0.50)) +
  theme(panel.background = element_blank(),
        legend.position = 'bottom') +
  ggsave('~/Documents/Research/boulder/scott_rescue/esa_2020/genvar_fig_highvar_only.pdf',
         width = 10, height = 5)

# Density independent only
ggplot(allgene %>% filter(!bottleneck & !ndd), aes(x = gen)) +
  geom_line(
    aes(y = vbar,
        group = interaction(ndd, n.pop0, bottleneck),
        colour = interaction(ndd, n.pop0))
  ) +
  geom_ribbon(
    aes(ymin = vbar - 2 * sqrt(vvar / n.trials),
        ymax = vbar + 2 * sqrt(vvar / n.trials),
        group = interaction(ndd, n.pop0, bottleneck),
        fill = interaction(ndd, n.pop0)),
    alpha = 0.2
  ) +
  scale_colour_manual(values = c('gray11', 'black'),
                      labels = c("Density independent, initially small",
                                 "Density independent, initially large"),
                      name = "Treatment") +
  scale_fill_manual(values = c('gray11', 'black'),
                    labels = c("Density independent, initially small",
                               "Density independent, initially large"),
                    name = "Treatment") +
  labs(x = "Generation", y = "Genetic variance") +
  guides(color = guide_legend(ncol = 1)) +
  lims(y = c(0.25, 0.50)) +
  theme(panel.background = element_blank(),
        legend.position = 'bottom') +
  ggsave('~/Documents/Research/boulder/scott_rescue/esa_2020/genvar_fig_highvar_black_only.pdf',
         width = 10, height = 5)
