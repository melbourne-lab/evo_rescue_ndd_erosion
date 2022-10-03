library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

setwd('')

### Load in parameters

global.pars = data.frame(
  end.time = 15,
  init.row = 1e4,
  n.loci = 25,
  w.max = 2,
  theta = 2.6,
  wfitn = sqrt(1 / 0.14 / 2),
  sig.e = 0.5,
  pos.p = 0.5,
  alpha = 0.0035
)

### Population size data

a0n = read.csv('simulations/outputs/thrice_bottlenecked_results/a000_evotrue_all_n.csv') %>%
  mutate(ndd = FALSE) %>%
  merge(y = read.csv('simulations/outputs/thrice_bottlenecked_results/a000_evotrue_all_demo.csv') %>%
          select(n.pop0, bottleneck, gen, n.trials))
a3n = read.csv('simulations/outputs/thrice_bottlenecked_results/a035_evotrue_all_n.csv') %>%
  mutate(ndd = TRUE) %>%
  merge(y = read.csv('simulations/outputs/thrice_bottlenecked_results/a035_evotrue_all_demo.csv') %>%
          select(n.pop0, bottleneck, gen, n.trials))

alln = rbind(a0n, a3n)

## All treatment figure
## This would be Figure 1 in manuscript (or at least first fig in results)

ggplot(alln %>% 
         mutate(init.size = factor(ifelse(n.pop0 < 50, 'Initially small', 'Initially large'))), 
       aes(x = gen)) +
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
  scale_y_log10(limits = c(9, 1180)) +
  facet_wrap(~ init.size) +
  theme(panel.background = element_blank(),
        legend.position = 'bottom') +
  ggsave('simulations/analysis_results/figure_drafts/draft_figs/fig_1.pdf')


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

a0demo = read.csv('simulations/outputs/thrice_bottlenecked_results/a000_evotrue_all_demo.csv') %>% 
  mutate(ndd = FALSE) %>%
  merge(y = read.csv('simulations/outputs/thrice_bottlenecked_results/a000_evotrue_all_n.csv') %>%
          select(n.pop0, bottleneck, gen, pext))
a3demo = read.csv('simulations/outputs/thrice_bottlenecked_results/a035_evotrue_all_demo.csv') %>% 
  mutate(ndd = TRUE) %>%
  merge(y = read.csv('simulations/outputs/thrice_bottlenecked_results/a035_evotrue_all_n.csv') %>%
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

a0gene = read.csv('simulations/outputs/thrice_bottlenecked_results/a000_evotrue_all_gene.csv') %>%
  mutate(ndd = FALSE) %>%
  merge(y = read.csv('simulations/outputs/thrice_bottlenecked_results/a000_evotrue_all_demo.csv') %>%
          select(n.pop0, bottleneck, gen, n.trials))
a3gene = read.csv('simulations/outputs/thrice_bottlenecked_results/a035_evotrue_all_gene.csv') %>%
  mutate(ndd = TRUE) %>%
  merge(y = read.csv('simulations/outputs/thrice_bottlenecked_results/a035_evotrue_all_demo.csv') %>%
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

# Re-plot genotypes but with G&H estimates

pred.ds = expand.grid(bottleneck = c(TRUE, FALSE),
                      n.pop0 = c(20, 60),
                      gen = 1:global.pars$end.time) %>%
  merge(y = allgene %>% 
          filter(gen %in% 1, ndd) %>%
          select(bottleneck, n.pop0, vbar),
        by = c('bottleneck', 'n.pop0')) %>%
  rename(v1 = vbar) %>%
  mutate(
       k = with(global.pars, (wfitn^2 + sig.e^2) / (wfitn^2 + sig.e^2 + v1)),
       dt = global.pars$theta * k^(gen-1),
       zt = global.pars$theta - dt
       )

pred.ds

ggplot(alldemo, aes(x = gen)) +
  geom_segment(
    aes(x = 1, xend = 15, y = 1, yend = 1),
    linetype = 3, colour = 'gray'
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
  geom_line(
    data = pred.ds,
    aes(x = gen, y = zt),
    linetype = 2
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  facet_wrap(n.pop0 ~ bottleneck) +
  theme(panel.background = element_blank(),
        legend.position = 'bottom')

# Well... G&H doesn't match at all.
# What the heck is happening here? 

ggplot(alldemo, aes(x = gen)) +
  geom_segment(
    aes(x = 1, xend = 15, y = 1, yend = 1),
    linetype = 3, colour = 'gray'
  ) +
  geom_line(
    aes(y = zbar,
        group = interaction(n.pop0, ndd, bottleneck),
        colour = ndd,
        linetype = bottleneck)
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
  facet_wrap(n.pop0 ~ .) +
  theme(panel.background = element_blank(),
        legend.position = 'bottom')
