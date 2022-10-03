library(cowplot)
library(dplyr)
library(tidyr)

bb = rbind(read.csv('simulations/outputs/test_gridded_v.csv') %>% mutate(ndd = TRUE), 
           read.csv('simulations/outputs/test_gridded_v2.csv') %>% mutate(ndd = FALSE))

head(bb)

bb.pos = bb %>%
  group_by(trial) %>%
  mutate(ext.next = extinct & (gen == ext.gen-1)) %>%
  filter(gen < 15) %>%
  group_by(p.fix.pos, n.pop0, ndd) %>%
  summarise(pexnex = mean(ext.next),
            n = n())

pos.plot = bb.pos %>%
  ggplot(aes(x = p.fix.pos)) +
  geom_line(
    aes(
      y = pexnex,
      group = interaction(ndd, n.pop0),
      linetype = factor(n.pop0),
      colour = ndd
      )
    ) +
  geom_ribbon(
    aes(
      ymin = pexnex - 2 * sqrt(pexnex * (1-pexnex) / n),
      ymax = pexnex + 2 * sqrt(pexnex * (1-pexnex) / n),
      group = interaction(ndd, n.pop0),
      fill = ndd
    ),
    alpha = 0.1
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  scale_linetype_manual(values = 2:3, 
                        labels = c("Small", "Large"),
                        'Size') +
  labs(x = 'Loci at negative fixation',
       y = 'Prob. of extinction') +
  ylim(-0.05, 0.2)

bb.neg = bb %>%
  group_by(trial) %>%
  mutate(ext.next = extinct & (gen == ext.gen-1)) %>%
  filter(gen < 15) %>%
  group_by(p.fix.neg, n.pop0, ndd) %>%
  summarise(pexnex = mean(ext.next),
            n = n())

neg.plot = bb.neg %>%
  ggplot(aes(x = p.fix.neg)) +
  geom_line(
    aes(
      y = pexnex,
      group = interaction(ndd, n.pop0),
      linetype = factor(n.pop0),
      colour = ndd
    )
  ) +
  geom_ribbon(
    aes(
      ymin = pexnex - 2 * sqrt(pexnex * (1-pexnex) / n),
      ymax = pexnex + 2 * sqrt(pexnex * (1-pexnex) / n),
      group = interaction(ndd, n.pop0),
      fill = ndd
    ),
    alpha = 0.1
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  scale_linetype_manual(values = 2:3, 
                        labels = c("Small", "Large"),
                        'Size') +
  labs(x = 'Loci at negative fixation',
       y = 'Prob. of extinction') +
  ylim(-0.05, 0.2)

ext.leg.se = get_legend(pos.plot)

plot_grid(pos.plot + theme(legend.position = 'none'), 
          neg.plot + theme(legend.position = 'none'), 
          ext.leg.se,
          nrow = 1,
          rel_widths = c(1, 1, 0.2))

save_plot(
  plot_grid(pos.plot + theme(legend.position = 'none'), 
            neg.plot + theme(legend.position = 'none'), 
            ext.leg.se,
            nrow = 1,
            rel_widths = c(1, 1, 0.2)),
  file = 'simulations/analysis_results/figure_drafts/test_figs/fix_ext_prob_err.pdf',
  base_width = 8
)

### Plot without errors.
          
pos.plot.ner = bb.pos %>%
  ggplot(aes(x = p.fix.pos)) +
  geom_line(
    aes(
      y = pexnex,
      group = interaction(ndd, n.pop0),
      linetype = factor(n.pop0),
      colour = ndd
    )
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_linetype_manual(values = 2:3, 
                        labels = c("Small", "Large"),
                        'Size') +
  labs(x = 'Loci at negative fixation',
       y = 'Prob. of extinction') +
  ylim(0, 0.2)

neg.plot.ner = bb.neg %>%
  ggplot(aes(x = p.fix.neg)) +
  geom_line(
    aes(
      y = pexnex,
      group = interaction(ndd, n.pop0),
      linetype = factor(n.pop0),
      colour = ndd
    )
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  scale_linetype_manual(values = 2:3, 
                        labels = c("Small", "Large"),
                        'Size') +
  labs(x = 'Loci at negative fixation',
       y = 'Prob. of extinction') +
  ylim(0, 0.2)

ext.leg.ner = get_legend(pos.plot.ner)

plot_grid(pos.plot.ner + theme(legend.position = 'none'), 
          neg.plot.ner + theme(legend.position = 'none'), 
          ext.leg.ner,
          nrow = 1,
          rel_widths = c(1, 1, 0.2))

save_plot(
  plot_grid(pos.plot.ner + theme(legend.position = 'none'), 
            neg.plot.ner + theme(legend.position = 'none'), 
            ext.leg.ner,
            nrow = 1,
            rel_widths = c(1, 1, 0.2)),
  file = 'simulations/analysis_results/figure_drafts/test_figs/fix_ext_prob_nerr.pdf',
  base_width = 8
)

##### Updated plots with recent (Sept ~20) batch of simulations.

library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

### Clear namespace

rm(list = ls())

### Global parameters

trials = 4000

pars = data.frame(
  end.time = 15,
  n.loci = 25,
  w.max = 2,
  theta = 2.75,
  wfitn = sqrt(1 / 0.14 / 2),
  sig.e = sqrt(0.5),
  pos.p = 0.5,
  alpha = 0
)

### Read in all data

all.data = read.csv('simulations/outputs/final_results/alldata_combined.csv')

### Do the gridding thing.

grid.y.res = 0.025

any.ext = all.data %>%
  # Now, bin variances and aggregate extinction probabilities
  mutate(vr = round((1/grid.y.res) * v) / (1/grid.y.res)) %>%
  group_by(gen, n.pop0, low.var, alpha, vr) %>%
  summarise(p.extinct = mean(extinct),
            n = n())

### Plot

any.ext %>%
  filter(gen < 15) %>%
  ungroup() %>%
  # mutate(ndd = factor(ndd, labels = c("Density independent", "Density dependent")),
  #        n.pop0 = factor(n.pop0, labels = c("Initially small", "Initially large"))) %>%
  ggplot() +
  geom_tile(
    aes(
      x = gen,
      y = vr,
      fill = p.extinct
    )
  ) +
  facet_grid(rows = vars(low.var, n.pop0), cols = vars(alpha)) +
  scale_fill_viridis_b(option = 'B') +
  labs(x = 'Generation', y = 'Genetic variance') +
  guides(fill = guide_colorbar('Extinction\nprobability', barwidth = unit(3, 'cm'))) +
  theme(panel.background = element_rect(fill = NA),
        legend.position = 'bottom') +
  ggsave(file = 'simulations/analysis_results/figure_drafts/test_figs/eg_extinction_prob_over_time.pdf')

### N.b. Ruth was curious once about *instantaneous extinction*

ins.ext = all.data %>%
  # Instantaneous extinction risk
  group_by(trial,  n.pop0, low.var, alpha) %>%
  mutate(extinct.now = (gen == max(gen)) & extinct) %>%
  # Now, bin variances and aggregate instantaneous risk
  mutate(vr = round((1/grid.y.res) * v) / (1/grid.y.res)) %>%
  group_by(gen, n.pop0, low.var, alpha, vr) %>%
  summarise(p.extinct = mean(extinct.now),
            n = n())

ins.ext %>%
  filter(gen < 15) %>%
  ungroup() %>%
  # mutate(ndd = factor(ndd, labels = c("Density independent", "Density dependent")),
  #        n.pop0 = factor(n.pop0, labels = c("Initially small", "Initially large"))) %>%
  ggplot() +
  geom_tile(
    aes(
      x = gen,
      y = vr,
      fill = p.extinct
    )
  ) +
  facet_grid(rows = vars(low.var, n.pop0), cols = vars(alpha)) +
  scale_fill_viridis_b(option = 'B') +
  labs(x = 'Generation', y = 'Genetic variance') +
  guides(fill = guide_colorbar('Extinction\nprobability', barwidth = unit(3, 'cm'))) +
  theme(panel.background = element_rect(fill = NA),
        legend.position = 'bottom') +
  ggsave(file = 'simulations/analysis_results/figure_drafts/test_figs/eg_ins_extinction_prob_over_time.pdf')

# Try separating out by NDD to get differences in extinction probability.

any.ext.alpha = any.ext %>%
  pivot_wider(names_from = alpha, values_from = c(p.extinct, n)) %>%
  mutate(p.diff = p.extinct_0 - p.extinct_0.0035)

any.ext.alpha %>%
  ggplot() +
  geom_tile(
    aes(
      x = gen,
      y = vr,
      fill = p.diff
    )
  ) +
  facet_grid(rows = vars(low.var, n.pop0)) +
  scale_fill_gradient2(high = 'blue', mid = 'white', midpoint = 0, low = 'red')
  
# Not great...

# Try a circle plot...

any.ext %>%
  filter(gen < 15) %>%
  ggplot() +
  geom_point(
    aes(
      x = gen,
      y = vr,
      colour = p.extinct,
      size = log(n)
    )
  ) +
  scale_colour_viridis_b(option = 'B') +
  facet_grid(rows = vars(low.var, n.pop0), cols = vars(alpha)) +
  labs(x = 'Generation', y = 'Genetic variance') +
  guides(fill = guide_colorbar('Extinction\nprobability', barwidth = unit(3, 'cm')),
         size = 'none') +
  theme(panel.background = element_rect(fill = NA),
        legend.position = 'bottom') +
  ggsave(file = 'simulations/analysis_results/figure_drafts/test_figs/eg_extinction_pt_size.pdf',
         height = 8, width = 5)

# ^ I kinda like this one.

# Instantaneous extinctions 

ins.ext %>%
  filter(gen < 15) %>%
  ggplot() +
  geom_point(
    aes(
      x = gen,
      y = vr,
      colour = p.extinct,
      size = log(n)
    )
  ) +
  scale_colour_viridis_b(option = 'B') +
  facet_grid(rows = vars(low.var, n.pop0), cols = vars(alpha)) +
  labs(x = 'Generation', y = 'Genetic variance') +
  guides(fill = guide_colorbar('Extinction\nprobability', barwidth = unit(3, 'cm')),
         size = 'none') +
  theme(panel.background = element_rect(fill = NA),
        legend.position = 'bottom') +
  ggsave(file = 'simulations/analysis_results/figure_drafts/test_figs/eg_ins_extinction_pt_size.pdf',
         height = 8, width = 5)
o