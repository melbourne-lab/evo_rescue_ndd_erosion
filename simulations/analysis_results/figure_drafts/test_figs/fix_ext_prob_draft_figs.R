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
