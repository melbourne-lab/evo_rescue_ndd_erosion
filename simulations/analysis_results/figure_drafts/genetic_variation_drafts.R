### Plotting genetic variaton in rescue simulations.
### Scripts: simulations/run_final_sims/*extinctions.R
### SN - run September 1 2020, plots September 2 2020

### Load packages

library(ggplot2)
library(dplyr)
library(tidyr)

### Case 1: means conditioned on being extant in given time

# Load in data
all.g = rbind(
  read.csv('simulations/outputs/final_results/a000_hivar_gene_all.csv')  %>%
    mutate(hi.var = TRUE, ndd = FALSE),
  read.csv('simulations/outputs/final_results/a035_hivar_gene_all.csv') %>%
    mutate(hi.var = TRUE, ndd = TRUE),
  read.csv('simulations/outputs/final_results/a000_lowvar_gene_all.csv') %>%
    mutate(hi.var = FALSE, ndd = FALSE),
  read.csv('simulations/outputs/final_results/a035_lowvar_gene_all.csv') %>%
    mutate(hi.var = FALSE, ndd = TRUE)
) %>%
  mutate(n.pop0 = factor(n.pop0, labels = c('Initially small', 'Initially large')),
         hi.var = factor(hi.var, labels = c('Low variation', 'High variation')),
         ndd = factor(ndd, labels = c('Density independent', 'Density dependent')))

# looka that
head(all.g)

# Try a plot of all genotypes at once
all.g %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = vbar,
      group = interaction(n.pop0, hi.var, ndd),
      linetype = paste(n.pop0, hi.var, sep = ', '),
      colour = ndd
    )
  ) +
  geom_ribbon(
    aes(
      ymin = vbar - 2 * sqrt(vvar / n),
      ymax = vbar + 2 * sqrt(vvar / n),
      group = interaction(n.pop0, hi.var, ndd),
      fill = ndd
    ),
    alpha = 0.1
  ) +
  scale_y_continuous(limits = c(0, 0.5)) +
  scale_fill_manual(values = c('black', 'purple')) +
  scale_colour_manual(values = c('black', 'purple')) +
  scale_linetype_manual(values = c(6, 4, 5, 2)) +
  labs(x = 'Generation', y = 'Genetic variation') +
  guides(linetype = guide_legend('', nrow = 2),
         colour = guide_legend("Growth form", nrow = 2),
         fill = guide_legend("Growth form", nrow = 2)) +
  theme(panel.background = element_rect(fill = 'white'),
        panel.grid = element_line(colour = 'gray88'),
        panel.border = element_rect(fill = NA),
        legend.position = 'bottom')

### Conditioned on being extinct/extant at end of simulation

ext.g = rbind(
  read.csv('simulations/outputs/final_results/a000_hivar_gene_ext.csv')  %>%
    mutate(hi.var = TRUE, ndd = FALSE),
  read.csv('simulations/outputs/final_results/a035_hivar_gene_ext.csv') %>%
    mutate(hi.var = TRUE, ndd = TRUE),
  read.csv('simulations/outputs/final_results/a000_lowvar_gene_ext.csv') %>%
    mutate(hi.var = FALSE, ndd = FALSE),
  read.csv('simulations/outputs/final_results/a035_lowvar_gene_ext.csv') %>%
    mutate(hi.var = FALSE, ndd = TRUE)
) %>%
  mutate(n.pop0 = factor(n.pop0, labels = c('Initially small', 'Initially large')),
         hi.var = factor(hi.var, labels = c('Low variation', 'High variation')),
         ndd = factor(ndd, labels = c('Density independent', 'Density dependent')))


ext.g %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = vbar,
      group = interaction(n.pop0, hi.var, ndd, extinct),
      linetype = extinct,
      colour = ndd
    )
  ) +
  geom_ribbon(
    aes(
      ymin = vbar - 2 * sqrt(vvar / n),
      ymax = vbar + 2 * sqrt(vvar / n),
      group = interaction(n.pop0, hi.var, ndd, extinct),
      fill = ndd
    ),
    alpha = 0.1
  ) +
  facet_wrap(n.pop0 ~ hi.var) +
  guides(linetype = guide_legend("", nrow = 2),
         colour = guide_legend("Growth form", nrow = 2),
         fill = guide_legend("Growth form", nrow = 2)) +
  scale_linetype(labels = c("Extant", "Extinct")) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  labs(x = 'Generation', y = 'Genetic variation') +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        panel.grid = element_line(colour = 'gray88'),
        legend.position = 'bottom')

###

gen.g = rbind(
  read.csv('simulations/outputs/final_results/a000_hivar_gene_gen.csv')  %>%
    mutate(hi.var = TRUE, ndd = FALSE),
  read.csv('simulations/outputs/final_results/a035_hivar_gene_gen.csv') %>%
    mutate(hi.var = TRUE, ndd = TRUE),
  read.csv('simulations/outputs/final_results/a000_lowvar_gene_gen.csv') %>%
    mutate(hi.var = FALSE, ndd = FALSE),
  read.csv('simulations/outputs/final_results/a035_lowvar_gene_gen.csv') %>%
    mutate(hi.var = FALSE, ndd = TRUE)
) %>%
  mutate(n.pop0 = factor(n.pop0, labels = c('Initially small', 'Initially large')),
         hi.var = factor(hi.var, labels = c('Low variation', 'High variation')),
         ndd = factor(ndd, labels = c('Density independent', 'Density dependent')))

gen.g %>%
  ggplot(aes(x = gen)) +
  geom_line(
    data = . %>% filter(ext.gen < 15),
    aes(
      y = vbar,
      group = interaction(n.pop0, hi.var, ndd, ext.gen),
      colour = ext.gen
    )
  ) +
  geom_ribbon(
    data = . %>% filter(ext.gen < 15),
    aes(
      ymin = vbar - 2 * sqrt(vvar / n),
      ymax = vbar + 2 * sqrt(vvar / n),
      group = interaction(n.pop0, hi.var, ndd, ext.gen),
      fill = ext.gen
    ),
    alpha = 0.1
  ) +
  geom_line(
    data = . %>% filter(ext.gen %in% 15),
    aes(
      y = vbar,
      group = interaction(n.pop0, hi.var, ndd)
    ),
    colour = 'black'
  )  +
  geom_ribbon(
    data = . %>% filter(ext.gen %in% 15),
    aes(
      ymin = vbar - 2 * sqrt(vvar / n),
      ymax = vbar + 2 * sqrt(vvar / n),
      group = interaction(n.pop0, hi.var, ndd)
    ),
    fill = 'black',
    alpha = 0.1
  ) +
  scale_colour_gradient(low = 'lightpink1', high = 'darkred') +
  scale_fill_gradient(low = 'lightpink1', high = 'darkred') +
  guides(colour = guide_legend("Extinction\ngeneration"),
         fill = guide_legend('Extinction\ngeneration')) +
  labs(x = 'Generation', y = 'Genetic variation') +
  facet_wrap(ndd ~ paste(hi.var, n.pop0, sep = ', '), ncol = 4) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        panel.grid = element_line(colour = 'gray88'),
        legend.position = 'bottom')

### Other drafts

gen.g %>%
  ggplot(aes(x = gen)) +
  geom_line(
    data = . %>% filter(ext.gen < 15),
    aes(
      y = vbar,
      group = interaction(n.pop0, hi.var, ndd, ext.gen),
      linetype = ndd,
      colour = ext.gen
    )
  ) +
  geom_ribbon(
    data = . %>% filter(ext.gen < 15),
    aes(
      ymin = vbar - 2 * sqrt(vvar / n),
      ymax = vbar + 2 * sqrt(vvar / n),
      group = interaction(n.pop0, hi.var, ndd, ext.gen),
      fill = ext.gen
    ),
    alpha = 0.1
  ) +
  geom_line(
    data = . %>% filter(ext.gen %in% 15),
    aes(
      y = vbar,
      group = interaction(n.pop0, hi.var, ndd),
      linetype = ndd,
    ),
    colour = 'black'
  )  +
  geom_ribbon(
    data = . %>% filter(ext.gen %in% 15),
    aes(
      ymin = vbar - 2 * sqrt(vvar / n),
      ymax = vbar + 2 * sqrt(vvar / n),
      group = interaction(n.pop0, hi.var, ndd)
    ),
    fill = 'black',
    alpha = 0.1
  ) +
  scale_colour_gradient(low = 'lightpink1', high = 'darkred') +
  scale_fill_gradient(low = 'lightpink1', high = 'darkred') +
  facet_wrap(~ interaction(n.pop0, hi.var), ncol = 2) +
  theme_bw()

gen.g %>%
  mutate(p = n / 4000) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    data = . %>% filter(ext.gen < 15),
    aes(
      y = vbar,
      group = interaction(n.pop0, hi.var, ndd, ext.gen),
      colour = ext.gen,
      size = p
    )
  ) +
  geom_ribbon(
    data = . %>% filter(ext.gen < 15),
    aes(
      ymin = vbar - 2 * sqrt(vvar / n),
      ymax = vbar + 2 * sqrt(vvar / n),
      group = interaction(n.pop0, hi.var, ndd, ext.gen),
      fill = ext.gen
    ),
    alpha = 0.1
  ) +
  geom_line(
    data = . %>% filter(ext.gen %in% 15),
    aes(
      y = vbar,
      group = interaction(n.pop0, hi.var, ndd),
      size = p
    ),
    colour = 'black',
  )  +
  geom_ribbon(
    data = . %>% filter(ext.gen %in% 15),
    aes(
      ymin = vbar - 2 * sqrt(vvar / n),
      ymax = vbar + 2 * sqrt(vvar / n),
      group = interaction(n.pop0, hi.var, ndd)
    ),
    fill = 'black',
    alpha = 0.1
  ) +
  scale_colour_gradient(low = 'lightpink1', high = 'darkred') +
  scale_fill_gradient(low = 'lightpink1', high = 'darkred') +
  scale_size(range = c(0.1, 2)) +
  facet_wrap(ndd ~ interaction(n.pop0, hi.var), ncol = 4) +
  theme_bw()

gen.g %>%
  mutate(p = n / 4000) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    data = . %>% filter(ext.gen < 15),
    aes(
      y = vbar,
      group = interaction(n.pop0, hi.var, ndd, ext.gen),
      colour = p
    )
  ) +
  geom_ribbon(
    data = . %>% filter(ext.gen < 15),
    aes(
      ymin = vbar - 2 * sqrt(vvar / n),
      ymax = vbar + 2 * sqrt(vvar / n),
      group = interaction(n.pop0, hi.var, ndd, ext.gen)
      #fill = p
    ),
    alpha = 0.1
  ) +
  geom_line(
    data = . %>% filter(ext.gen %in% 15),
    aes(
      y = vbar,
      group = interaction(n.pop0, hi.var, ndd),
      colour = p
    ),
    size = 2
  )  +
  geom_ribbon(
    data = . %>% filter(ext.gen %in% 15),
    aes(
      ymin = vbar - 2 * sqrt(vvar / n),
      ymax = vbar + 2 * sqrt(vvar / n),
      group = interaction(n.pop0, hi.var, ndd)
      #fill = p
    ),
    alpha = 0.1
  ) +
  scale_colour_viridis_c(option = 'B') +
  scale_fill_viridis_c(option = 'B') +
  scale_size(range = c(0.1, 2)) +
  facet_wrap(ndd ~ interaction(n.pop0, hi.var), ncol = 4) +
  theme_bw()

gen.g %>%
  mutate(p = n / 4000) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    data = . %>% filter(ext.gen < 15),
    aes(
      y = vbar,
      group = interaction(n.pop0, hi.var, ndd, ext.gen),
      colour = ext.gen
    )
  ) +
  geom_ribbon(
    data = . %>% filter(ext.gen < 15),
    aes(
      ymin = vbar - 2 * sqrt(vvar / n),
      ymax = vbar + 2 * sqrt(vvar / n),
      group = interaction(n.pop0, hi.var, ndd, ext.gen),
      fill = ext.gen
    ),
    alpha = 0.1
  ) +
  geom_line(
    data = . %>% filter(ext.gen %in% 15),
    aes(
      y = vbar,
      group = interaction(n.pop0, hi.var, ndd)
    ),
    colour = 'black'
  )  +
  geom_ribbon(
    data = . %>% filter(ext.gen %in% 15),
    aes(
      ymin = vbar - 2 * sqrt(vvar / n),
      ymax = vbar + 2 * sqrt(vvar / n),
      group = interaction(n.pop0, hi.var, ndd)
    ),
    fill = 'black',
    alpha = 0.1
  ) +
  geom_rug(
    data = . %>% filter(gen == ext.gen),
    aes(
      size = p
    ),
    lineend - 'square',
    sides = 'b'
  ) +
  scale_colour_gradient(low = 'lightpink1', high = 'darkred') +
  scale_fill_gradient(low = 'lightpink1', high = 'darkred') +
  facet_wrap(ndd ~ interaction(n.pop0, hi.var), ncol = 4) +
  theme_bw()

gen.g %>%
  filter(gen == ext.gen) %>%
  mutate(p = n / 4000) %>%
  ggplot(aes(x = gen)) +
  geom_rug(
    aes(
      size = p
    ),
    sides = 'b'
  )

### Combining (with cowplot)

library(cowplot)

comp.pool.plot = all.g %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = vbar,
      group = interaction(n.pop0, hi.var, ndd),
      colour = ndd
    )
  ) +
  geom_ribbon(
    aes(
      ymin = vbar - 2 * sqrt(vvar / n),
      ymax = vbar + 2 * sqrt(vvar / n),
      group = interaction(n.pop0, hi.var, ndd),
      fill = ndd
    ),
    alpha = 0.1
  ) +
  scale_y_continuous(limits = c(0, 0.5)) +
  scale_fill_manual(values = c('black', 'purple')) +
  scale_colour_manual(values = c('black', 'purple')) +
  guides(colour = guide_legend("Growth form", nrow = 1),
         fill = guide_legend("Growth form", nrow = 1)) +
  labs(x = '', y = 'Genetic variation') +
  facet_wrap( ~ paste(hi.var, n.pop0, sep = ', '), ncol = 4) +
  theme(panel.background = element_rect(fill = 'white'),
        panel.grid = element_line(colour = 'gray88'),
        panel.border = element_rect(fill = NA),
        legend.position = 'none')

comp.pool.plot

gent.pool.plot = gen.g %>%
  ggplot(aes(x = gen)) +
  geom_line(
    data = . %>% filter(ext.gen < 15),
    aes(
      y = vbar,
      group = interaction(n.pop0, hi.var, ndd, ext.gen),
      colour = ext.gen
    )
  ) +
  geom_ribbon(
    data = . %>% filter(ext.gen < 15),
    aes(
      ymin = vbar - 2 * sqrt(vvar / n),
      ymax = vbar + 2 * sqrt(vvar / n),
      group = interaction(n.pop0, hi.var, ndd, ext.gen),
      fill = ext.gen
    ),
    alpha = 0.1
  ) +
  geom_line(
    data = . %>% filter(ext.gen %in% 15),
    aes(
      y = vbar,
      group = interaction(n.pop0, hi.var, ndd)
    ),
    colour = 'black'
  )  +
  geom_ribbon(
    data = . %>% filter(ext.gen %in% 15),
    aes(
      ymin = vbar - 2 * sqrt(vvar / n),
      ymax = vbar + 2 * sqrt(vvar / n),
      group = interaction(n.pop0, hi.var, ndd)
    ),
    fill = 'black',
    alpha = 0.1
  ) +
  scale_colour_gradient(low = 'lightpink1', high = 'darkred') +
  scale_fill_gradient(low = 'lightpink1', high = 'darkred') +
  guides(colour = guide_legend("Extinction\ngeneration"),
         fill = guide_legend('Extinction\ngeneration')) +
  labs(x = 'Generation', y = 'Genetic variation') +
  facet_wrap(ndd ~ paste(hi.var, n.pop0, sep = ', '), ncol = 4) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        panel.grid = element_line(colour = 'gray88'),
        legend.position = 'none')

gent.pool.plot

comb.plots = plot_grid(comp.pool.plot, NULL, gent.pool.plot, 
                       ncol = 1, rel_heights = c(1, -0.1, 2))

comb.plots
