library(ggplot2)
library(dplyr)
library(tidyr)

bns = read.csv('analyze_model/bottleneck_test_pop_sizes.csv')

head(bns)
nrow(bns)

n.trials = max(bns$trial)

bns %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = n, group = trial, colour = extinct), 
            size = 0.025) +
  geom_line(aes(y = nbar, group = interaction(n.pop0, bottleneck)),
            size = 2) +
  theme(panel.background = element_blank()) +
  scale_colour_manual(values = c('gray', 'red')) +
  scale_y_log10() +
  facet_wrap(factor(n.pop0) ~ bottleneck)

bns %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = nbar, group = interaction(n.pop0, bottleneck),
                linetype = bottleneck)) +
  geom_point(aes(y = nbar, fill = pext),
             shape = 21, stroke = 1, size = 3) +
  theme(panel.background = element_blank()) +
  scale_colour_manual(values = c('gray', 'red')) +
  scale_fill_gradient(low = 'red', high = 'black') +
  scale_y_log10()
  
bns %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = nbar, group = interaction(n.pop0, bottleneck),
                linetype = bottleneck)) +
  geom_ribbon(aes(ymin = nbar - 2 * sqrt(nvar / n.trials),
                  ymax = nbar + 2 * sqrt(nvar / n.trials),
                  group = interaction(n.pop0, bottleneck)),
               alpha = 0.2) +
  geom_point(aes(y = nbar, fill = 1-pext),
             stroke = 1, shape = 21, size = 3) +
  geom_segment(aes(x = 1, xend = 15,
                   y = n.pop0,
                   yend = n.pop0),
               linetype = 2,
               size = 0.1) +
  theme(panel.background = element_blank(),
        legend.position = 'bottom') +
  scale_fill_viridis_c(option = 'B') +
  scale_y_log10() +
  guides(linetype = guide_legend("Bottleneck", title.position = 'top'),
         fill = guide_legend("Prob. extinction", title.position = 'top')) +
  facet_wrap(~ factor(n.pop0))
