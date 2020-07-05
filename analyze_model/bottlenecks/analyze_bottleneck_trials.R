library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

setwd('analyze_model/bottlenecks')

stdrd = read.csv('bottleneck_test_pop_sizes.csv')
negdd = read.csv('bottleneck_ndd_test_pop_sizes.csv')
nevol = read.csv('bottleneck_nevolve_test_pop_sizes.csv')

ncol(stdrd)
ncol(negdd)
ncol(nevol)

stdrd = stdrd %>% mutate(alpha = 0, evolve = TRUE)
negdd = negdd %>% mutate(evolve = TRUE)
nevol = nevol %>% mutate(evolve = FALSE)

stdrd = stdrd %>% select(names(negdd))
nevol = nevol %>% select(names(nevol))

### Put it all together

all.trials = rbind(stdrd, negdd, nevol)

nrow(all.trials)

all.summs = all.trials %>%
  mutate(ndd = alpha > 0) %>%
  distinct(n.pop0, bottleneck, evolve, ndd, gen, .keep_all = TRUE) %>%
  select(n.pop0, bottleneck, evolve, ndd, extinct,
         gen, nbar, nvar, pext)

all.summs %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = nbar, 
                group = interaction(n.pop0, bottleneck, evolve, ndd))) +
  scale_y_log10()

all.summs %>%
  filter(evolve) %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = nbar, 
                group = interaction(n.pop0, bottleneck, ndd))) +
  facet_grid(rows = vars(n.pop0), cols = vars(bottleneck)) +
  scale_y_log10()

all.summs %>%
  filter(evolve) %>%
  ggplot(aes(x = gen)) +
  geom_segment(aes(x = 1, y = n.pop0,
                   xend = 15, yend = n.pop0),
               linetype = 3,
               colour = 'gray',
               size = 0.5) +
  geom_ribbon(aes(ymin = nbar - 2 * sqrt(nvar / 8000),
                  ymax = nbar + 2 * sqrt(nvar / 8000),
                  group = interaction(n.pop0, bottleneck, ndd)),
              alpha = 0.5) +
  geom_line(aes(y = nbar, 
                linetype = bottleneck,
                colour = ndd,
                group = interaction(n.pop0, bottleneck, ndd)),
            size = 2) +
  geom_point(aes(y = nbar,
                 fill = 1 - pext),
             shape = 21, stroke = 1,
             size = 3) +
  scale_y_log10() +
  scale_colour_manual(values = c('black', 'purple')) +
  scale_fill_viridis_c(option = 'A') +
  facet_grid(cols = vars(n.pop0)) +
  labs(x = 'Generation', y = 'Mean population size',
       fill = 'Prob. extant') +
  theme(panel.background = element_blank(),
        legend.position = 'bottom') +
  ggsave('2k_fourway_fig_2020-07-03.pdf')
  
