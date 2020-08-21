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
  
all.summs %>%
  filter(evolve | nbar > 1) %>%
  ggplot(aes(x = gen)) +
  geom_segment(aes(x = 1, y = n.pop0,
                   xend = 15, yend = n.pop0),
               linetype = 3,
               colour = 'gray',
               size = 0.5) +
  geom_ribbon(aes(ymin = nbar - 2 * sqrt(nvar / 8000),
                  ymax = nbar + 2 * sqrt(nvar / 8000),
                  group = interaction(n.pop0, bottleneck, ndd, evolve)),
              alpha = 0.5) +
  geom_line(aes(y = nbar, 
                linetype = bottleneck,
                colour = ndd,
                group = interaction(n.pop0, bottleneck, ndd, evolve)),
            size = 2) +
  geom_point(aes(y = nbar,
                 fill = 1 - pext,
                 shape = evolve),
             shape = 21, stroke = 1,
             size = 3) +
  scale_y_log10() +
  scale_colour_manual(values = c('black', 'purple')) +
  scale_fill_viridis_c(option = 'A') +
  facet_grid(rows = vars(bottleneck), cols = vars(ndd)) +
  labs(x = 'Generation', y = 'Mean population size',
       fill = 'Prob. extant') +
  theme(panel.background = element_blank(),
        legend.position = 'bottom')


##### Comparing mean fitness, growth rates, etc.

# The best way to compare mean population sizes seems to be log-ratios.
# log(N_b / N_a) compares population sizes of b and a
# Shortcoming of this approach is that it allows only two-way comparison.

### Look at difference in population size over time

# Compare bottlenecks
all.summs %>%
  select(-c(extinct, nvar, pext)) %>%
  mutate(bottleneck = paste0('bn', bottleneck)) %>%
  spread(key = bottleneck, value = nbar) %>%
  mutate(log.bn = log(bnTRUE / bnFALSE, base = 2)) %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = log.bn,
                group = interaction(n.pop0, evolve, ndd),
                colour = evolve)) +
  facet_grid(rows = vars(n.pop0), cols = vars(ndd))

all.summs %>%
  select(-c(extinct, nvar, pext)) %>%
  mutate(bottleneck = paste0('bn', bottleneck)) %>%
  spread(key = bottleneck, value = nbar) %>%
  mutate(log.bn = log(bnTRUE / bnFALSE, base = 2)) %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = log.bn,
                group = interaction(n.pop0, evolve, ndd),
                colour = evolve,
                linetype = ndd)) +
  geom_point(aes(y = log.bn,
                 shape = factor(n.pop0)))

# Effects of bottlenecking:
#   - with no NDD and evolution, the "gap" between BN and non-BN 
#     populations grows wth time. Likely this reflects differences
#     in the max growth rates, probably due to the lost fitness!
#   - The non-evolving treatments... is this noise? Truly non-zero?
#   - Over time the differences in evolving populations with NDD
#     seems to flatten over time. Is this because NDD constrains
#     the non-BN population?
#   - There is an interaction between NDD and initial size:
#     smaller populations tend to have less penalty of bottlenecking...


# Compare initial size
all.summs %>%
  select(-c(extinct, nvar, pext)) %>%
  mutate(n.pop0 = ifelse(n.pop0 > 50, 'big', 'small')) %>%
  spread(key = n.pop0, value = nbar) %>%
  mutate(log.n0 = log(big/ small, base = 5)) %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = log.n0,
                group = interaction(bottleneck, evolve, ndd),
                colour = evolve)) +
  facet_grid(rows = vars(bottleneck), cols = vars(ndd))

all.summs %>%
  select(-c(extinct, nvar, pext)) %>%
  mutate(n.pop0 = ifelse(n.pop0 > 50, 'big', 'small')) %>%
  spread(key = n.pop0, value = nbar) %>%
  mutate(log.n0 = log(big/ small, base = 5)) %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = log.n0,
                group = interaction(bottleneck, evolve, ndd),
                colour = interaction(evolve, ndd),
                linetype = bottleneck))

# Interpretation:
# (n.b., if number is less than 1, it means that the gap between the populations
# is shrinking)
# - The gap shrinks most (and most quickly) for evolving, density dependent
# populations (likely because the density dependence acts more on the large
# population thatn the small one), although this levels out over time
# - For the evolving, the gap stays constant but then gap shrinks - initially
# these populations decline at similar rates, but then the gap shrinks... this
# means the population rebound is slower inthe large population? Think about
# this... it's counter-intuitive
# Not sure what's going on with non-evolving treatment... probably just showing
# us basic expectation of populations without NDD? NOTE: this is why the
# controls are helpful! ahhh...
# Little discernible effect of bottlenecking... except perhaps in the evolving+ndd?

