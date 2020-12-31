# Quick and dirty plot looking at differences in genotypic change between
# density dependent and density independent populations.
# I ran 1100 simulations of each alpha = 0 and alpha > 0 trials at two different
# levels of initial genetic variation.
# Note that to minimize extinctions all populations are initialized at N = 120.
# At each time step I saved genotype (g_i), phenotype (z_i), and fitness (w_i).
# Here I am plotting results.
# Simulation code: simulations/sim_thrice_bottlenecked/sim_a***_gt_big_init.R


### Load packages

library(ggplot2)
library(dplyr)
library(tidyr)

### Load and combine data

# density independent trials
bign.a000 = read.csv('simulations/outputs/thrice_bottlenecked_results/a000_gt_big_init.csv')
# density dependent trials
bign.a035 = read.csv('simulations/outputs/thrice_bottlenecked_results/a035_gt_big_init.csv')

# combine
bign0.gt = rbind(
  bign.a000 %>% mutate(ndd = FALSE),
  bign.a035 %>% mutate(ndd = TRUE)
) %>%
  # Add column for probablity/proportion of populations extant at this time
  mutate(p.extant = n.trials / max(n.trials))

### Plot results

bign0.gt %>%
  mutate(p.extant = n.trials / max(n.trials)) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(y = gbar, group = interaction(bottleneck, ndd))
  ) +
  geom_ribbon(
    aes(ymin = gbar - 2 * sqrt(gvar / n.trials),
        ymax = gbar + 2 * sqrt(gvar / n.trials),
        group = interaction(bottleneck, ndd)),
    alpha = 0.2
  ) +
  geom_point(
    aes(y = gbar,
        fill = ndd,
        alpha = p.extant),
    shape = 22, size = 3
  ) +
  scale_alpha_continuous(range = c(0, 1))
# Here, plotting genotypes (g_i) over time
# Colored by density dependence (NDD = TRUE is density dependent)

bign0.gt %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(y = wbar, 
        linetype = ndd,
        group = interaction(bottleneck, ndd)),
    size = 2
  ) +
  geom_ribbon(
    aes(ymin = wbar - 2 * sqrt(wvar / n.trials),
        ymax = wbar + 2 * sqrt(wvar / n.trials),
        group = interaction(bottleneck, ndd)),
    alpha = 0.2
  ) +
  geom_point(
    aes(y = wbar,
        fill = 1 - p.extant),
    shape = 22, size = 4
  ) +
  geom_segment(
    aes(x = 0, xend = 10, y = 1, yend = 1),
    linetype = 3, colour = 'gray22'
  ) +
  scale_fill_viridis_c(option = 'B') +
  theme(panel.background = element_blank(),
        legend.position = 'bottom') +
  labs(x = 'Generation', y = 'Mean intrisic fitness',
       linetype = 'Density\ndependence',
       fill = 'Extinction\nprobability')
# Here, plotting fitness (w, wbar, w_i)
# Darker color = fewer populations extinct (higher probability of a population being extant).

bign0.gt %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(y = gbar, 
        linetype = ndd,
        group = interaction(bottleneck, ndd)),
    size = 2
  ) +
  geom_ribbon(
    aes(ymin = gbar - 2 * sqrt(gvar / n.trials),
        ymax = gbar + 2 * sqrt(gvar / n.trials),
        group = interaction(bottleneck, ndd)),
    alpha = 0.2
  ) +
  geom_point(
    aes(y = gbar,
        fill = 1 - p.extant),
    shape = 22, size = 4
  ) +
  geom_text(
    x = 2.5, y = 0.5, label = 'High variation'
  ) +
  geom_text(
    x = 8, y = 0.25, label = 'Low variation'
  ) +
  scale_x_continuous(breaks = 2*(1:5)) +
  scale_fill_viridis_c(option = 'B') +
  labs(x = 'Generation', y = 'Mean genotype',
       linetype = 'Density\ndependence',
       fill = 'Extinction\nprobability') +
  theme(panel.background = element_blank(),
        legend.position = 'bottom')
### Conclusions

# Here, actually, there is no difference in rate of genotypic change.
# Huh... think this is because populations never get small enough to lose
# genetic variation.
# You can definitely see the effect of genetic variation.
# The difference in extinctions (need a better way to visualize this) shows NDD
# selecting for more fit populations.

# Differences in genetic variation are caused by populations being small. But
# then those populations go extinct. This raises the average genotype of the
# surviving population, partially negating the effect of eroded genetic
# variation.

