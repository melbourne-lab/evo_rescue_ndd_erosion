library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

yes = data.frame(
  x = rep(0:9, each = 2) + rep(1:2, times = 10),
  y = rep(0:4, each = 4) + rep(1:2, times = 10)
)

plot(yes, type = 'l')

yes %>%
  ggplot(aes(x = x, y = y)) +
  geom_line()

# Try for rizzle

all.data = read.csv('simulations/outputs/final_results/alldata_combined.csv')

# Want to plot two populatons.

test.two = all.data %>% 
  filter(trial %in% 1:2, !alpha, n.pop0 < 50, low.var) %>%
  rbind(., .) %>%
  arrange(trial, gen) %>%
  group_by(trial) %>%
  mutate(v2 = c(NA, v[1:(n()-1)])) %>%
  ungroup() %>%
  filter(!is.na(v2))
  #ungroup()

head(test.two)

test.two %>%
  ggplot() +
  geom_path(aes(x = n, y = v2, group = trial)) +
  geom_point(aes(x = n, y = v2, colour = factor(gen), shape = factor(trial)), size = 3) +
  scale_color_brewer('BuPu')

###

steps.all = all.data %>% 
  rbind(., .) %>%
  arrange(trial, alpha, low.var, n.pop0, gen) %>%
  group_by(trial, alpha, low.var, n.pop0) %>%
  mutate(v2 = c(NA, v[1:(n()-1)])) %>%
  ungroup() %>%
  filter(!is.na(v2))

steps.all %>% filter(trial %in% 1:2, !alpha, n.pop0 < 50, low.var)

steps.all %>%
  filter(trial %in% 1:100) %>%
  ggplot(aes(x = n, y = v2)) +
  geom_path(aes(colour = extinct, group = trial), alpha = 0.5) +
  scale_x_log10() +
  scale_color_manual(values = c('black', 'red')) +
  facet_wrap(n.pop0 + low.var ~ alpha, nrow = 2)

### Try with summary

all.summ = all.data %>%
  group_by(n.pop0, low.var, alpha, trial) %>%
  arrange(desc(gen)) %>%
  distinct(n.pop0, low.var, alpha, trial, .keep_all = TRUE) %>%
  uncount(weight = 15 - gen) %>%
  mutate(n = 0, gen = gen + (1:(15-gen[1]))) %>%
  ungroup() %>% 
  rbind(all.data) %>%
  group_by(alpha, low.var, n.pop0, extinct, gen) %>%
  summarise(n.bar = mean(n),
            n.var = var(n),
            v.bar = mean(v),
            v.var = var(v),
            w.bar = mean(wbar),
            w.var = var(wbar),
            nobs = n()
  ) %>%
  ungroup() %>%
  filter(!(extinct & gen %in% 15))

all.summ %>% # filter(!low.var) %>%
  ggplot(aes(x = n.bar, y = v.bar)) +
  geom_path(
    aes(
      group = interaction(alpha, extinct), 
      linetype = extinct, 
      colour = factor(alpha)
      )
    ) +
  geom_point(
    aes(
      colour = factor(alpha)
    ),
    size = 5,
    shape = 21,
    fill = 'white'
  ) +
  geom_text(
    aes(
      colour = factor(alpha),
      label = gen
    )
  ) +
  # geom_ribbon(
  #   aes(
  #     group = interaction(alpha, extinct), 
  #     linetype = extinct, 
  #     colour = factor(alpha)
  #   )
  # ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_x_log10() +
  facet_wrap(n.pop0 ~ low.var, ncol = 4)

all.summ %>%
  ggplot(aes(x = n.bar, y = w.bar)) +
  geom_path(
    aes(
      group = interaction(alpha, extinct), 
      linetype = extinct, 
      colour = factor(alpha)
    )
  ) +
  geom_point(
    aes(
      colour = factor(alpha)
    ),
    size = 5,
    shape = 21,
    fill = 'white'
  ) +
  geom_text(
    aes(
      colour = factor(alpha),
      label = gen
    )
  ) +
  # geom_ribbon(
  #   aes(
  #     group = interaction(alpha, extinct), 
  #     linetype = extinct, 
  #     colour = factor(alpha)
  #   )
  # ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_x_log10() +
  facet_wrap(n.pop0 ~ low.var, ncol = 4)

all.summ %>%
  ggplot(aes(x = v.bar, y = w.bar)) +
  geom_path(
    aes(
      group = interaction(alpha, extinct), 
      linetype = extinct, 
      colour = factor(alpha)
    )
  ) +
  geom_point(
    aes(
      colour = factor(alpha)
    ),
    size = 5,
    shape = 21,
    fill = 'white'
  ) +
  geom_text(
    aes(
      colour = factor(alpha),
      label = gen
    )
  ) +
  # geom_ribbon(
  #   aes(
  #     group = interaction(alpha, extinct), 
  #     linetype = extinct, 
  #     colour = factor(alpha)
  #   )
  # ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_x_log10() +
  facet_wrap(n.pop0 ~ low.var, ncol = 4)

### Trying one plot with confints...

all.summ %>%
  filter(tau > 0) %>%
  ggplot() +
  geom_path(
    aes(
      x = v.bar,
      y = w.bar,
      group = interaction(alpha, extinct), 
      linetype = extinct, 
      colour = factor(alpha)
    )
  ) +
  geom_path(
    aes(
      x = v.bar - 2 * sqrt(v.var / sqrt(nobs)),
      y = w.bar - 2 * sqrt(w.var / sqrt(nobs)),
      group = interaction(alpha, extinct),
      colour = factor(alpha)
    ),
    linetype = 3,
    size = 0.4
  ) +
  geom_path(
    aes(
      x = v.bar + 2 * sqrt(v.var / sqrt(nobs)),
      y = w.bar + 2 * sqrt(w.var / sqrt(nobs)),
      group = interaction(alpha, extinct),
      colour = factor(alpha)
    ),
    linetype = 3,
    size = 0.4
  ) +
  geom_point( 
    aes(
        x = v.bar,
        y = w.bar,
        colour = factor(alpha)
    )
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  facet_wrap(n.pop0 ~ low.var, ncol = 4)

### Try putting everything together.

library(cowplot)

plot.nv = all.summ %>%
  ggplot(aes(x = n.bar, y = v.bar)) +
  geom_path(
    aes(
      group = interaction(alpha, extinct), 
      linetype = extinct, 
      colour = factor(alpha)
    )
  ) +
  geom_point(
    aes(
      colour = factor(alpha)
    )
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_x_log10() +
  facet_wrap(low.var ~ n.pop0, ncol = 4) +
  theme(legend.position = 'none')

plot.nw = all.summ %>%
  ggplot(aes(x = n.bar, y = w.bar)) +
  geom_path(
    aes(
      group = interaction(alpha, extinct), 
      linetype = extinct, 
      colour = factor(alpha)
    )
  ) +
  geom_point(
    aes(
      colour = factor(alpha)
    )
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_x_log10() +
  facet_wrap(low.var ~ n.pop0, ncol = 4) +
  theme(legend.position = 'none')

plot.vw = all.summ %>%
  ggplot(aes(x = v.bar, y = w.bar)) +
  geom_path(
    aes(
      group = interaction(alpha, extinct), 
      linetype = extinct, 
      colour = factor(alpha)
    )
  ) +
  geom_point(
    aes(
      colour = factor(alpha)
    )
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  facet_wrap(low.var ~ n.pop0, ncol = 4) +
  theme(legend.position = 'none')

plot_grid(plot.nv, plot.nw, plot.vw, 
          nrow = 3, rel_heights = c(1, 1, 1))

### Try Pr(extinct | x) plots

binned.pext = all.data %>%
  mutate(n.round = ifelse(n > 100, 101, n),
         v.round = round(v / .05) * .05,
         w.round = round(wbar / .05) * .05) %>%
  group_by(n.pop0, low.var, alpha, n.round, v.round, w.round) %>%
  summarise(pext = mean(extinct),
            n = n())

binned.nv = all.data %>%
  mutate(n.round = ifelse(n > 100, 101, 5 * round(n/5)),
         v.round = round(v / .05) * .05) %>%
  group_by(n.pop0, low.var, alpha, n.round, v.round) %>%
  summarise(pext = mean(extinct),
            n = n()) %>%
  uncount(weights = 5) %>%
  group_by(n.pop0, low.var, alpha, v.round) %>%
  mutate(n.round = n.round + 0:4)

binned.vw = all.data %>%
  mutate(v.round = round(v / .05) * .05,
         w.round = round(wbar / .05) * .05) %>%
  group_by(n.pop0, low.var, alpha, v.round, w.round) %>%
  summarise(pext = mean(extinct),
            n = n())

binned.nw = all.data %>%
  mutate(n.round = ifelse(n > 100, 101, 5 * round(n/5)),
         w.round = round(wbar / .05) * .05) %>%
  group_by(n.pop0, low.var, alpha, n.round, w.round) %>%
  summarise(pext = mean(extinct),
            n = n()) %>%
  uncount(weights = 5) %>%
  group_by(n.pop0, low.var, alpha, w.round) %>%
  mutate(n.round = n.round + 0:4)

binned.nv %>%
  ggplot() +
  geom_raster(
    aes(
      x = v.round,
      y = n.round,
      fill = pext #> 0.5
    )
  ) +
  scale_fill_viridis_b(option = 'B') +
  #scale_fill_manual(values = c('black', 'red')) +
  facet_wrap(paste(n.pop0, low.var, sep = ', ') ~ alpha, nrow = 4)

binned.nw %>%
  ggplot() +
  geom_raster(
    aes(
      x = w.round,
      y = n.round,
      fill = pext #> 0.5
    )
  ) +
  scale_fill_viridis_b(option = 'B') +
  #scale_fill_manual(values = c('black', 'red')) +
  facet_wrap(paste(n.pop0, low.var, sep = ', ') ~ alpha, nrow = 4)

binned.vw %>%
  ggplot() +
  geom_raster(
    aes(
      x = w.round,
      y = v.round,
      fill = pext #> 0.5
    )
  ) +
  scale_fill_viridis_b(option = 'B') +
  #scale_fill_manual(values = c('black', 'red')) +
  facet_wrap(paste(n.pop0, low.var, sep = ', ') ~ alpha, nrow = 4)

### What would capture hte vortex more is... arrows!

# First need to add in ones:

all.data.ones = all.data %>%
  select(n.pop0, low.var, alpha, trial, gen, n, v, wbar, extinct) %>%
  group_by(n.pop0, low.var, alpha, trial) %>%
  arrange(desc(gen)) %>%
  distinct(n.pop0, low.var, alpha, trial, .keep_all = TRUE) %>%
  filter(extinct & (n > 1)) %>%
  mutate_at(vars(wbar, v), list(~ NA)) %>%
  mutate(n = 1, gen = gen + 1) %>%
  ungroup() %>%
  rbind(all.data %>% select(n.pop0, low.var, alpha, trial, gen, n, v, wbar, extinct)) %>%
  arrange(n.pop0, low.var, alpha, trial, gen)

all.deltas = all.data.ones %>%
  group_by(n.pop0, low.var, alpha, trial) %>%
  mutate(n_tp1 = c(n[2:n()], NA),
         v_tp1 = c(v[2:n()], NA),
         w_tp1 = c(wbar[2:n()], NA))# c(exp(diff(log(v))), NA)
# this takes a long time to run... never seen anything like this.
# takes literally minutes!!! Wht the fuck
# there must be a faster way... oh well.

nv.slopes = all.deltas %>%
  mutate(n.round = ifelse(n > 100, 101, 5 * round(n / 5)),
         v.round = round(v / .05) * .05) %>%
  group_by(n.pop0, low.var, alpha, n.round, v.round) %>%
  summarise(ntp1bar = mean(n_tp1, na.rm = TRUE),
            vtp1bar = mean(v_tp1, na.rm = TRUE)) 

nw.slopes = all.deltas %>%
  mutate(n.round = ifelse(n > 100, 101, 5 * round(n / 5)),
         w.round = round(wbar / .05) * .05) %>%
  group_by(n.pop0, low.var, alpha, n.round, w.round) %>%
  summarise(ntp1bar = mean(n_tp1, na.rm = TRUE),
            wtp1bar = mean(w_tp1, na.rm = TRUE)) 

vw.slopes = all.deltas %>%
  mutate(v.round = round(v / .05) * .05,
         w.round = round(wbar / .05) * .05) %>%
  group_by(n.pop0, low.var, alpha, v.round, w.round) %>%
  summarise(vtp1bar = mean(v_tp1, na.rm = TRUE),
            wtp1bar = mean(w_tp1, na.rm = TRUE)) 

nv.slopes %>%
  ggplot() +
  geom_segment(aes(x = n.round, xend = ntp1bar,
                   y = v.round, yend = vtp1bar,
                   colour = ntp1bar > n.round),
               size = 0.8, 
               arrow = arrow(length = unit(0.01, "npc"))) +
  scale_color_manual(values = c('blue', 'red')) +
  scale_x_continuous(limits = c(0, 100)) +
  facet_wrap(n.pop0 + low.var ~ alpha, ncol = 4) +
  ggsave('~/Dropbox/rescue_ndd_paper_2020/figures/ex_n-v_slopes.pdf')

nw.slopes %>%
  ggplot() +
  geom_segment(aes(x = n.round, xend = ntp1bar,
                   y = w.round, yend = wtp1bar,
                   colour = ntp1bar > n.round),
               size = 0.8, 
               arrow = arrow(length = unit(0.01, "npc"))) +
  scale_color_manual(values = c('blue', 'red')) +
  scale_x_continuous(limits = c(0, 100)) +
  facet_wrap(n.pop0 + low.var ~ alpha, ncol = 4) +
  ggsave('~/Dropbox/rescue_ndd_paper_2020/figures/ex_n-w_slopes.pdf')
  
vw.slopes %>%
  ggplot() +
  geom_segment(aes(x = v.round, xend = vtp1bar,
                   y = w.round, yend = wtp1bar),
               size = 0.8, 
               arrow = arrow(length = unit(0.01, "npc"))) +
  facet_wrap(n.pop0 + low.var ~ alpha, ncol = 4)

### Okay try actual semi-serious plots in cowplot

nv.hivar = all.summ %>%  
  filter(!low.var) %>%
  mutate(n.pop0 = factor(n.pop0, labels = c("Small", "Large"))) %>%
  ggplot(aes(x = n.bar, y = v.bar)) +
  geom_path(
    aes(
      group = interaction(alpha, extinct), 
      linetype = extinct, 
      colour = factor(alpha)
    )
  ) +
  geom_point(
    aes(
      colour = factor(alpha),
      shape = extinct
    ),
    size = 4.5,
    fill = 'white'
  ) +
  geom_text(
    aes(
      colour = factor(alpha),
      label = gen
    ),
    size = 3
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_shape_manual(values = 21:22) +
  scale_x_log10(breaks = c(0.1, 10, 100), labels = c(0.1, 10, 1000)) +
  labs(x = 'Population size', y = 'Genetic variance') +
  facet_wrap(n.pop0 ~ ., ncol = 2) +
  theme(legend.position = 'none',
        panel.background = element_blank())

nw.hivar = all.summ %>%  
  filter(!low.var) %>%
  mutate(n.pop0 = factor(n.pop0, labels = c("Small", "Large"))) %>%
  ggplot(aes(x = n.bar, y = w.bar)) +
  geom_path(
    aes(
      group = interaction(alpha, extinct), 
      linetype = extinct, 
      colour = factor(alpha)
    )
  ) +
  geom_point(
    aes(
      colour = factor(alpha),
      shape = extinct
    ),
    size = 4.5,
    fill = 'white'
  ) +
  geom_text(
    aes(
      colour = factor(alpha),
      label = gen
    ),
    size = 3
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_shape_manual(values = 21:22) +
  scale_x_log10(breaks = c(0.1, 10, 100), labels = c(0.1, 10, 1000)) +
  labs(x = 'Population size', y = 'Intrinsic fitness') +
  facet_wrap(n.pop0 ~ ., ncol = 2) +
  theme(legend.position = 'none',
        panel.background = element_blank())

vw.hivar = all.summ %>%  
  filter(!low.var) %>%
  mutate(n.pop0 = factor(n.pop0, labels = c("Small", "Large"))) %>%
  ggplot(aes(x = v.bar, y = w.bar)) +
  geom_path(
    aes(
      group = interaction(alpha, extinct), 
      linetype = extinct, 
      colour = factor(alpha)
    )
  ) +
  geom_point(
    aes(
      colour = factor(alpha),
      shape = extinct
    ),
    size = 4.5,
    fill = 'white'
  ) +
  geom_text(
    aes(
      colour = factor(alpha),
      label = gen
    ),
    size = 3
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_shape_manual(values = 21:22) +
  labs(x = 'Additive genetic variance', y = 'Intrinsic fitness') +
  facet_wrap(n.pop0 ~ ., ncol = 2) +
  theme(legend.position = 'none',
        panel.background = element_blank())

save_plot('~/Dropbox/rescue_ndd_paper_2020/figures/eg_vortex_hivar.pdf',
          plot_grid(nv.hivar, nw.hivar, vw.hivar, nrow = 3),
          base_height = 11, base_width = 8)

# Low variance vortex

nv.lovar = all.summ %>%  
  filter(low.var) %>%
  mutate(n.pop0 = factor(n.pop0, labels = c("Small", "Large"))) %>%
  ggplot(aes(x = n.bar, y = v.bar)) +
  geom_path(
    aes(
      group = interaction(alpha, extinct), 
      linetype = extinct, 
      colour = factor(alpha)
    )
  ) +
  geom_point(
    aes(
      colour = factor(alpha),
      shape = extinct
    ),
    size = 4.5,
    fill = 'white'
  ) +
  geom_text(
    aes(
      colour = factor(alpha),
      label = gen
    ),
    size = 3
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_shape_manual(values = 21:22) +
  scale_x_log10(breaks = c(0.1, 10, 100), labels = c(0.1, 10, 1000)) +
  labs(x = 'Population size', y = 'Genetic variance') +
  facet_wrap(n.pop0 ~ ., ncol = 2) +
  theme(legend.position = 'none',
        panel.background = element_blank())

nw.lovar = all.summ %>%  
  filter(low.var) %>%
  mutate(n.pop0 = factor(n.pop0, labels = c("Small", "Large"))) %>%
  ggplot(aes(x = n.bar, y = w.bar)) +
  geom_path(
    aes(
      group = interaction(alpha, extinct), 
      linetype = extinct, 
      colour = factor(alpha)
    )
  ) +
  geom_point(
    aes(
      colour = factor(alpha),
      shape = extinct
    ),
    size = 4.5,
    fill = 'white'
  ) +
  geom_text(
    aes(
      colour = factor(alpha),
      label = gen
    ),
    size = 3
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_shape_manual(values = 21:22) +
  scale_x_log10(breaks = c(0.1, 10, 100), labels = c(0.1, 10, 1000)) +
  labs(x = 'Population size', y = 'Intrinsic fitness') +
  facet_wrap(n.pop0 ~ ., ncol = 2) +
  theme(legend.position = 'none',
        panel.background = element_blank())

vw.lovar = all.summ %>%  
  filter(low.var) %>%
  mutate(n.pop0 = factor(n.pop0, labels = c("Small", "Large"))) %>%
  ggplot(aes(x = v.bar, y = w.bar)) +
  geom_path(
    aes(
      group = interaction(alpha, extinct), 
      linetype = extinct, 
      colour = factor(alpha)
    )
  ) +
  geom_point(
    aes(
      colour = factor(alpha),
      shape = extinct
    ),
    size = 4.5,
    fill = 'white'
  ) +
  geom_text(
    aes(
      colour = factor(alpha),
      label = gen
    ),
    size = 3
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_shape_manual(values = 21:22) +
  labs(x = 'Additive genetic variance', y = 'Intrinsic fitness') +
  facet_wrap(n.pop0 ~ ., ncol = 2) +
  theme(legend.position = 'none',
        panel.background = element_blank())

save_plot('~/Dropbox/rescue_ndd_paper_2020/figures/eg_vortex_lowvar.pdf',
          plot_grid(nv.lovar, nw.lovar, vw.lovar, nrow = 3),
          base_height = 11, base_width = 8)

### Try to get slope fields

nv.slopes = all.deltas %>%
  mutate(n.round = ifelse(n > 100, 101, 5 * round(n / 5)),
         v.round = round(v / .05) * .05) %>%
  group_by(n.pop0, low.var, alpha, n.round, v.round) %>%
  summarise(nbar = mean(n),
            vbar = mean(v),
            ntp1bar = mean(n_tp1, na.rm = TRUE),
            vtp1bar = mean(v_tp1, na.rm = TRUE)) 

nw.slopes = all.deltas %>%
  mutate(n.round = ifelse(n > 100, 101, 5 * round(n / 5)),
         w.round = round(wbar / .1) * .1) %>%
  group_by(n.pop0, low.var, alpha, n.round, w.round) %>%
  summarise(nbar = mean(n),
            wbar = mean(wbar),
            ntp1bar = mean(n_tp1, na.rm = TRUE),
            wtp1bar = mean(w_tp1, na.rm = TRUE)) 

nv.slope.hivar = nv.slopes %>%
  ungroup() %>%
  filter(!low.var) %>%
  mutate(n.pop0 = factor(n.pop0, labels = c("Large", "Small")),
         ndd = ifelse(alpha > 0, "dens. dep.", "dens. dep.")) %>%
  ggplot() +
  geom_segment(aes(x = nbar, xend = ntp1bar,
                   y = vbar, yend = vtp1bar),
               size = 0.8, 
               arrow = arrow(length = unit(0.02, "npc"))) +
  scale_x_continuous(limits = c(0, 100), oob = scales::squish) +
  facet_wrap(n.pop0 ~ alpha, ncol = 2) 

nw.slope.hivar = nw.slopes %>%
  ungroup() %>%
  filter(!low.var) %>%
  mutate(n.pop0 = factor(n.pop0, labels = c("Large", "Small"))) %>%
  ggplot() +
  geom_segment(aes(x = nbar, xend = ntp1bar,
                   y = wbar, yend = wtp1bar),
               size = 0.4, 
               arrow = arrow(length = unit(0.02, "npc"))) +
  scale_x_continuous(limits = c(0, 100), oob = scales::squish) +
  facet_wrap(n.pop0 ~ alpha, ncol = 2)  

plot_grid(nv.slope.hivar, nw.slope.hivar, nrow = 2)

nv.slope.hivar = nv.slopes %>%
  ungroup() %>%
  filter(low.var) %>%
  mutate(n.pop0 = factor(n.pop0, labels = c("Large", "Small"))) %>%
  ggplot() +
  geom_segment(aes(x = nbar, xend = ntp1bar,
                   y = vbar, yend = vtp1bar),
               size = 0.8, 
               arrow = arrow(length = unit(0.02, "npc"))) +
  scale_x_continuous(limits = c(0, 100), oob = scales::squish) +
  facet_wrap(n.pop0 ~ alpha, ncol = 2) 

nw.slope.hivar = nw.slopes %>%
  ungroup() %>%
  filter(low.var) %>%
  mutate(n.pop0 = factor(n.pop0, labels = c("Large", "Small"))) %>%
  ggplot() +
  geom_segment(aes(x = nbar, xend = ntp1bar,
                   y = wbar, yend = wtp1bar),
               size = 0.4, 
               arrow = arrow(length = unit(0.02, "npc"))) +
  scale_x_continuous(limits = c(0, 100), oob = scales::squish) +
  facet_wrap(n.pop0 ~ alpha, ncol = 2) 

plot_grid(nv.slope.hivar, nw.slope.hivar, nrow = 2)
