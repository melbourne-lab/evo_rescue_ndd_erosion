library(ggplot2)
library(dplyr)
library(tidyr)

all.data = read.csv('simulations/outputs/final_results/alldata_combined.csv')

# Get extinction probs
pext = all.data %>%
  distinct(trial, n.pop0, alpha, low.var, extinct) %>%
  group_by(n.pop0, alpha, low.var) %>%
  summarise(p.ext = mean(extinct))

p.ext.v = all.data %>%
  mutate(vr = round(v * 50) / 50) %>%
  group_by(n.pop0, low.var, alpha, gen, extinct) %>%
  mutate(n.ex = n()) %>%
  group_by(n.pop0, low.var, alpha, gen, extinct, vr) %>%
  summarise(n.vr.ex = n(),
            n.ex = n.ex[1],
            p.vr.ex = n.vr.ex / n.ex) %>%
  merge(y = pext) %>%
  group_by(n.pop0, low.var, alpha, gen, vr) %>%
  filter(any(extinct) & any(!extinct)) %>%
  summarise(p1 = (p.vr.ex[extinct] * p.ext[1]),
            p2 = (p.vr.ex[!extinct] * (1 - p.ext[1])),
            p.ex.vr = p1 / (p1 + p2),
            n.ex = n.vr.ex[extinct],
            n.sv = n.vr.ex[!extinct])

head(test2)  

ggplot(test2 %>% filter(gen %in% 1:10)) +
  geom_tile(
    aes(
      x = gen, y = vr, fill = p.ex.vr
    )
  ) +
  scale_fill_viridis_b(option = 'B') +
  facet_wrap(n.pop0 + low.var ~ alpha, ncol = 4)

ggplot(p.ext.v %>% filter(gen %in% (1 + (0:3)*2), n.ex > 10 & n.sv > 10)) +
  geom_line(
    aes(
      x = vr, y = p.ex.vr,
      group = interaction(gen, alpha),
      colour = factor(alpha)
    )
  ) +
  facet_wrap(n.pop0 + low.var ~ gen, nrow = 4)

#### By intrinsic fitness w

p.ext.w = all.data %>%
  mutate(wr = round(wbar * 10) / 10) %>%
  group_by(n.pop0, low.var, alpha, gen, extinct) %>%
  mutate(n.ex = n()) %>%
  group_by(n.pop0, low.var, alpha, gen, extinct, wr) %>%
  summarise(n.w.e = n(),
            n.e = n.ex[1],
            p.w.e = n.w.e / n.e) %>%
  merge(y = pext) %>%  
  group_by(n.pop0, low.var, alpha, gen, wr) %>%
  filter(any(extinct) & any(!extinct)) %>%
  summarise(p1 = (p.w.e[extinct] * p.ext[1]),
            p2 = (p.w.e[!extinct] * (1 - p.ext[1])),
            p.e.w = p1 / (p1 + p2),
            n.ex = n.w.e[extinct],
            n.sv = n.w.e[!extinct])

ggplot(p.ext.w %>% filter(gen %in% 1:6)) +
  geom_tile(
    aes(
      x = gen, y = wr, fill = p.e.w
    )
  ) +
  scale_fill_viridis_b(option = 'B') +
  facet_wrap(n.pop0 + low.var ~ alpha, ncol = 4)

ggplot(p.ext.w %>% filter(gen %in% (1 + (0:3)*2), n.ex > 10 & n.sv > 10)) +
  geom_line(
    aes(
      x = wr, y = p.e.w,
      group = interaction(gen, alpha),
      colour = factor(alpha)
    )
  ) +
  facet_wrap(n.pop0 + low.var ~ gen, nrow = 4)

ggplot(p.ext.w %>% filter(gen %in% (1 + (0:2)*3), n.ex > 10 & n.sv > 10),
       aes(x = wr, y = p.e.w)) +
  geom_line(
    aes(
      group = interaction(gen, alpha),
      colour = factor(alpha)
    ),
    size = 0.75
  ) +
  geom_point(
    aes(
      colour = factor(alpha),
      shape = factor(gen)
    ),
    size = 2
  ) +
  scale_shape_manual(values = c(16, 17, 15)) +
  scale_color_manual(values = c('black', 'purple')) +
  facet_wrap(n.pop0 ~ low.var, nrow = 2)

ggplot(p.ext.w %>% filter(gen %in% (1 + (0:3)*2), n.ex > 10 & n.sv > 10)) +
  geom_line(
    aes(
      x = wr, y = p.e.w,
      group = interaction(gen, alpha),
      colour = factor(alpha)
    )
  ) +
  facet_wrap(n.pop0 + low.var ~ gen, nrow = 4)

# ---

p.ext.g = all.data %>%
  mutate(gr = round(gbar * 10) / 10) %>%
  group_by(n.pop0, low.var, alpha, gen, extinct) %>%
  mutate(n.ex = n()) %>%
  group_by(n.pop0, low.var, alpha, gen, extinct, gr) %>%
  summarise(n.g.e = n(),
            n.e = n.ex[1],
            p.g.e = n.g.e / n.e) %>%
  merge(y = pext) %>%  
  group_by(n.pop0, low.var, alpha, gen, gr) %>%
  filter(any(extinct) & any(!extinct)) %>%
  summarise(p1 = (p.g.e[extinct] * p.ext[1]),
            p2 = (p.g.e[!extinct] * (1 - p.ext[1])),
            p.g.w = p1 / (p1 + p2),
            n.ex = n.g.e[extinct],
            n.sv = n.g.e[!extinct])

ggplot(p.ext.g %>% filter(gen %in% (1 + (0:2)*3), n.ex > 10 & n.sv > 10),
       aes(x = gr, y = p.g.w)) +
  geom_line(
    aes(
      group = interaction(gen, alpha),
      colour = factor(alpha)
    ),
    size = 0.75
  ) +
  geom_point(
    aes(
      colour = factor(alpha),
      shape = factor(gen)
    ),
    size = 2
  ) +
  scale_shape_manual(values = c(16, 17, 15)) +
  scale_color_manual(values = c('black', 'purple')) +
  facet_wrap(n.pop0 ~ low.var, nrow = 2)

ggplot(p.ext.g %>% filter(gen %in% (1 + (0:3)*2), n.ex > 10 & n.sv > 10)) +
  geom_line(
    aes(
      x = gr, y = p.g.w,
      group = interaction(gen, alpha),
      colour = factor(alpha)
    )
  ) +
  facet_wrap(n.pop0 + low.var ~ gen, nrow = 4)

p.ext.g.overall = all.data %>%
  mutate(gr = round(gbar * 10) / 10) %>%
  group_by(n.pop0, low.var, alpha, extinct) %>%
  mutate(n.ex = n()) %>%
  group_by(n.pop0, low.var, alpha, extinct, gr) %>%
  summarise(n.g.e = n(),
            n.e = n.ex[1],
            p.g.e = n.g.e / n.e) %>%
  merge(y = pext) %>%  
  group_by(n.pop0, low.var, alpha, gr) %>%
  filter(any(extinct) & any(!extinct)) %>%
  summarise(p1 = (p.g.e[extinct] * p.ext[1]),
            p2 = (p.g.e[!extinct] * (1 - p.ext[1])),
            p.g.w = p1 / (p1 + p2),
            n.ex = n.g.e[extinct],
            n.sv = n.g.e[!extinct])

ggplot(p.ext.g.overall %>% filter(n.ex > 10 & n.sv > 10)) +
  geom_line(
    aes(
      x = gr, y = p.g.w,
      group = alpha,
      colour = factor(alpha)
    )
  ) +
  facet_wrap(n.pop0 ~ low.var)

###

# I was also interested in P()