# Estimating extinction probabilities, conditioned on state vars, using Bayes's rule

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
  mutate(vr = round(v * 20) / 20) %>%
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

p.ext.v %>% 
  ungroup() %>%
  filter(gen %in% c(1, 5), 
         n.ex > 10 & n.sv > 10) %>%
  mutate(n.pop0 = factor(ifelse(n.pop0 > 50, 'large', 'small')),
         low.var = factor(ifelse(low.var, 'low diversity', 'high diversity'))) %>%
  ggplot(aes(x = vr, y = p.ex.vr)) +
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
  scale_x_reverse() +
  scale_shape_manual(values = c(16, 17, 15),
                     'Generation') +
  scale_color_manual(values = c('black', 'purple'),
                     'Alpha') +
  facet_wrap(~ paste(n.pop0, low.var, sep = ', '), nrow = 2) +
  labs(x = 'Additive genetic variance', y = 'Prob. extinct') +
  theme_bw() +
  theme(legend.position = 'bottom') #+
  ggsave('~/Dropbox/rescue_ndd_paper_2020/figures/p_extinct_sig_genvar.pdf',
         width = 6, height = 6)

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

p.ext.w %>% 
  ungroup() %>%
  filter(gen %in% (1 + (0:2)*3), 
         n.ex > 10 & n.sv > 10) %>%
  mutate(n.pop0 = factor(ifelse(n.pop0 > 50, 'large', 'small')),
         low.var = factor(ifelse(low.var, 'low diversity', 'high diversity'))) %>%
  ggplot(aes(x = wr, y = p.e.w)) +
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
  scale_shape_manual(values = c(16, 17, 15),
                     'Generation') +
  scale_color_manual(values = c('black', 'purple'),
                     'Alpha') +
  facet_wrap(~ paste(n.pop0, low.var, sep = ', '), nrow = 2) +
  labs(x = 'Intrinsic fitness (W)', y = 'Prob. extinct') +
  theme_bw() +
  theme(legend.position = 'bottom') #+
  ggsave('~/Dropbox/rescue_ndd_paper_2020/figures/p_extinct_w_fitn.pdf',
         width = 6, height = 6)

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

p.ext.g %>% 
  ungroup() %>%
  filter(gen %in% (1 + (0:2)*3), 
         n.ex > 10 & n.sv > 10) %>%
  mutate(n.pop0 = factor(ifelse(n.pop0 > 50, 'large', 'small')),
         low.var = factor(ifelse(low.var, 'low diversity', 'high diversity'))) %>%
  ggplot(aes(x = gr, y = p.g.w)) +
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
  scale_shape_manual(values = c(16, 17, 15),
                     'Generation') +
  scale_color_manual(values = c('black', 'purple'),
                     'Alpha') +
  facet_wrap(~ paste(n.pop0, low.var, sep = ', '), nrow = 2) +
  labs(x = 'Genotype (g)', y = 'Prob. extinct') +
  theme_bw() +
  theme(legend.position = 'bottom') +
  ggsave('~/Dropbox/rescue_ndd_paper_2020/figures/p_extinct_g_geno.pdf',
         width = 6, height = 6)

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

# ----

# Fixations?

p.ext.fp = all.data %>%
  mutate(pfp = 25 * p.fix.pos) %>%
  group_by(n.pop0, low.var, alpha, gen, extinct) %>%
  mutate(n.ex = n()) %>%
  group_by(n.pop0, low.var, alpha, gen, extinct, pfp) %>%
  summarise(n.w.fp = n(),
            n.e = n.ex[1],
            p.w.fp = n.w.fp / n.e) %>%
  merge(y = pext) %>%  
  group_by(n.pop0, low.var, alpha, gen, pfp) %>%
  filter(any(extinct) & any(!extinct)) %>%
  summarise(p1 = (p.w.fp[extinct] * p.ext[1]),
            p2 = (p.w.fp[!extinct] * (1 - p.ext[1])),
            p.e.fp = p1 / (p1 + p2),
            n.ex = n.w.fp[extinct],
            n.sv = n.w.fp[!extinct])

p.ext.fp %>% 
  ungroup() %>%
  filter(gen %in% (1 + (0:2)*3), 
         n.ex > 10 & n.sv > 10) %>%
  mutate(n.pop0 = factor(ifelse(n.pop0 > 50, 'large', 'small')),
         low.var = factor(ifelse(low.var, 'low diversity', 'high diversity'))) %>%
  ggplot(aes(x = pfp, y = p.e.fp)) +
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
  scale_shape_manual(values = c(16, 17, 15),
                     'Generation') +
  scale_color_manual(values = c('black', 'purple'),
                     'Alpha') +
  facet_wrap(~ paste(n.pop0, low.var, sep = ', '), nrow = 2) +
  labs(x = 'Loci at fixation (positive)', y = 'Prob. extinct') +
  theme_bw() +
  theme(legend.position = 'bottom') +
  ggsave('~/Dropbox/rescue_ndd_paper_2020/figures/p_extinct_pf_posfix.pdf',
         width = 6, height = 6)
  
p.ext.fn = all.data %>%
  mutate(pfn = 25 * p.fix.neg) %>%
  group_by(n.pop0, low.var, alpha, gen, extinct) %>%
  mutate(n.ex = n()) %>%
  group_by(n.pop0, low.var, alpha, gen, extinct, pfn) %>%
  summarise(n.w.fn = n(),
            n.e = n.ex[1],
            p.w.fn = n.w.fn / n.e) %>%
  merge(y = pext) %>%  
  group_by(n.pop0, low.var, alpha, gen, pfn) %>%
  filter(any(extinct) & any(!extinct)) %>%
  summarise(p1 = (p.w.fn[extinct] * p.ext[1]),
            p2 = (p.w.fn[!extinct] * (1 - p.ext[1])),
            p.e.fn = p1 / (p1 + p2),
            n.ex = n.w.fn[extinct],
            n.sv = n.w.fn[!extinct])

p.ext.fn %>% 
  ungroup() %>%
  filter(gen %in% (1 + (0:2)*3), 
         n.ex > 10 & n.sv > 10) %>%
  mutate(n.pop0 = factor(ifelse(n.pop0 > 50, 'large', 'small')),
         low.var = factor(ifelse(low.var, 'low diversity', 'high diversity'))) %>%
  ggplot(aes(x = pfn, y = p.e.fn)) +
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
  scale_shape_manual(values = c(16, 17, 15),
                     'Generation') +
  scale_color_manual(values = c('black', 'purple'),
                     'Alpha') +
  facet_wrap(~ paste(n.pop0, low.var, sep = ', '), nrow = 2) +
  labs(x = 'Loci at fixation (negatve)', y = 'Prob. extinct') +
  theme_bw() +
  theme(legend.position = 'bottom') +
  ggsave('~/Dropbox/rescue_ndd_paper_2020/figures/p_extinct_pn_negfix.pdf',
       width = 6, height = 6)

###

p.ext.n = all.data %>%
  mutate(nr = round(n / 5) * 5) %>%
  group_by(n.pop0, low.var, alpha, gen, extinct) %>%
  mutate(n.ex = n()) %>%
  group_by(n.pop0, low.var, alpha, gen, extinct, nr) %>%
  summarise(n.n.ex = n(),
            n.n = n.ex[1],
            p.n.ex = n.n.ex / n.n) %>%
  merge(y = pext) %>%
  group_by(n.pop0, low.var, alpha, gen, nr) %>%
  filter(any(extinct) & any(!extinct)) %>%
  summarise(p1 = (p.n.ex[extinct] * p.ext[1]),
            p2 = (p.n.ex[!extinct] * (1 - p.ext[1])),
            p.n.n = p1 / (p1 + p2),
            n.ex = n.n.ex[extinct],
            n.sv = n.n.ex[!extinct])

p.ext.n %>% 
  ungroup() %>%
  filter(gen %in% (c(2, 5, 8)), 
         n.ex > 5 & n.sv > 5) %>%
  mutate(n.pop0 = factor(ifelse(n.pop0 > 50, 'large', 'small')),
         low.var = factor(ifelse(low.var, 'low diversity', 'high diversity'))) %>%
  ggplot(aes(x = nr, y = p.n.n)) +
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
  scale_shape_manual(values = c(16, 17, 15),
                     'Generation') +
  scale_color_manual(values = c('black', 'purple'),
                     'Alpha') +
  facet_wrap(~ paste(n.pop0, low.var, sep = ', '), nrow = 2) +
  labs(x = 'Loci at fixation (negatve)', y = 'Prob. extinct') +
  theme_bw() +
  theme(legend.position = 'bottom')

### Another attempt incorporating population size

# Population size is a latent variable in the above.
# But, we can estimate the joint probabilities of extinction, conditioned on population size
# using Bayes's Rule in a similar way to the above.

# A problem is the population size. As with above, I'll bin the population size variable
# (the binning here may change the output)
# In particular, bin on a log scale or on the natural scale? I try this here with the log scale
# (log 2) but the code would work with a natural scaling too.

# Aux wrapper function (tested and proven) for binning
# x is the vector of numbers to bin, nums are the bin boundaries
numbin = function(x, nums) {
  expand.grid(input = x, ns = rev(nums)) %>%
    mutate(i = rep(1:length(x), times = length(nums))) %>%
    filter(input >= ns) %>%
    distinct(i, input, .keep_all = TRUE) %>%
    arrange(i) %>%
    select(ns) %>%
    unlist()
}

# The numeric scales used
# numscale = c(1, 2, 3, 5, 10, 20, 30, 50, 100, 200, 300, 500)
# numscale = c(1, 2, 5, 10, 20, 50, 100, 200, 500)
numscale = c(1, 5, 10, 20, 40, 80, 160) # nice log 2 scale

# Generate the distributions
p.ext.w = all.data %>%
  mutate(wr = round(wbar * 10) / 10,
         nr = numbin(n, numscale)) %>%
  group_by(n.pop0, low.var, alpha, gen, extinct) %>%
  mutate(n.ex = n()) %>%
  group_by(n.pop0, low.var, alpha, gen, extinct, wr, nr) %>%
  summarise(n.w.e = n(),
            n.e = n.ex[1],
            p.w.e = n.w.e / n.e) %>%
  merge(y = pext) %>%  
  group_by(n.pop0, low.var, alpha, gen, wr, nr) %>%
  filter(any(extinct) & any(!extinct)) %>%
  summarise(p1 = (p.w.e[extinct] * p.ext[1]),
            p2 = (p.w.e[!extinct] * (1 - p.ext[1])),
            p.e.w = p1 / (p1 + p2),
            n.ex = n.w.e[extinct],
            n.sv = n.w.e[!extinct])

# This looks like the best plot thus far
p.ext.w %>% 
  ungroup() %>%
  filter(gen %in% c(2, 6), nr > 2, nr < 100, # note - excludes populations above size 160
         n.ex > 10 & n.sv > 10) %>%
  mutate(n.pop0 = factor(ifelse(n.pop0 > 50, 'large', 'small')),
         low.var = factor(ifelse(low.var, 'low diversity', 'high diversity'))) %>%
  ggplot(aes(x = wr, y = p.e.w)) +
  geom_line(
    aes(
      group = interaction(gen, alpha, nr),
      linetype = factor(alpha),
      colour = factor(nr)
    ),
    size = 1
  ) +
  geom_point(
    aes(
      colour = factor(nr),
      shape = factor(alpha)
    ),
    size = 3
  ) +
  scale_shape_manual(values = c(16, 17), #, 15),
                     'Generation') +
  scale_color_brewer(palette = 'Spectral',
                     'Alpha') +
  facet_wrap(gen ~ paste(n.pop0, low.var, sep = ', ') , nrow = 2) +
  labs(x = 'Intrinsic fitness (W)', y = 'Prob. extinct') +
  theme_bw() +
  theme(legend.position = 'bottom') 

# may want to try the above plot using fill instead of colour,
# and shapes with a black packground (to make things pop)

# A less likeable plot, faceting on size and featuring several generations
p.ext.w %>% 
  ungroup() %>%
  filter(gen %in% c(2, 6),
         n.ex > 10 & n.sv > 10) %>%
  mutate(n.pop0 = factor(ifelse(n.pop0 > 50, 'large', 'small')),
         low.var = factor(ifelse(low.var, 'low diversity', 'high diversity'))) %>%
  ggplot(aes(x = wr, y = p.e.w)) +
  geom_line(
    aes(
      group = interaction(gen, alpha, nr),
      linetype = factor(alpha),
      colour = factor(nr)
    ),
    size = 0.75
  ) +
  geom_point(
    aes(
      colour = factor(nr),
    ),
    size = 2
  ) +
  scale_shape_manual(values = c(16, 17, 15),
                     'Generation') +
  scale_color_brewer(palette = 'YlOrRd',
                     'Alpha') +
  # facet_wrap(paste(n.pop0, low.var, sep = ', ') ~ nr) +
  facet_grid(rows = vars(n.pop0, low.var), cols = vars(nr)) +
  labs(x = 'Intrinsic fitness (W)', y = 'Prob. extinct') +
  theme_bw() +
  theme(legend.position = 'bottom') 
