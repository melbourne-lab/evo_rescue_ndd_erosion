# Analyze first time steps in model
# Ran 2000 first-steps for a param combo on Brett's server.

library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

# Parameters used

params = data.frame(end.time = 2,
                    init.row = 1e4,
                    n.loci = 25, 
                    n.pop0 = 20,
                    w.max = 2, 
                    theta = 2.6, 
                    wfitn = sqrt(1 / 0.14 / 2),
                    sig.e = 0)

# The simulation data
firsts = read.csv('analyze_model/genotypic changes/first_timers.csv')

# Look at change in allele frequencies 
first.long = firsts %>%
  select(-c(z_i, fem)) %>%
  gather(locus, gval, -c(i, g_i, w_i, r_i, trial, gen)) %>%
  mutate(lcopy = gsub('\\d', '', locus),
         locus = gsub('[ab]', '', locus))

# Allele frequencies - wide form, both generations
p.all.wide = first.long %>%
  group_by(locus, trial, gen) %>%
  summarise(mean.val = mean(gval > 0)) %>%
  mutate(gen = paste0('gen', gen)) %>%
  spread(gen, mean.val)
  
p.all.wide %>%
  ggplot() +
  geom_line(aes(x = gen1, y = gen2, group = trial), size = 0.02)

p.all.wide %>%
  ggplot() +
  geom_point(aes(x = gen1, y = gen2), 
             position = position_jitter(height = 0.01, width = 0.01),
             alpha = 0.01)

p.all.wide %>% lm(formula = gen2 ~ gen1) %>% summary()
# 0.015629  pm 2 *  0.003291 is the change over time

with(params %>% 
       mutate(p = 0.5, q = 1 - p,
              wA = exp(-(theta - 1/sqrt(n.loci))^2 / (2*wfitn^2)),
              w0 = exp(-(theta)^2 / (2*wfitn^2)),
              wa = exp(-(theta + 1/sqrt(n.loci))^2 / (2*wfitn^2))),
     (p * (wA * p*(p+1) + w0 * q*(2*p+1) + wa * q^2)) / (2*(wA * p^2 + w0*2*p*q + wa * q^2)) - p)
# change in p should be # 0.01811, which is close

dt(0.01811731 - 0.015629, df = 4998)
# approx 0.4 so based on our data fail to reject null hypothesis
# that the observed value is the true value above

p.all.wide %>%
  filter(trial %in% 1:3) %>%
  gather(gen, pbar, -c(locus, trial)) %>%
  ggplot() +
  geom_line(aes(x = gen, y = pbar, group = interaction(trial, locus), colour = factor(trial)))

### What could influence the rate of change?

# Think mechanistically.
# How many offspring does an AA individual have?
# Aa? aa?

first = firsts %>% filter(trial %in% 1)

head(first)

fy = first %>% 
  select(i, g_i, z_i, w_i, r_i, fem, gen, a1, b1) %>%
  mutate(l1t1pos = (a1 > 0) + (b1 > 0),
         l1t2pos = r_i * l1t1pos,
         l1t1w = w_i * l1t1pos)

cov(fy$l1t1pos[fy$fem], fy$g_i[fy$fem])

sum(fy$l1t1pos) / (2*nrow(fy))
sum(fy$l1t2pos) / (2*sum(fy$r_i))


das.covars = firsts %>% 
  select(i, g_i, z_i, w_i, r_i, fem, gen, a1, b1, trial) %>%
  group_by(trial, gen) %>%
  summarise(covars = cov((a1 > 0) + (b1 > 0), z_i))

das.covars %>% filter(gen %in% 1) %>% ggplot() + geom_histogram(aes(x = covars), binwidth = 0.05)
# So clearly there's some positive covariance here.
# Not sure why it's so g-darn hard to find it though...

### Looking at the Lande expression.

gen.vars = first.long %>%
  group_by(trial, gen, locus) %>%
  summarise(p = mean(gval > 0),
            z = mean(g_i)) %>%
  mutate(local.var = 2 * p * (1 - p) / params$n.loci) %>%
  group_by(trial, gen) %>%
  summarise(p = mean(p),
            v = sum(local.var),
            z = mean(z))

is.this.it = gen.vars %>%
  group_by(trial) %>%
  summarise(k = params$wfitn^2 / (params$wfitn^2 + v[1]),
            g = params$theta - z[1],
            p1 = p[1],
            eg = g * k,
            hg = params$theta - z[2],
            z2 = z[2])


head(is.this.it, 20)
plot(hg ~ eg, is.this.it, asp = 1); abline(0,1)
plot(z2 ~ p1, is.this.it)

is.this.it %>%
  ggplot() +
  geom_point(aes(x = p1, hg))
is.this.it %>%
  ggplot() +
  geom_point(aes(x = p1, z2))

eg = gen.vars %>%
  mutate(k = params$wfitn^2 / (params$wfitn^2 + v),
         g = params$theta - z) %>%
  group_by(trial) %>%
  summarise(delta.g = g[1] * k[1],
            edelt.g = g[2],
            v = v[1]) 

eg %>%
  ggplot() +
  geom_point(aes(x = edelt.g, y = delta.g, colour = v), 
             alpha = 0.5) +
  geom_segment(aes(x = 1, xend = 3.75, y = 1, yend = 3.75),
               colour = 'blue') +
  scale_color_viridis_c()

eg %>%
  ggplot() +
  geom_point(aes(x = v, y = delta.g), 
             alpha = 0.5)

gen.vars %>%
  select(trial, gen, p) %>%
  mutate(gen = paste0('gen', gen)) %>%
  spread(gen, p) %>%
  ggplot() +
  geom_point(aes(x = gen1, y = gen2))

gen.vars %>%
  select(trial, gen, z) %>%
  mutate(z = params$theta - z,
         gen = paste0('gen', gen)) %>%
  spread(gen, z) %>%
  merge(y = gen.vars %>% filter(gen < 2) %>% select(trial, v)) %>%
  filter(v > 0.475) %>%
  ggplot() +
  geom_point(aes(x = gen1, gen2, fill = v),
             shape = 21, colour = 'gray55',
             alpha = 0.85) +
  geom_segment(aes(x = 2, xend = 3, y = 2, yend = 3),
               colour = 'darkblue') +
  scale_fill_viridis_c()

# Suggests... no relationship between overall genetic diversity and change.
# Strange.

# Go allele by allele?

hmm = first.long %>%
  group_by(locus, trial, gen) %>%
  summarise(p = mean(gval > 0)) %>%
  group_by(trial, gen) %>%
  mutate(p.all = mean(p)) %>%
  ungroup()

hmm %>%
  filter(trial < 500) %>%
  group_by(trial) %>%
  mutate(pall1 = p.all[1]) %>%
  ungroup() %>%
  mutate(t.jitter = gen + rnorm(nrow(.), 0, 0.0125)) %>%
  ggplot() +
  geom_line(aes(x = t.jitter, y = p, colour = pall1,
                group = interaction(trial, locus)),
            size = 0.0125) +
  geom_point(aes(x = t.jitter, y = p, colour = p.all),
             alpha = 0.5)

# Man what the hell??

try.this = hmm %>%
  select(-p.all) %>%
  group_by(trial, locus) %>%
  summarise(p1.var = 2 * p[1] * (1 - p[1]) / params$n.loci,
            del.p = p[2] - p[1]) %>%
  group_by(trial) %>%
  mutate(all.var = sum(p1.var),
         all.dp = mean(del.p))

try.this %>%
  ggplot(aes(x = p1.var, y = del.p)) +
  geom_point(position = position_jitter(width = 0.0001),
             alpha = 0.5) +
  stat_smooth(method = 'lm')
# what the hell...?

cor(try.this$del.p, try.this$p1.var)
# this is very, very small.

mean(try.this$del.p)
mean(try.this$del.p > 0)
# So only 53% of alleles increase???

p.overall = first.long %>%
  group_by(trial, gen) %>%
  summarise(p = mean(gval > 0),
            n = length(unique(i))) 

p.overall %>%
  ggplot(aes(x = gen, y = p)) +
  geom_line(aes(group = trial), size = 0.2) +
  stat_smooth(method = 'lm')


