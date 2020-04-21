# Script to look at change in allele frequencies

library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

#### First, try an expression for change in allele frequencies over time.

## Expression depends on p, d, w, m

# Set w
w = sqrt(1 / 0.14 / 2)

# Make a grid of values
# and calculate change in allele frequency
a = expand.grid(d = (0:260) / 100,
                p = (1:99) / 100,
                m = c(20, 40, 80, 160)) %>%
  mutate(q = 1 - p,
         w_A = exp(-(d - 1/sqrt(m))^2 / (2*w^2)),
         w_0 = exp(-d^2 / (2*w^2)),
         w_a = exp(-(d + 1/sqrt(m))^2 / (2*w^2)),
         N_p = p * (w_A*p*(p+1) + w_0*q*(2*p+1) + w_a * q^2),
         N   = 2 * (p^2 * w_A + 2*p*q*w_0 + q^2 * w_a),
         p_2 = N_p / N,
         del = p_2 - p)

# Plot results
ggplot(a) +
  geom_tile(aes(x = p, y = d, fill = del)) +
  facet_wrap(~ m) +
  theme_classic()

# Looks right - most change when m is small, d is large, p = 0.5

### Now try to iterate out changes in p and d

# Initialize
# m = number of alleles
# wsq = fitness term used in model fitting
# ts = timestep
# pt = allele frequency over time
#   init at 1/2
# dt = distance from phenotypic optimum over time
#   init at 2.6 (used in sims)
m = 25
ts = 1:15
pt = 1:15
dt = 1:15
pt[1] = 0.5
dt[1] = 2.6

# Iterate forward pt, dt
for (tt in ts[-1]) {
  p = pt[tt-1]
  q = 1 - p
  d = dt[tt-1]
  wA = exp(-(d - 1/sqrt(m))^2 / (2*w^2))
  w0 = exp(-d^2 / (2*w^2))
  wa = exp(-(d + 1/sqrt(m))^2 / (2*w^2))
  Np = p * (wA*p*(p+1) + w0*q*(2*p+1) + wa * q^2) # this is correct... just numerical error
  N  = 2 * (p^2 * wA + 2*p*q*w0 + q^2 * wa)
  pt[tt] = Np / N
  dt[tt] = dt[1] - sqrt(m) * (pt[tt] - (1 - pt[tt]))
}

# Plot results.
plot(ts, pt, type = 'l')
plot(ts, dt, type = 'l')

# Any signal of logginess in d?
plot(ts, dt, type = 'l', log = 'y')
# Sort of linear...
dt %>% log() %>% diff() %>% exp() %>% plot(type = 'l')
# If it was geometric, this curve would be flat.
# So, sub-geometric?

#### Now compare with simulation results

# Read in simulation results
egk = read.csv('analyze_model/eg_1k_n0_20.csv')

head(egk)

# Plot genotypic change against time
#   black = observed in data
#   purple = predicted above
ggplot(egk, aes(x = gen)) +
  geom_line(aes(y = 2.6 - g.bar)) +
  geom_ribbon(aes(ymin = 2.6 - (g.bar + 2 * sqrt(g.var / n.trials)),
                  ymax = 2.6 - (g.bar - 2 * sqrt(g.var / n.trials))),
              alpha = 0.2) +
  geom_line(aes(x = ts, y = dt),
            data = data.frame(ts = ts, dt = dt),
            colour = 'purple') +
  geom_point(aes(y = 2.6 - g.bar, colour = (n.trials / 1000))) +
  scale_color_viridis_c(option = 'B')

ggplot(egk %>% filter(gen < 4), aes(x = gen)) +
  geom_line(aes(y = 2.6 - g.bar)) +
  geom_ribbon(aes(ymin = 2.6 - (g.bar + 2 * sqrt(g.var / n.trials)),
                  ymax = 2.6 - (g.bar - 2 * sqrt(g.var / n.trials))),
              alpha = 0.2) +
  geom_line(aes(x = ts, y = dt),
            data = data.frame(ts = ts[1:3], dt = dt[1:3]),
            colour = 'purple') +
  geom_point(aes(y = 2.6 - g.bar, colour = (n.trials / 1000))) +
  scale_color_viridis_c(option = 'B')

# Genotypic change as modeled above is too fast.


