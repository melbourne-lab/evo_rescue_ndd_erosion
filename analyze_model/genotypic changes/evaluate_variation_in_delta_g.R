# A script looking at variation in phenotypic change and its effects on change
# in fitness.

# Wrapper for calculating fitness (relative to max)
wex = function(x, w) exp(-x^2 / (2*w^2))

# some constants
# selection strength
w = sqrt(1/0.14/2)
# mean phenotypic change
Z = -0.1
# number of draws from normal (genotypc changes) to simulate
n.draws = 1000

# Grid with
# - coefficient of variation in genotypic change (sd(delta z) / mean(delta z))
# - z0 (keeping in mind that effect may depend on intial phenotype)
zgrid = expand.grid(cvz = (0:50)/10, z = (0:30)/10) %>% mutate(edw = NA)

set.seed(636)

for (k in 1:nrow(zgrid)) {
  dws = 1:ntrials
  # set coeffcient of varation
  cvz = zgrid$cvz[k]
  # set initial genotype
  z1  = zgrid$z[k]
  # draw updated phenotype, z1 + [~ N(Z, cvz * Z)]
  z2 = z1 + rnorm(n.draws, Z, -cvz * Z)
  # calculate delta Ws (changes in fitness)
  dws = wex(z2, w) - wex(z1, w)
  # Take the mean of this all of the delta Ws
  zgrid$edw[k] = mean(dws)
  print(k)
}

ggplot(zgrid) +
  geom_tile(aes(x = z, y = cvz, fill = edw)) +
  scale_fill_gradient2(low = 'red', mid = 'white', midpoint = 0, high = 'blue')
            
ggplot(zgrid) +
  geom_raster(aes(x = z, y = cvz, fill = edw)) +
  scale_fill_gradient2(low = 'red', mid = 'white', midpoint = 0, high = 'royalblue') +
  theme_bw()

# Normalize the expected change by the chnge expected in the deterministic setting:
zgrid = zgrid %>% mutate(edw.norm = edw / (wex(z + Z, w) - wex(z, w)))

ggplot(zgrid %>% filter(z > 0.25)) +
  geom_raster(aes(x = z, y = cvz, fill = edw.norm)) +
  scale_fill_gradient2(low = 'red', mid = 'white', midpoint = 0, high = 'royalblue') +
  theme_bw()

# Hmm okay this is the same.

# Now, mean change in geno is a function of distance

kgrid = expand.grid(cvz = (0:50)/10, z = (0:30)/10) %>% mutate(edw = NA)

set.seed(636)

k = 0.95

for (i in 1:nrow(zgrid)) {
  dws = 1:ntrials
  cvz = kgrid$cvz[i]
  z1  = kgrid$z[i]
  # draw a random value of dws ~ N(Z, cvz * Z)
  z2 = rnorm(n.draws, k*z1, cvz * (1-k)*z1)
  # calculate dw
  dws = wex(z2, w) - wex(z1, w)
  kgrid$edw[i] = mean(dws)
  print(i)
}

ggplot(kgrid) +
  geom_raster(aes(x = z, y = cvz, fill = edw)) +
  scale_fill_gradient2(low = 'red', mid = 'white', midpoint = 0, high = 'royalblue') +
  theme_bw()
