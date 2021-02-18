library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
library(rstan)     #for extract()
library(rstanarm)  #Bayesian multilevel: stan_lmer(), stan_glmer() etc
options(mc.cores = parallel::detectCores())
theme_set(theme_grey()) #rstanarm overrides default ggplot theme: set it back

expmt.sum = read.csv('experiment/expmt_bottleneck24/out/sim_results.csv')

# Record extinct and create factors
expmt.dat <- expmt.sum %>%
  group_by(trial) %>%
  mutate(extinct = max(gen) < 15) %>%
  ungroup() %>%
  mutate(isize      = factor(n.pop0 > 20, labels = c('small', 'large')),
         ddep       = factor(alpha  > 0,  labels = c('dens indep', 'dens dep')),
         block      = factor(block),
         repct      = factor(repct),
         sourcepops = factor(sourcepops),
         trial      = factor(trial)
         )

# Plot with treatments in panels
expmt.dat %>%    
  ggplot(aes(x = gen, y = n)) +
  geom_line(
    aes(
      group = trial,
      linetype = evolve,
      size = evolve,
      colour = extinct
    )
  ) +
  scale_color_manual(values = c('black', 'red')) +
  scale_linetype_manual(values = 2:1) +
  scale_size_manual(values = c(0.25, 0.75)) +
  #scale_y_log10() +
  ylim(0,200) +
  facet_wrap(~ paste(isize, ddep, sep = ', ')) +
  theme(legend.position = 'bottom')


# Analysis of extinction - bottlenecked data

# Extract extinction data
expmt.extinct.dat <- expmt.dat %>%
    filter(gen ==1 ) %>%
    select(trial,isize,ddep,evolve,block,sourcepops,extinct)


glm1 <- glmer(extinct ~ isize + (1|block) + (1|sourcepops), 
              family=binomial, data=expmt.extinct.dat)
summary(glm1)

glm2 <- glmer(extinct ~ ddep + (1|block) + (1|sourcepops),
              family=binomial, data=expmt.extinct.dat)
summary(glm2)

glm3 <- glmer(extinct ~ evolve + (1|block) + (1|sourcepops),
              family=binomial, data=expmt.extinct.dat)
summary(glm3)

# The max likelihood fit fails to converge (as expected)
glm6 <- glmer(extinct ~ isize*ddep*evolve + (1|block) + (1|sourcepops),
              family=binomial, data=expmt.extinct.dat)

# Bayesian fit is good
glm6 <- stan_glmer(extinct ~ isize*ddep*evolve + (1|block) + (1|sourcepops),
              family=binomial, data=expmt.extinct.dat)
print(summary(glm6)[,c("mean","sd","n_eff","Rhat")],digits=3)

