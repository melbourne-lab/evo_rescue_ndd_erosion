library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
library(rstan)     #for extract()
library(rstanarm)  #Bayesian multilevel: stan_lmer(), stan_glmer() etc
options(mc.cores = parallel::detectCores())
theme_set(theme_grey()) #rstanarm overrides default ggplot theme: set it back

expmt.sum <- read.csv('experiment/expmt_bottleneck24/out/sim_results.csv')

# Record extinct and create factors (with explicit level ordering)
expmt.dat <- expmt.sum %>%
  group_by(trial) %>%
  mutate(extinct = max(gen) < 15) %>%
  ungroup() %>%
  mutate(isize      = factor(ifelse(n.pop0 > 20,"large","small"), 
                             levels=c('large','small')),
         ddep       = factor(ifelse(alpha == 0,"DI","DD"),
                             levels=c('DI','DD')),
         evolve     = factor(ifelse(evolve,"evolving","not_evolving"), 
                             levels=c("evolving","not_evolving")),
         block      = factor(block),
         rep      = factor(repct),
         sourcepop = factor(sourcepops),
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
  scale_color_manual(values=c('black', 'red')) +
  scale_linetype_manual(values = 1:2) +
  scale_size_manual(values = c(0.75, 0.25)) +
  coord_cartesian(ylim=c(0,200)) +
  facet_grid(rows=vars(isize),cols=vars(ddep)) +
  labs(x="Generation",y="Abundance") +
  theme_bw() +
  theme(legend.position = 'bottom')


# Analysis of extinction - bottlenecked data

# Extract extinction data
expmt.extinct.dat <- expmt.dat %>%
    filter(gen == 1) %>%
    select(trial,isize,ddep,evolve,block,sourcepop,extinct)

# The max likelihood fit is finicky but converges with factors ordered as above.
# The 3 way interaction is not significant.
glm_full <- glmer(extinct ~ isize*ddep*evolve + (1|block) + (1|sourcepop),
              family=binomial, data=expmt.extinct.dat)
summary(glm_full)

# The full additive model
glm_add <- glmer(extinct ~ isize + ddep + evolve + (1|block) + (1|sourcepop), 
              family=binomial, data=expmt.extinct.dat)

# Main effects are all significant (need to confirm via bootstrap)
glm_drop_evolve <- glmer(extinct ~ isize + ddep + (1|block) + (1|sourcepop), 
              family=binomial, data=expmt.extinct.dat)
anova(glm_add,glm_drop_evolve)

glm_drop_ddep <- glmer(extinct ~ isize + evolve + (1|block) + (1|sourcepop), 
              family=binomial, data=expmt.extinct.dat)
anova(glm_add,glm_drop_ddep)

glm_drop_isize <- glmer(extinct ~ ddep + evolve + (1|block) + (1|sourcepop), 
              family=binomial, data=expmt.extinct.dat)
anova(glm_add,glm_drop_isize)

# And we could look at two-way interactions by adding to additive model


# Bayesian fit is good
glm_full_bayes <- stan_glmer(extinct ~ isize*ddep*evolve + (1|block) + (1|sourcepop),
              family=binomial, data=expmt.extinct.dat)
print(summary(glm_full_bayes)[,c("mean","sd","n_eff","Rhat")],digits=3)

# Treatment means and credible intervals
isize_nd <- factor(rep(rep(c("small","large"),each=2),2))
ddep_nd <- factor(rep(c("DD","DI"),4))
evolve_nd <- factor(rep(c("not_evolving","evolving"),each=4))

newd <- data.frame(isize=isize_nd,ddep=ddep_nd,evolve=evolve_nd)
pmu <- posterior_epred(glm_full_bayes, re.form=NA, newdata = newd)
mnmu <- colMeans(pmu)
cred_intervals <- t(apply(pmu,2,quantile,prob=c(0.025,0.975)))
colnames(cred_intervals) <- c("mulo95","muhi95")

mean_prob_df <- cbind(trt=8:1,newd,mnmu,cred_intervals)

# Base plot
plot(mnmu~trt,data=mean_prob_df, ylim=c(0,1),ann=FALSE,axes=FALSE)
segments(mean_prob_df$trt,mean_prob_df$mulo95,mean_prob_df$trt,mean_prob_df$muhi95)
axis(1,at=mean_prob_df$trt,labels=as.character(mean_prob_df$ddep))
axis(1,at=c(2.5,6.5),labels=c("Evolving","Not evolving"),tick=FALSE,line=2)
axis(2,at=seq(0,1,0.25),las=2)
axis(3,at=c(1.5,3.5,5.5,7.5),labels=rep(c("Large","Small"),2),tick=FALSE)
abline(v=c(2.5,6.5),lty=2,col="gray")
abline(v=4.5)
box()
mtext("Extinction probability",2,line=3)

# ggplot
# While facet_grid ordering follows the factor level order, ggplot() doesn't, so
# ddep ordering does not follow the factor level order but is alphabetical. We
# need a fix for that.
# e.g. https://www.r-bloggers.com/2016/11/ordering-categories-within-ggplot2-facets/
mean_prob_df %>%
  ggplot(aes(x=ddep,y=mnmu)) +
  geom_linerange(aes(ymin=mulo95,ymax=muhi95)) +
  geom_point() +
  facet_grid(cols=vars(evolve,isize)) +
  ylim(0,1) +
  labs(x="Density dependence",y="Extinction probability") +
  theme_bw()
