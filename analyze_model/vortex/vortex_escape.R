# What allows populations to escape the extinction vortex?
# Here, looking at all of the populations below a size threshold
# and looking for what varies by extinction/survival.
# SN 29 Jun 2021 - quick and dirty script

library(ggplot2)
library(dplyr)
library(tidyr)

allpops = read.csv('simulations/outputs/final_results/alldata_combined.csv')

head(allpops)

size.thresh = 11

gdiffs = allpops %>%
  group_by(n.pop0, low.var, alpha, trial) %>%
  filter(any(n < size.thresh))  %>%
  mutate(gdiff = c(diff(gbar), NA)) %>% 
  filter(n < size.thresh) %>%
  filter(!is.na(gdiff)) %>%
  distinct(n.pop0, low.var, alpha, trial, .keep_all = TRUE) %>%
  group_by(n.pop0, low.var, alpha, extinct) %>%
  summarise(mean.gdiff = mean(gdiff),
            sdev.gdiff = sd(gdiff),
            n = n())

gdiffs %>%
  ggplot(aes(x = paste(n.pop0, low.var))) +
  geom_point(aes(y = mean.gdiff, group = interaction(extinct, alpha), 
                 shape = factor(alpha), colour = extinct),
             size = 2,
             position = position_dodge(width = 0.25)) +
  geom_linerange(aes( 
                   ymin = mean.gdiff - 2 * sqrt(sdev.gdiff / n),
                   ymax = mean.gdiff + 2 * sqrt(sdev.gdiff / n),
                   group = interaction(extinct, alpha),
                   colour = extinct),
               position = position_dodge(width = 0.25)) +
  scale_color_manual(values = c('black', 'red'))

wdiffs = allpops %>%
  group_by(n.pop0, low.var, alpha, trial) %>%
  filter(any(n < size.thresh))  %>%
  mutate(wdiff = c(diff(wbar), NA)) %>% 
  filter(n < size.thresh) %>%
  filter(!is.na(wdiff)) %>%
  distinct(n.pop0, low.var, alpha, trial, .keep_all = TRUE) %>%
  group_by(n.pop0, low.var, alpha, extinct) %>%
  summarise(mean.wdiff = mean(wdiff),
            sdev.wdiff = sd(wdiff),
            n = n())

wdiffs %>%
  ggplot(aes(x = paste(n.pop0, low.var))) +
  geom_point(aes(y = mean.wdiff, group = interaction(extinct, alpha), 
                 shape = factor(alpha), colour = extinct),
             size = 2,
             position = position_dodge(width = 0.25)) +
  geom_linerange(aes( 
    ymin = mean.wdiff - 2 * sqrt(sdev.wdiff / n),
    ymax = mean.wdiff + 2 * sqrt(sdev.wdiff / n),
    group = interaction(extinct, alpha),
    colour = extinct),
    position = position_dodge(width = 0.25)) +
  scale_color_manual(values = c('black', 'red'))


posfixdiffs = allpops %>%
  group_by(n.pop0, low.var, alpha, trial) %>%
  filter(any(n < size.thresh))  %>%
  mutate(pfdiff = c(diff(p.fix.pos), NA)) %>% 
  filter(n < size.thresh) %>%
  filter(!is.na(pfdiff)) %>%
  distinct(n.pop0, low.var, alpha, trial, .keep_all = TRUE) %>%
  group_by(n.pop0, low.var, alpha, extinct) %>%
  summarise(mean.pfdiff = mean(pfdiff),
            sdev.pfdiff = sd(pfdiff),
            n = n())

posfixdiffs %>%
  ggplot(aes(x = paste(n.pop0, low.var))) +
  geom_point(aes(y = mean.pfdiff, group = interaction(extinct, alpha), 
                 shape = factor(alpha), colour = extinct),
             size = 2,
             position = position_dodge(width = 0.25)) +
  geom_linerange(aes( 
    ymin = mean.pfdiff - 2 * sqrt(sdev.pfdiff / n),
    ymax = mean.pfdiff + 2 * sqrt(sdev.pfdiff / n),
    group = interaction(extinct, alpha),
    colour = extinct),
    position = position_dodge(width = 0.25)) +
  scale_color_manual(values = c('black', 'red'))
 
# very much not what I expected...
# expected black to be higher than red i.e. more fixations in surviving populations...

# next to try.., just remove rows preceding size below size thresh
# and plot means...

#%>%
  ggplot(aes(x = gen, y = gbar)) +
  geom_line(aes(group = trial, colour = extinct), size = 0.1) +
  geom_line(aes(group = trial, colour = extinct, alpha = !extinct), size = 0.1) +
  facet_wrap(n.pop0 + low.var ~ alpha, nrow = 2) +
  scale_color_manual(values = c('black', 'red'))

allpops %>%
  group_by(n.pop0, low.var, alpha, trial) %>%
  filter(trial < 100) %>%
  ggplot(aes(x = gen, y = n)) +
  geom_line(aes(group = trial, color = extinct)) +
  facet_wrap(n.pop0 + low.var ~ alpha, nrow = 2) +
  scale_color_manual(values = c('black', 'red')) +
  scale_y_log10()
  
