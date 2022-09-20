library(ggplot2)
library(dplyr)

rbo = rainbow(10)

zdf = data.frame(z = (-150:400)/100) 

z0 = zdf %>% mutate(pz0 = dnorm(z, 0, 0.4))

z1 = zdf %>% 
  mutate(pz0 = dnorm(z, 0, 0.4)) %>%
  mutate(pz1 = dnorm(z, 2.75, 0.4))

z1 %>%
  ggplot(aes(x = z)) +
  geom_line(aes(y = pz0), colour = 'red') +
  geom_line(aes(y = pz1), colour = 'white', linetype = 2) +
  labs(x = '', y = '') +
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

ggsave(filename = '~/Documents/Research/boulder/fifth_semester/presentation/rescue_norms1.pdf',
       height = 8, width = 8)

z1 %>%
  ggplot(aes(x = z)) +
  geom_line(aes(y = pz0), colour = 'red') +
  geom_line(aes(y = pz1), colour = 'blue', linetype = 2) +
  labs(x = '', y = '') +
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

ggsave(filename = '~/Documents/Research/boulder/fifth_semester/presentation/rescue_norms2.pdf',
       height = 8, width = 8)


expand.grid(z = zdf$z, t = 0:12) %>%
  mutate(zbar = 2.75 * (1 - 0.95^t)) %>%
  mutate(pzt = dnorm(z, zbar, 0.4)) %>%
  #mutate(col = rbo[t]) %>%
  ggplot(aes(x = z)) +
  geom_line(aes(y = pzt, group = t, colour = t)) +
  scale_color_gradient(low = 'red', high = 'blue') +
  # geom_segment(aes(x = 0, xend = 0, y = 0, yend = 1), colour = 'red') +
  # geom_segment(aes(x = 2.75, xend = 2.75, y = 0, yend = 1), linetype = 2, colour = 'blue') +
  geom_line(data = z1, aes(x = z, y = pz1), linetype = 2, colour = 'blue') +
  labs(x = '', y = '') +
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

ggsave(filename = '~/Documents/Research/boulder/fifth_semester/presentation/rescue_norms3.pdf',
       height = 8, width = 8)

####

tt = 1:100

k1 = 0.05
s1 = 0.15

theta1 = k1 * tt + rnorm(length(tt), 0, s1)
plot(tt, theta1, type = 'l', ylim = c(-1, 5),
     xlab = 't', ylab = expression(theta[t]))
pdf('~/Documents/Research/boulder/fifth_semester/presentation/theta_fast.pdf')
dev.off()

k2 = 0.025
s2 = 0.5

theta2 = k2 * tt + rnorm(length(tt), 0, s2)
plot(tt, theta2, type = 'l', ylim = c(-1, 5),
     xlab = 't', ylab = expression(theta[t]))
pdf('~/Documents/Research/boulder/fifth_semester/presentation/theta_noisy.pdf')
dev.off()

ggplot(data = data.frame(t = tt, theta = theta1),
       aes(x = t, y = theta)) +
  geom_line(colour = 'red') +
  scale_x_continuous(breaks = NULL) +
  labs(x = 't', y = expression(theta[t])) +
  lims(y = c(-1, 5.5)) +
  theme(panel.background = element_blank())

ggsave('~/Documents/Research/boulder/fifth_semester/presentation/theta_fast.pdf',
       height = 3, width = 3)

ggplot(data = data.frame(t = tt, theta = theta2),
       aes(x = t, y = theta)) +
  geom_line(colour = 'blue') +
  scale_x_continuous(breaks = NULL) +
  labs(x = 't', y = expression(theta[t])) +
  lims(y = c(-1, 5.5)) +
  theme(panel.background = element_blank())

ggsave('~/Documents/Research/boulder/fifth_semester/presentation/theta_noisy.pdf',
       height = 3, width = 3)

#### Single panel plots

gen.cond %>%
  filter(n.pop0 < 50, !low.var, !alpha) %>%
  ggplot(aes(x = gen)) +
  geom_ribbon(
    aes(
      ymin = vbar - 2 * sqrt(vvar / nobs),
      ymax = vbar + 2 * sqrt(vvar / nobs),
      group = interaction(ext.gen, extinct),
      fill = extinct
    ),
    alpha = 0.1
  ) +
  geom_line(
    aes(
      y = vbar,
      group = interaction(ext.gen, extinct),
      colour = extinct
    )
  ) +
  geom_line(
    aes(
      y = vbar,
      group = interaction(ext.gen, extinct),
      colour = extinct,
      alpha = extinct
    )
  ) +
  scale_alpha_manual(values = 1:0, '') +
  scale_fill_manual(values = c('black', 'red'), '') +
  scale_color_manual(values = c('black', 'red'), '') +
  labs(x = "Generation", y = "Additive genetic variance") +
  # facet_wrap(dd ~ paste(n0, lv, sep = ', '), ncol = 4) +
  theme(
    legend.position = 'none',
    panel.grid.major = element_line(colour = 'gray88'),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour = 'black')
  )

ggsave('~/Documents/Research/boulder/fifth_semester/presentation/genvar.pdf', width = 5, height = 5)

gen.cond %>%
  filter(n.pop0 < 50, !low.var, !alpha) %>%
  mutate(fixp = 25 * fixp) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = fixp,
      group = interaction(ext.gen, extinct),
      colour = extinct
    )
  ) +
  geom_line(
    aes(
      y = fixp,
      group = interaction(ext.gen, extinct),
      colour = extinct,
      alpha = extinct,
    )
  ) +
  scale_alpha_manual(values = 1:0, '') +
  scale_fill_manual(values = c('black', 'red'), '') +
  scale_color_manual(values = c('black', 'red'), '') +
  scale_y_continuous(breaks = (0:3)*3) +
  labs(x = "Generation", y = "Loci at fixation (positive allele)") +
  # facet_wrap(dd ~ paste(n0, lv, sep = ', '), ncol = 4) +
  theme(
    legend.position = 'none',
    panel.grid.major = element_line(colour = 'gray88'),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour = 'black')
  )

ggsave('~/Documents/Research/boulder/fifth_semester/presentation/fixpos.pdf', width = 5, height = 5)
