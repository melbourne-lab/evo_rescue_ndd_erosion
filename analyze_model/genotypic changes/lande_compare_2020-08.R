# Another script where I compare phenotypic predictions with expectations in Lande, 1976.
# August 20 2020

set.seed(4199)

parpar = data.frame(n.pop0 = 100,
                    end.time = 15,
                    init.row = 1e4,
                    n.loci = 25,
                    w.max = 2,
                    theta = 3,
                    wfitn = sqrt(1 / 0.14 / 2),
                    sig.e = 0.5,
                    pos.p = 0.5,
                    alpha = 0.000)

testsim = sim(a = c(1/2, -1/2), params = parpar)

table(testsim$gen)
table(testsim$gen) %>% plot()

tsgen1 = testsim %>% filter(gen %in% 1)

tsgen1 %>% mutate(zwt = z_i * w_i / mean(w_i)) %>% summarise(zbar = mean(z_i), zwbart = parpar$theta - mean(zwt), sig2 = var(z_i))
parpar$wfitn^2 / (parpar$wfitn^2 + 0.908) * (parpar$theta - 0.07168)

(
tsgen1.summ = testsim %>%
  filter(gen %in% 1:2) %>%
  group_by(gen) %>%
  mutate(d_i = parpar$theta - z_i) %>%
  summarise(dbar = mean(d_i),
            sig2 = var(g_i),
            dwbar = mean(d_i * w_i / mean(w_i)),
            predzbar = dbar * with(parpar, (wfitn^2 + sig.e^2) / (wfitn^2 + sig.e^2 + sig2)))
)

(
  tsgen.summ = testsim %>%
    group_by(gen) %>%
    mutate(d_i = parpar$theta - z_i) %>%
    summarise(dbar = mean(d_i),
              # sig2 = var(g_i),
              dwbar = mean(d_i * w_i / mean(w_i)),
              preddw = dbar * with(parpar, wfitn^2 / (wfitn^2 + var(z_i))),
              preddbar = dbar * with(parpar, (wfitn^2 + sig.e^2) / (wfitn^2 + sig.e^2 + var(g_i))))
)

plot(tsgen.summ$preddbar[1:14], tsgen.summ$dbar[2:15])
abline(0,1)
# as we've seen, prediction is biased
plot(dwbar ~ preddw, tsgen.summ)
abline(0,1)
# but the dwbar is correct... this means Lande's z_w(t) is correct

# Maybe heritability in our model is not actual h^2?
# Test heritability in model tomorrow.
# Also wonder if the sex does anything - should increase variance but not bias results...
