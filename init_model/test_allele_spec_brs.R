n.trials = 1000

liszt = vector(mode = 'list', length = n.trials)

set.seed(1581823)

for (trial in 1:n.trials) {
    liszt[[trial]] = sim(
        a = c(-1/2, 1/2),
        params = data.frame(end.time = 15,
                            init.row = 1e4,
                            n.loci = 40, 
                            n.pop0 = 20,
                            w.max = 2, 
                            theta = 2.6, 
                            w.sqd = 1 / 0.14)
    )
    print(trial)
}
    

all.sims = unroller(liszt) %>%
    select(trial, names(.)[!grepl('trial', names(.))])

sumry = all.sims %>%
    select(-c(r_i, fem)) %>%
    gather(key = loc.copy, value = val, -c(trial, i, b_i, w_i, gen)) %>%
    mutate(loc = gsub('[ab]', '', loc.copy) %>% as.numeric()) %>%
    group_by(trial, gen, loc) %>%
    summarise(n = length(unique(i)),
              m = mean(val),
              v = var(val, na.rm = TRUE),
              b = mean(b_i),
              w = mean(w_i),
              fixp = all(val > 0),
              fixn = all(val < 0)) %>%
    group_by(trial, gen) %>%
    summarise(mean.var = mean(v, na.rm = TRUE),
              pop.size = mean(n),
              n.fixp = sum(fixp),
              n.fixn = sum(fixn),
              mean.b = mean(b),
              mean.w = mean(w)) %>%
    group_by(gen) %>%
    summarise(n.trials = n(),
              v.bar = mean(mean.var, na.rm = TRUE),
              v.var = var(mean.var, na.rm = TRUE),
              n.bar = mean(pop.size),
              n.var = var(pop.size),
              b.bar = mean(mean.b),
              b.var = var(mean.b),
              w.bar = mean(mean.w),
              w.var = var(mean.w),
              p.fixp = mean(n.fixp) / 40,
              p.fixn = mean(n.fixn) / 40)

sumry %>%
    write.csv('model_building/eg_1k_n0_20.csv',
              row.names = FALSE)

g.n = sumry %>%
    mutate_all(function(x) ifelse(is.na(x), 0, x)) %>%
    ggplot(aes(x = gen)) +
    geom_line(aes(y = n.bar)) +
    geom_ribbon(aes(ymin = n.bar - 2 * sqrt(n.var / n.trials),
                    ymax = n.bar + 2 * sqrt(n.var / n.trials)),
                alpha = 0.2) +
    scale_y_log10()
g.n

g.r = sumry %>%
    mutate_all(function(x) ifelse(is.na(x), 0, x)) %>%
    ggplot(aes(x = gen)) +
    geom_line(aes(y = r.bar)) +
    geom_ribbon(aes(ymin = r.bar - 2 * sqrt(r.var / n.trials),
                    ymax = r.bar + 2 * sqrt(r.var / n.trials)),
                alpha = 0.2)
g.r

g.f = sumry %>%
    ggplot(aes(x = gen)) +
    geom_line(aes(y = p.fixp)) +
    geom_ribbon(aes(ymin = p.fixp - 2 * sqrt(40 * p.fixp * (1 - p.fixp) / n.trials),
                    ymax = p.fixp + 2 * sqrt(40 * p.fixp * (1 - p.fixp) / n.trials)),
                alpha = 0.2)
g.f # hmm... why does this start off with non-zero amount of fixation?
    
g.v = sumry %>%
    ggplot(aes(x = gen)) +
    geom_line(aes(y = v.bar)) +
    geom_ribbon(aes(ymin = v.bar - 2 * sqrt(v.var / n.trials),
                    ymax = v.bar + 2 * sqrt(v.var / n.trials)),
                alpha = 0.2) +
    scale_y_log10()
g.v

g.n
g.r
g.f
g.v

# What is the mean allele frequency?