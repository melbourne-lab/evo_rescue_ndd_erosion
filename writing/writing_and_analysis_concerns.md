Concerns on writing and presentation

### Extinction differences

I'm actually not so sure that the long-term extinction risk is that different. Plotting the instatanous extinction risk (i.e., Pr(extinct in time t | extant in time t-1)) shows little difference in extinction probability later in the experiment.

Then again, just plotting the instantaneous extinction probabilities at any time step shows little difference in actual extinction probability... but there is a cumulative difference. Is this plot even compelling?

But, specifically with respect to long-term risk - the logic here is not great. First, popluations are quite large - even if their growth rate stayed lower, they would not be incredibly at risk of extinction at stochasticity. Second, even if the population size *did* fall, the growth rate would increase (by definition of NDD) so extinction would still be less likely. Of course, extinction risk is going to be much greater than for a population exploding to infinite size. Of course, with evolution *specifically in the Ricker model*, faster adaptation means the growth rate grows more quickly. Here NDD-induced genetic erosion would still influence population dynamics by depressing the rate at which the carrying capacity grows. But the case here still seems kind of weak to me.

One thing that may be happening is that we aren't seeing enough extinction happening to observe a difference. Perhaps differences in extinction probability would become more stark if both populations were at greater risk of extinction initially - e.g., d_0 greater, or w^2 smaller. 

###### Solutions to consider

First, proposed Fig. 2 (extinction probabilities) would likely be better served by including more trials. Here, we wouldn't need to save every individual - we could run a greater number of trials and for each trial just save the population size at each step.

We could also increase the extinction risk to see if that makes differences in extinction probabilities more drastic. One logical way to do this is just to increase d_0 - the current d_0 used gives us a w_0 approx 0.77, whereas the proposal specifies w_0 approx 0.7. Using d_0 = 2.75 instead of 2.6 would give us w_0 approx 0.69. Increasing the extinction risk comes with its own downsides (see below).

### Rate of genotypic change

I'm very weirded out by the fact that there is close to zero difference in genotypes between the high and low variation treatments. The Lande expression says that this effect should be quite large! What gives? For a triple-bottlenecked population, which has only something like 40% of the genetic variation of a non-bottlenecked population, the rate of adaptation should be a lot slower. So why isn't it?

A (very) likely candidate is selection. Less-fit populations go extinct; this raises the mean genotype of the remaining (surviving) populations. Perhaps we are just stuck in a paramter-range where the degree of selection is just enough to offset the difference? This may explain why there are divergences in the later time steps, when extinctions are much more rare. 

If selection is what is causing this disparity, is this figure even worth reporting? Perhaps it is - it's worth noting that we are deviating from the genotypic predictions of G&H. On the other hand, it makes the bit about NDD slowing down adaptation a little more difficult to believe. There is some confounding going on here.

Note that if looking at only the first two timesteps in simulations (here there are close to zero extinctions occurring) that we do see some weak evidence that there is a difference in the rate of adaptation. Consider the following: open up script `simulations/analysis_results/manuscript_figures/generate_manuscript_figs.R` and load in the data frame `alldemo`. Now run the following command

```
alldemo %>% 
	filter(gen %in% 1:2) %>% 
	ggplot(aes(x = gen)) + 
	geom_line(aes(y = gbar, group = interaction(n.pop0, bottleneck, ndd), colour = bottleneck))
```

(requires `dplyr` and `ggplot2`). It's sloppy, but it shows that the bottlenecking (variation) does have a slight effect on the slope of the lines (one can do linetype by initial population size to see additional effects from slope). 

###### Solutions to consider

Remove extinctions from consideration - run a separate batch of trials where populations are initially large enough that they do not go extinct. This perhaps can be reported in supplements? Do this to confirm the shape of the change in genotypes over time - is it geometric? Is it sub-geometric? To handle size-issues, consider the approach from above summarizing within a trial and saving those results.

Running in a different part of parameter space may also cause th lines to converge. Right now, potential effects may be masked by offsetting forces (selection versus low variation); running in a different part of parameter space may cause change the magnitude of one or both of these in a way that sums up to a bigger difference.

A thing to consider: if we up the extinction probability (as suggested above) we further confound any differences induced by low population variation. 

**Update August 20** Trying to run trials with a large initial population size actually will erase diferences the model is trying to show. Large initial population size means populations will not be as small at their lowest points, meaning the deterioration will be small. One solution is perhaps to decrease the genetic distance (initial load, d_0).  

Another problem is that our mean population change deviates substantially from the G&H estimates and I am not sure why. Obviously the model is dfferent, but that difference is apparently huge. Specifically in the high variation models, genotypic change is much slower than what is predicted in the models. It also has a strange shape - not very geometric looking in the first several timesteps. I checked the estimates going in to the Lande/G&H expression and I am pretty sure they are correct. What gives?

### True sampling variation

Not a novel insight, but the bottlenecked trials are not independent of each other (big effect of population line). This should increase the sampling variation and requires some thought for calculating standard error. Is it really as simple as thinking hierarchically - summing together within- and among-population variation for each time step? Ask Brett.

###### Solutions to consider

Ask Brett.

