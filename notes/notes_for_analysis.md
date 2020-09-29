# Ideas for analysis of simulations:

### Time to extinction (09/21)

Will survival analysis be required here...? 

Intuitive and expected result: under NDD, shorter time to extinction (all else equal).

But - keep in mind there will be *more* extinctions. Fig. 2 will show where these extinctions will be happening. 

What if time to extinction is higher? Could be that more extinctions are happening *later* in the simulation. The temporal aspect here coudl be a little bit confusing.

### Look at relative influences of above/below average (median?) genotypic change in first several time steps and risks on extinction (09/21)

Idea here comes from Peniston et al. 2020 (I think). Consider: how predictive of ultimate extinction risk are above/below average changes in the first several timesteps? E.g., say you have above-average genotypic change, but below average population growth rate. Which is more predictive of final outcome?

Likely problems here: (1) these will be correlated (genotypic change -> phenotypic change -> population growth; genetic variation -> genotypic change) although degree of correlation is unclear, stochasticity may weaken the correlations (2) what time step to use? Time steps 2-4 seem to me to be good times. Is temporal autocorrelation an issue here?

Issue (2) makes me think - there's a difference between the *overall average* and the average *in that time step*. Change within a time step should be Markovian, no? E.g., a lower-than-average change in some of these should be independent of a lower-than-average decrease in the very next time step. The exception here should be genetic variation I think... Hmm. This would at least perhaps take care of autocorrelation. You could get a sequence of something like "LHLLL" (e.g., lower higher lower lower lower) which could be helpful. Think about this more.

### Inbreeding coefficient (then look at inbreeding vs. extinction) (09/21)

It would be good to see how much inbreeding contributes both to the overall genetic load and also to extinction risk. The pedigreed sims would be essential here. I don't know a neat way to calculate inbreeding coefficients.

### In general: growth rate versus genetic diversity (09/21)

This seems like the pair of variables with the weakest correlation, but there is some connction. More diversity -> more genetic change -> more phenotypic change. Because the phenotypic change is expected to be positive, more phenotypic change should mean greater increase in population size. But, there's stochasticity at every step of this process. This is interesting because the premise of these models is that increased genetic diversity directly contributes to population growth rates. It's also central to the vortex.

### Completng the loop: pop size/NDD on variation in genotypic change (09/29)

Lower population size (and through this means, NDD) may not necessarily change the *mean* rate of phenotypic/genotypic change, but it *would* increase the sampling variance in genotypic change.

What is the net effect of this change? The sampling variation in phenotypic change is not biased, i.e., the expected phenotypic change *should* still be the same. But, the effect of this phenotypic change on fitness may be biased.

Say the current genotypic distance from the optimum is $d$, and the expected genotype is $Dd$. Under our parameterization $Dd < d$. The expected change in fitness, then, is positive as exp(-Dd) < exp(-d). Note that the expected change in fitness is by definition the change in fitness with in a deterministic model, i.e., a model with no variance in the genotypic changes. However, say that there is $e > 0$ variance in the genotypic change. If  $W''(d) < 0$ (i.e. if the fitness landscape is concave down at $d$), then for any $0 < x < e$, $W(d+x) - W(d) < W(d) - W(d-x)$, i.e, a population with greater-than-average change will have a marginal fitness advantage that is less than the marginal fitness disadvantage of a population with less-than-average change. The result will be flipped if the concavity is positive.

Because added variation in the phenotypic change changes the expected change in fitness, it will change the population growth rate. The concavity of the fitness landscape ($W(z)$) determines whether the change is positive or negative. For a Gaussian fitness landscape, this means that sampling variation would *increase* fitness when populations are far from the optimum but decrease it when they are closer to the optimum. 

This alone would not necessarily change the rate of adaptation, just the (instantaneous) population growth rate. Perhaps a more evolutionarily-relevant and straightforward way that sampling variation would affect the rate of genotypic variation is to cause genetic drift, which will lower the rate of adaptation. Argh... but this is what we've been having so much trouble demonstrating!

