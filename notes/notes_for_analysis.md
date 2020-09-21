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

