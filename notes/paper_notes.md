### Reznick, 2016. Hard and soft selection revisited: How evolution by natural selection works in the real world. Journal of Heredity, 2016, 3 - 14.

##### Haldane and the "classical hypothesis"
- Models assume there is one best genotype (w/ highest fitness)
	- difference between this fitness and (mean) observed fitness is *genetic load*
- Assume that adapting to optimum occurs through slow allelic replacement (one at a time)
	- this is the "classical hypothesis"
	- somehow (?) it follows from this that genetic variation is limited (presumably because it is selected against)
		- although, we know that there is large natural variation

##### Hard versus soft selecton
- Haldane's selection is hard
- Soft selection: environment has vacancies to be filled, w/ genotypes competing to fill them; it will be filled either way
	- this somehow could explain genetic variation...?
	- in a sense these models assume population growth is exponential and ignore resource limitations

##### Guppy example
- Naturally occurring high predation (HP) and low predation (LP) environments
- Predation risk higher for larger individuals
- In HP environments, smaller reproductive individuals, more reproduction
- However, no evidence of age-dependent mortality changing across predator-environments (i.e., no age x predator interaction on mortality)
	- the population ecology (size distribution) explains this: in LP environments, resource stress induced by unproliferated growth not present in HP environments
- Bassar et al., 2010, 2012, 2013, experimentally found that HP guppies have higher fitness than LP at low population density, but not at high density
	- i.e., adaptation to density influences adapting to LP environments, because guppies in LP environments are less sensitive to density

### Shaw and Shaw, 2014. Quantitative genetic study of the adaptive process. Heredity, 112: 13 - 20.

##### Fisher's Fundamental Theorem of Natural Selection
- Change in mean fitness is additive genetic variance divided by mean absolute fitness
	- V_A(W) / \bar{W}
	- doesn't account for mutations (short time scale)
- Idea hasn't been used empirically very much
- Charlesworth, 1987, suggests VA should be slight at equilibrium

##### Quantitative genetic models
- Turelli and Barton, 1984 demonstrate that with strong selection, dist'n of breeding values is approximately normal (we should confirm)
- Under Gaussian selection (like our models), phenotype additive environmental and genetic (e + a = z)
	- Distribution of genotypic effects on absolute fitness, E(W|a)
	- check appendix for PDFs of W?
- "the inclusion in the model of environmental variance as a contributor to variation in the trait causes the genotypic distribution for fitness to be shifted toward intermediate values. In that case, no genotype has an expected fitness that coincides with the maximal fitness, and genotypes with low fitness values in the absence of higher variation have higher expected fitness when it is included"
	- so environmental variance speeds up when far away but slows down when you get closer (because it's adding noise)
- I'm really lost by the second paragraph of the "Fitness under Gaussian selection" paragraph...
	- What is V_A(W) and what is V_G(W)?
		- V_G(W) is the variance of the PDF of W, f(W)
		- V_A(W) is the additive genetic variance for fitness ?
	- Fig 2: increasing sig_a (additive genetic variance) increases V_G(W)
- Pease et al., 1989, has a population experiencing gaussian selection with a moving target
- Burger and Lynch, 1995, look at persistence of populations using environmental change
- Ronce et al., 2009, has a model very similar to ours, as well as a framework for incorporating inbreeding
	- this model allows genetic variance to evolve freely
	- includes mutations
	- here, a rapid increase in trait V_G(W) with a subsequent change in sig^2_a... is this due to the moving target? What is going on here.
	- also calculate V_A(W) somehow... but how?
	- take home of this is that V_A(W) can increase when the selection regime changes

### Ronce, Shaw, Rousset, and Shaw. Is inbreeding depression lower in maladapted populations? A quantitative genetics model. Evolution, 63(7): 2009.

##### How do ecological scenarios and influence inbreeding depression?
- A widely-accepted model of inbreeding: two alleles, one is unconditionally deleterious, at least partially recessive
	- see Charlesworth and Charlesworth, 1999
	- Acc'd to this model, inbreding should purge deleterious aleles
- These models directly link genotype to fitness, no traits, unlike qGen
- This paper looks at inbreeding w/ stabilizing selection on quantitative traits
	- see Lande and Schemske, 1985, for something simpler

##### Model
- trait, z, has additive action of n loci, diploid, panmictic
	- z = mu + sum_i=1^n (am_i + af_i) + e
	- am_i, af_i are allelic effects in herited from parents, mu is mean
- trait under Gaussian selection
- I = 1 - (\bar{W}_s) / (\bar{W}_r)
	- I = inbreeding depression
	- selfing versus random mating fitesses (selfing: self is parents?)
- Analytic model follows Lande, 1976
	- here, an expression for inbreeding depression! interesting

Come back to this paper later. This is going in some directions beyond what we need right now.

### Osmond and de Mazancourt, 2013. How competition affects evolutionary rescue. Phil Trans R Soc B 368: 20120085.

##### Competition
- Reduces population abundance (extinction more likely)
- But, also increase selective pressure, increase rate of evolution (refs herein)
- What is the influence of competition (intra and interspecific) on rescue?

##### One species model
- Asexual population: density- and frequency-dependent growth
- Continuous traits? Continuous time
- Growth of # individuals w/ a trait is logistic
	- features per capita intrinsic growth rate
	- carrying capacity for each trait is a defined by gaussian
- Intraspecific competition: same trait value compete most strongly
- Mutation is rare, s/t demography is at equilibrium, evolution is slow
- Population initialized at constant trait value, perfectly adapted, reaching optimal carrying capacity
	- environment changes, carrying capacity changes
- N_c is extinction threshold
- for a given N_c, get an expression for expected time (t_R) to surpass N_c
	- t_R increases with increasing environmental change and N_c
	- unimodial w/r/t environmental tolerance (sigma_K)/selection strength
		- tol too high, and population is never at extinction risk
		- tol too low and population immediately extinct
- One term in here is the analogue to the rate of change in the trait value from G&H
	- in O&dM model, evolution stops with absence of mutations

##### Discussion
- Competition which increases the strength of selection may speed up adaptation... where is this coming from? Ref Schluter et al., 1985.
- Refs: "intraspecfic competition often has relatively little impact on selectiv pressures" (refs)
- Impacts of NDD:
	- Boulding and Hay: populations stay small, where NDD doesn't act.
	- Chevin and Lande: NDD *makes* populations small, where extinction is more likely
- This model assumes "evolution is slow, and hence, on the time-scale of evolution, populations are always at carrying capacity"
	- this paper is different from our study b/c we're interested in *rapid adaptation* and the stochasticity is important
	- "we assume greater abundances lead to faster evolution because they lead to greater mutational inputs" - population size is important for mutation
	- "intrinsic growth rate is therefore irrelevant for evolutionary rescue in our mode" (due to the demographic equilibrium)
		- this is evolutionary rescue under K-selection rather than r-selection

So... this paper is all about how mutations drive rescue and how that interacts with density dependence. Model is fundamentally different because populations are assumed to be at carrying capacity for most of the time... but then how to explain Fig. 2? Not sure.

### Vinton, A.C., and Vasseur, D.A. Evolutionary rescue is determined by differntial selection on demographic rates and density dependence. bioRxiv preprint.

##### Evolutionary rescue
- Four factors influencing potential for rescue (Bell and Gonzalez, 2009)
	- Initial population size (N_0)
	- Genetic variability (standing variation and mutations)
	- genetic variability due to dispersal (interesting...)
	- extent/severity of environmental change
	- (n.b., three of these are in G&H's equations)

##### Demography
- Models typically only focus on $r$, ignore separate births nad deaths
- e.g., Holt, 1990: populations with high birth and death vs. low birth and death
	- these have the same $r$, but different rates of turnover, mutation, etc.

##### Model
- Continuous time, ndividual-based logistic growth
- N' / n = b(gamma0, gamma1) - d(gamma0, gamma1) where
	- gamma0: density independent contributions
	- gamma1: density dependent contributions
	- b(), d() depend on rates, carrying capacity K
	- carrying capacity assumes all individuals well adapted
- environmental effect: mu_0(t,f), varies sinusoidally over t
	- e_mu = abs(mu - mu_0) is a measure of maladaptation
- Birth and death models:
	- Case 1a: environment alters b(), density independent modification
	- Case 1b: environment alters b(), density dependent modification
	- i.e., the degree to which b() is lowered does/does not depend on environment
	- Case 2a, 2b: density modifies d(), density independent/dependent
	- Density dependence: only see decreases/increases in rates when pop'n size is large
		- e.g., little maladaptation at size one
- Stochastic implementation via the Gillespie algorithm as birth/death process
	- good description of Gillespie here!
	- includes mutation, randomly drawn from U(-0.3,0.3), added to parential trait

##### Results
- r is differece of birth and death, but stochasticity is determined by sum (?)
	- i.e., variance is due to the sum of variances in birth and death?
	- refs herein
- death models have higher demographic stochasticity
	- there's a reason here but I don't understand it.
- due to model construction, density dependence means *less* variable population size
	- (by construction this means you maintain maladapted individuals at low size)
	- notably, this is is relative to a DI-*effect* on rates, not on DI itself
	- this is a case with phenotype x growth rate interaction
- increasing persistence time with a faster-changing environment?

(didn't quite finish, but worth re-reading)

### Ferriere, R., and Legendre, S. 2013. Eco-evolutionary feedbacks, adaptive dynamcs and evolutionary rescue theory. RTSB 368: 20120081.

##### Framing
- The "optimization principle"
	- evolution should optimize phenotype by maximizing fitness
- Two approaches to evolutionary rescue
	- Pop gen modeling: how do mutations reduce extinction risk
	- Quantitative genetics: how does a population track a moving optimum
	- Differences:
		- Environmental change (abrupt vs. gradual)
		- Threat (demographic defecit vs. risk of defecit)
		- Rescue pathway (achieving demographic balance vs. avoiding imablance) (?)
- Haldane outlined that adaptation doesn't necessarily always improve population performance

##### Frequency dependence and the optimization principle
- Frequency dependence means departures from optimization
- Frequency dependence does not occur iff the trait effects fitness in a 1D monotonic way
- Realistic eco-evo feedbacks are expected to rely on freq. dependence
- "Even when the adaptive process does optimize, phenotypic evolution maximizes the population's intrinsic growth rate [...] or size [...] only if the eco-evolutionary feedback loop posseses special features regarding how density dependence effects transmissions in the life cycle, and how traits influence invasion fitness" refs
- "And in the rare cases where adaptive evolution gradually tends to improve population viability, one cannot exclude that the process ends with the population's sudden demise" refs

##### Adaptive dynamics framework
- Three ingredients:
	- Descrption of an individual phenotype (e.g., quant gen)
	- Relates the traits to population properties (e.g., growth rate)
	- model of trait inheritance
	- n.b. we have all three
- "Evolutionary singularity" - phenotype for which local fitness gradient vanishes
	- a fixed point...? described as being stable or unstable
- Pop gen and Q-gen models ignore how population changes can vary selection pressures acting on heritable variation
	- the adaptive dynamics closed-loop framework suggests that only under "exceptional circumstances" is pop'n size/viability maximized at an attractor
- "evolutionary trapping" - an evolutionary attractor in the unviable region
	- eolutionary suicide: trait substitution driven by mutation and selection in a constant environment pushes a population to inviability
	- in short, populations can tank their own viability through evolution (inbreeding an example?)
- Examples of trapping
	- altruism...?, dispersal behaviors, etc.
	- (none of these look relevant to our fitnesses)

Beginning of this paper is useful for the framing insisting that including population's effects on selection itself is important; G&H model, e.g., does not include the fact that population size may influence the adaptation itself (fixation or stochasticity). But most of this paper is talking about "threats of evolution itself" on population viability (see: conclusions) but this isn't what we want unless our model becomes much more complicated (e.g., including a dispersal-competition trade off).

### Kopp and Matuszewski, 2013. Rapid evolution of quantitative traits: theoretical perspectives. Evolutionary Applications.

Some good info: *"Haldanes"* - a proposed measure of rate of phenotypic change; one phenotypic standard deviation per generation. There is some empirical evidence that evolutionary rates above 0.1 haldanes are not uncommon (references herein).

##### Overview of genetic models
- Environmental change: can be sudden (like our models), gradual unidirectional change, or random fluctuations of the optimum around a stationary optimum or linear trned.
- Phenotypic variance: sigma2_g is genetic variance (e.g., additive genetic variance, like our models), sigma2_e, summed to get sigma2_p (phenotypic variance)
	- heritbility is sigma2_g / sigma2_p
	- Measure phenotypes in units of environmental variance - this allows sigma_e = 1.
- The Lande equation (Lande, 1976) has delta zbar_t = sigma^2g \beta_t, where \beta_t is derivative in log fitness w/r/t zbar
	- i.e., rate of phenotypic change is proportional to the variance and how much the (log) of fitness changes w/ pheno change
	- divde both sides by sigma^2_g to get haldanes
- Phenotypic plasticity is characterized as a reaction norm, i.e., how a phenotype changes with a given environmental variable
	- if different genotypes produce different reaction norms, selection can act on plasticity
	- typically models feature linear rxn norms w/ quantitative slopes and intercepts

##### Overview of population dynamic models
- typically N_{t+1} = \bar{w}_t N_t
	- \bar{w}_t, mean fitness, reduced by *selection load* (distance from phenotypic optimum) and *standing load* due to phenotypic variance (n.b. this shows up in the G&H equations!)
	- Note: population dynamics depends on relative fitness but evolutionary change depends on the fitness gradient/selection landscape

##### Syntheses

###### Sudden environmental change with one trait
- Gomulkiewicz and Holt, 1995 case
- Is evolution fast enough for change to occur to move a population out of the \bar{w} < 1 zone?

###### Gradual environmental change with one trait
- Lynch et al., Lynch and Lande, etc. studied cases where optimum increases at linear rate k
- Adaptation is initially slow, leading to an increase in selection load, which increases the rate of adaptation until dynamic equilibrium is reached where phenotypic rate of change is proportional to k
	- persistence depends on the selection load, which is proportional to k, so there is a critical k value which is the maximal rate of change a population can handle
	- persistence requires positive growth at all times? A population which has "fallen behind" in the race w/ the optimum will have no chance to catch up.
- critical k depends on additive genetic vriance and (square root of) the max. population growth rate
	- fastest rate of adapatation for intermediate width of fitness function; extremes are bad due to the two types of selection load (this is interesting!)

###### On additive genetic variance
- All models assume it is true
	- this approximately true on short timescales, but not on longer timescales
- Before environmental change, populations may be assumed to be a mutation-selection-drift balance
- Variance may increase in larger populations, as there are more individauls and therefore more mutations
	- But this increase takes time and may come too late to be of use
- (No mention of potential decreases as selection is occurring...)

###### Fluctuating selection
- Uncorrelated (white noise) fluctuations is akin to recurring sudden-change but won't create long-term patterns (other than increased extinction)
	- generally produces higher selection loads b/c adaptations in one generation may not be helpful in the next
- Autocorrelated fluctuations are more similar to a gradual-change scenario which a population can track

###### Multiple correlated traits
- Genetic correlations will bias the phenotypic response towards the leading eigenvector of the covariance matrix
- Some effects here where a population may have permanently maladapted traits or trait combinations (flying kite effect)

##### Max sustainable rates of evolution
- Barrett and Hendry, 2012, conclude that most rates of change are below 0.1 haldanes (this is based in part on a metaanalysis)
	- They include the warning that "many unrealistic assumptions, such as perpetual persistence under constant environmental change"
	- likewise "critical rates for natural populations over time frames of conservation interest could be very different"
- Detecting a difference of 0.1 sds between two populations requires a lot of samples...

##### Adaptive potential
- Check Gomulkiewicz and Houle (2009) - demographic constraints requiring critical amounts of genetic variance IOF rescue to occur
- For multiple traits, the G matrix can have eigenvalue zero, in which case selection won't happen
- Ploidy: advantages of diploidy may include heterozygote advantage (smaller phenotypic effects in heterozygotes)
	- This heterozygote advantage would maintain polymorphism and increase standing variation available for adaptive potential
- Theory suggests that more complex genetic architectures typically hinder evolvabilty as mutations are more likely to have negative pleiotropic effects
	- modular genetc organization can help reduce "effective complexity"
- More adaptive potential in populations in fluctuating environments
	- more genetic variance due to pretesting in environment
	- species may evolve plasticity and/or flexible architecture

Some thoughts in here about *adaptation in space* too (incl. gene flow).

##### Other species
- Presence of competitors reduce population size (and therefore genetic variance and mutational input) decreasing adaptive potential

*Conclusions* there's a great review in here of different theoretical results. It's worth keeping in mind shifting targets (theta) in mind for our models. The point about additive variance being constant on short time scales is worth pointing out, as it is what we are testing. 

### Gomulkiewicz and Holt. 1995. When does evolution by natural selection prevent extinction? Evolution, 49: 201 - 207.

##### Motivation
- Pease et al., 1989, various work by Lynch, follow populations tracking a continuously moving environment; here look at a single, abrupt shift.
- A feature of this model (SN - not sure how this is different from other models) is that populations may fall below a critical population size during which it is vulnerable to stochasticity-induced extinction
- Prior works (Goodman, 1987, Lande, 1993) talk about extinction thresholds
- Adapting to a novel environment as a race between two processes: demographic and evolutionary (see Maynard Smith 1989)

##### Model
- Simple population growth w/ per-capita growth rate \bar{W}_t
- Evolutionary polygenic quantitiatve trait z
	- Gaussian selection w/ strength w (wsq)
	- Phenotypic variance P, mean d_t
- In simplest case, N_t = N_0 \prod_{i=0}^{t-1} \bar{W}_i
	- but, with density selection, N_{t+1} = N_t \bar{W}_t g_t(N_t)
	- typically g_t(N_t) < 1, such that the density independent case establishes an upper bound on persistence capacity
- In this model, phenotypes geometrically approach optimum at rate described by Lande, 1976
	- From Lande, an expression of the mean fitness
	- n.b. this equation features the selection load (the exponential term) and the standing load (the \hat{W}) term
	- (interesting trade off... increasing your variation means your selection load decreases more quickly, but you have a lower ceiling due to standing load)
- Analytical derivation of times where population size is equal to the extinction threshold
	- can we get insight here from a non-constant g_t()?

##### Analysis
- Fig 3: initial conditions below the curve persist deterministically, while those above curve deterministically fall below N_c
	- i.e., lower curve is, the larger the danger zone is
- Fig 4: demonstrates that there is a threshold of initial maladaptation (\beta^*) below which an evolving population will not deterministically fall below the extinction threshold
	- For larger degrees of adaptation, the difference between an evolving case and a non-evolving case are small (i.e., similar times to reach the extinction threshold)
	- This may be the "near-threshold" effect of evolution - only close to a critical value does evolution influence the time to the threshold
- Initial degree of maladaptation and initial population size (both properly scaled) determine time to threshold
- What influences extinction risk in this model?
	- Influences on evolutionary capacity: k and \hat{W} (resp. the rate of phenotypic change and the maximum fitness rate)
	- The initial state of the population: \beta_0 and \nu_0 (resp. scaled squared degree of maladaptation and population distance to optimum phenotype)
- Fig 5: the "time at high risk" (time below the threshold) increases linearly with maladaptation, but only changes slightly with initial density (b/c most of the change happens for N_c/N_0 approx. 0)

##### Discussion
- Interesting connection made here to colonizers in a new environment: a large, high max. growth, and/or well-adapted pool of propagules is more likely to have a successful colonization
- One thing mentioned here which I (SN) had never thought of is that the max growth rate, perhaps intrinsic to a species? (or life history, at the very least a population) is something that also inhibits its growth
- Another reason the expressions here are an upper bound: evolutionary forces may increase maladaptation (mutations, drift) or reduce genetic variance (drift)
- For assisted migration, (1/t_E) is a rough guide for frequency of migration episodes (relevant to the drift)

##### Points for discussion, future thought
- This does make figuring out that "extinction threshold" via simulations more appealing (first step is to read Lande, 1993)
- This has likely been done in some form, but, thinking about the trade off between selection load and standing load, there needs to be some sort of optimal rate correct? Does the optimal rate for population growth guarantee the highest probabilty of rescue/persistence?
- One thing this model assumes is that the population dynamics and genetics operate independently. I think this may arise from the fact that the quantitative genetics model has infinitely many loci. But a finite-locus model I think will capture losses in genetic diversity more faithfully.
- In an infinite-locus model, is there a functional form of g(N) which will give an analytically tractable expression for N_t?

### Barfield and Holt, 2016. Evolutionary rescue in novel environments: towards improving predictability. Evolutionary Ecology Research, 17: 771  786.

##### Framing
- Good description of rescue in first paras
- "Whether evolutionary rescue happens or not depends upon both demographic and genetic factors (Gomulkiewcz and Houle, 2009)
- "Not all selective processes enhance overall population fitness" e.g., strong frequency dependence, drift, stochasticity when populations are smaller, variable amount of genetic variation
	- "This raises challenging questions of how predictable evolutionary rescue might be"
- Genetic architectures, including "polygenic models where may loci contribute..." (H&G 2004, Holt et al. 2005, Knight et al. 2008, Gienpp et al., 2013)
- This paper: compares original G&H95 to a heuristic model w/ changing heritability (see Lande and Barrowclough, 1987 and Ellstrand and Elam, 1993)

##### Novel model
- Small populations are expected to lose genetic variation (drift and inbreeding) and also by reduced heritability [n.b. connection here not clear or explicit]
- G&H modified such that h^2 = h^2_{max} N / (N + N_{0.5})
	- N_{0.5} is population size where heritability is half of max (h^2_{max})
	- Fig. 1 shows that with increasing N_{0.5}, lower depth, longer to depth
	- [n.b. this means that heritability can recover? mutations?]
- Fig. 2: extinct populations - heritability, maladaptation, etc. still track over time...? (I suppose "extinct" simply means below an extinction threshold)
- "even though fitness always increases, it is possible for it to be unable to reach 1 in the limit, if the slowdown in evolution is sufficiently rapid [...]. In this case, the population size asymptotically approaches 0. This model illustrates one possible cause of the 'extinction vortex' identified by Gilpin and Soule (1986)'
- Varying initial population size also has dramatic effect: lower size increases extinction risk
	- "This model thus reveals a kind of bi-stability - initially small populations decline inexorably towards extinction, whereas larger populations with the sae population parameters eventually rebound in numbers and persist" (see Figs. 3, 4)
- "This model should be viewed as a heuristic exercise, since the relationship between population size and genetic variation is unlikely to be as tightly bound as we have assumed."

##### Novel IBM
- Based on Burger and Lynch (1995); see Holt et al., 2003
- $n$ diploid loci, randomly segregated w/ no dominance or epistasis
	- mutations at rate $n\mu$ (add gaussian RV to allele)
	- alleles take any real value...?
- demography: hermpaphroditic, stochastc mate choice, deterministic number of offspring ($f$)

###### Analysis of IBM
- Is genetics at the time of environmental change important?
	- Initialized 100 populations, simulated for 1000 generations, then used these in a changed environment
	- Compared probability of persistence with a null distribution to see if initial genetics had influence on probability of persistence
	- It did (i.e., initial genetic makeup influences persistence risk)
- Genetic mean and genetic variance of starting population: significantly influences persistence probabilities but low R2
	- n.b. initial mean was consistent, low variation
	- mean + 2*sd (upper tail of the genetic distribution) was significant, R2 approx 68%
- Extinction risk over time: most extinctions occurred in initial stages of simulations
- Figs. 9, 10 shows pop. size and heritability over time for extant (>100 time steps) and extinct (surviving to time step 25 but not to 100)
	- This is good - shows average trajectory of each
	- Heritability - this doesn't look right at all. Why does heritability increase in extant populations? Mutations? They say it's "complex"...
	- "Therefore, failure of evolutionary rescue could occur in part because sometmes heritability does drop at low population size"
	- "nitially, however, in contrast to what was assumed in the determnstic model, heritability actually rises after the environmental change, even though population size is declining. A likely reason for this is that there were alleles segregating at low frequency in the pre-disturbance environment that now become advantageous in the new envronment. As their frequency rises, there can be a transient increase in genetic variance, as since it takes a while for previously favored alleles to be selected out of the population. But then genetic variance declines again as the newly favoured alleles increase towards fixation (without getting there) in the new environment."
- "The bottom-line of the individual-based model is that evolutionary rescue becomes more predictable if one includes information about several aspects of the entire distribution of genotypic values, not just the mean or variance."

This paper presents a two models that are quite similar to ours. The first model has a heuristic (non-mechanistic) heritability but is otherwise deterministic. The second is an IBM with a couple of sources of stochasticity and a finite number of loci. The second genetics in this model are smilar to ours, except in this model the alleles take on any real value rather than ours. Ours might have more analytic tractabilty, although we don't exploit it. We also can explore genetic aspects of this model a little bit more (fixations). Notably, the genetic model here features mutations, which allow heritability to increase.

They highlight feedback loops and the extinction vortex here. This removes some of the novelty from our models. Shoot. Although, we can still mention that NDD increases extinction risk by steepening the vortex. There are some similar results for us to think about and compare with (e.g., mutations happening during the first several timesteps).

The importance of initial genetic variation is highlighted here. Using mean + 2sd is a good idea. One thing it makes me think is that there's some confounding with the bottleneck - initial population in bottlenecks whereas all trials are independent in the non-bottlenecked populations.

### Lyberger et al. 2020. Evolution in response to extreme events: good or bad? bioRxiv preprint.

##### Framing
- Pulse versus push disturbances (instant vs. sustained)
- This paper: look at how severity and duration of pulse events (?) influences population size and extinction risk

##### Model
- Individual based
- Gaussian fitness defined identically to ours
- Rndom mating pairs, poisson number of ofspring with mean $\Lambda$
- Habitat with at most $K$ individuals, randomly sampling $K$ to survive
- Quantitative genetics population: infinite loci, normal draw with segregational variance (constant) but genotypic variance does change
- Extreme event of length $\tau$; increasing optimum trait value (linearly? instantaneously?)
- Analytical work - recursion relationship for genotypic change and change in genotypic variation
- Analytically found a genetic variance that should be approached regardless of trait or environmental dynamics
- Also some branching theory work to show extinction risk

##### Results
- For one generation changes, heritability doesn't effect population size, but phenotypic variation does
	- high variation has smaller population size for low-severity event (standing load) but higher size for high-severity event
	- "High phenotypic variance therefore reduces both mean fitness within a generation and the variance in fitness across generations - a form of short-term bet-hedging [...] The negative impact of load is more apparent in the case of a mild event, when the population remains near the fitness peak and variation means that more individuals are spread further from the peak."
- Heritability influences recovery: heritability dampens growth rate in subsequent generations because high heritability means maladaptation when the environment returns to normal
- Higher phenotypic variation means lower extinction risk (more individuals survive an extreme event)
- For fixed severity, Heritability lowers extinction risk as the duration of events increases
- Calculating lag load and its contributions oer time: cumulative long-term lag load expression a function of "evolvability", squared severty of event
	- as duration increases, lag load as a function of heritability goes from positive (duration low) to negative (duration high)

Very interesting paper, thankfully quite different from ours. Idea here is that *if the environment returns to its pre-disturbance state* on a short timescale, heritability leads to maladaptation. This is mentioned in the discussion, but this is a question of environmental autocorrelation. The classic rescue model is essentially very high autocorrelation, this in a sense is kind of low (at least in one direction - a ratio of transition probs is high/low). This paper demonstrates how helpful heritability is in these. Thinking about it for two seconds, heritability will of course be helpful overall if the autocorrelation his high!

### Stewart, G.S., et al. 2017. The power of evolutionary rescue is constrained by genetic load. Evolutionary Applications.

##### Framing
- Evolutionary rescue in diploid sexual reproduction understudied
- Recessive and deleterious mutations can be masked by dominance but can be revealed at low population size
- Genetic load: population mean deviance from an optimal phenotoype
	- Increased as drift and inbreeding increase
- Genetic load = drift load + inbreeding load
	- Drift load: fitness of outcross (btwn pop) offspring vs. fitness of within-population randomly mated offspring
	- Inbreeding load: fitness of mating between relatives vs. random mating within pop
	- high outcross population is heterosis
- Founding populations in novel environments: evolving but also possibly could accrue genetic load

##### Experiment
- Beetles in to either a favorable or unfavorable environment
- Initial population size: 10 or 40
- Evolving vs. non-evolving: one-for-one replacement in non-evolving
- After generation 8, mix populations to assess heterosis
	- Treatments: mixed evolving, unmixed evolving, non-evolving control
- Evaluate: 
	- Did Pr(extinction) vary by treatments
	- Time to extinction
	- Growth rate
- Founding size and evolution treatment had similar effects on extinction rates
	- n.s. evolution x founding size term
- In poor environment, small populations went extinct sooner
	- in small poulations, evolution was 1.2x faster in evolving pops
		- (likely reached homozygosity quickly)
	- in large populations, evolution was 1.2 times slower
- In favorable environment, evolving populatons grew but then slightly shrank while non-evolving was more consistent in size
	- growth rates in favorable: evolvers started high but then declined
	- in all environments, population growth rates declined in final few time steps
- outbreeding populations have higher growth rates rel. to close populations in both environments (also better than controls)

##### Discussion
- Even small founding poopulations in poor environments had reduced extinction risk if they could evolve
- The fact control populations could grow in favorable environment means that habitat quality is important
- "Strong selection on traits conferring adaptation could have reduced overall genetic diversity and increased homozygosty in our experimental populations, and at the same time, genetic drift due to small population size could also have contributed to the rapid fixation of deleterious alleles." on why there is genetic load
- Mixture even in the poor environment confers adaptation

Here: as the title says, rescue is constrained by genetic load. Accumulation of deleterious mutations due to small population sizes (n.b. comparable to our populations) likely candidate.

### Jain, K., and W. Stephan. 2017. Rapid adaptation of a polygenic trait after a sudden environmental shift. Genetics. 

##### Model
- Polygenic trait with sudden shift in optimum: diallelic loci ($l$) of unqual effects
	- see also Jain and Stephan 2015, Wright 1935, Barton 1986, de Vladar and Barton, 2014
- Infinitely large population of diploids
- positve allele has (g_i)/2 at freq. p_i and negative is (-g_i)/2
- trait c1 is the sum across all l loci
	- n.b. phenotypes are not necessarily gaussian
- quadratic selection: w(z) = 1 - (s/2)(z-z')^2
	- somehow produces average fitness which is gaussian-looking
- analytic expression for dp_i/dt
	- one term to stabilize the phenotypic mean to optimum value
	- one term to produce fixation of (one) allele
	- final term changes between pos/neg alleles (mutations)
- equilibrium frequencies of p_i are roots of a cubic polynomial
	- three roots when mean (phenotype?) is optimum
		- of these, p = 1/2 is stable if effect is strongr than a threshold: 2 * sqrt(2 mu / s)
		- if selection is weaker than mutation, eqm. freq. is 1/2
- when most effects are small, equilibrated to a pheno of 0, p approx 1/2
	- then, dp/dt = S_i p_i q_i for S_i = -s g_i (c1 - zf)
- closed for expr. for trait mean over time

(resume on p. 392)

### Peniston, J.H., et al. 2020. Environmental fluctuations can promote evolutionary rescue in high-extinction scenarios. Proc. R. Soc. B.

##### Previous stuff on environmental variance (temporal)
- Environmental fluctuations understudied rel. to magn. of environmental change
- In general, temporal autocorrelation (typically positive) can increase or decrease extinction probabilities (refs. herein)
- Prior work suggests that greater magnitude of variation -> decreased rescue probability (decreasing expected population size and increased variation in pop. size)
	- But, autocorrelation: pop. size makes it easier to track change sin optimal phenotype so more (pos.) autocorrelation is good, although it increases variance and can cause a run of bd years.

##### Models
- *Baseline risk* is probability of a population going extinct following abrupt environmental change (w/o fluctuations)
	- a function of initial size, mutation rate, rate of environmental change
###### *Polygenic model*
- Quantitative trait ($z$) determining survival prob., with Gaussian selection
	- additive genetic component $g$ w/ normally distributed noise
	- amount of genotypic and phenotypic variation assumed constant
	- no density dependence
- Phenotyipe changes through a gaussian random sequence with autocorrelation $\rho$, and mean $\mu$, s.d. $\theta$.
- Population dynamics smulated over 1000 generations
- Also simulated an IBM with ten freely recombining additive loci
- G&H gives a baseline extinction risk w/o variability
	- if populaton will go below threshold (deterministically) high ext. risk
###### *Monogenic model*
- Adaptation only through novel mutations (no reliance on standing variation)
- Mutation: low probability (u) for each newborn, one possible mutation, no back-mutation
	- Mutations: probability of adult dying before adulthood
- Number of births evolving over time in relation to gaussian random sequence
- Here: rescue req's a single mutant to produce a lineage that persists
	- Here: branching process

##### Results
###### Polygenic model
- mild env. changes (baseline risk low): variation in optimal pheno means more extinction due to runs of bad years
- harsh env. changes (baseline risk high): variation in opt. pheno means less extinction due to runs of good years
- The first change was very important: mean of opt. phenotype predicted outcome with 96% accuracy
###### Monogenic model
- Variation in fecundity increases mean cumulative # births before extinction
	- More mutations
- But, variation in fecundity does not affect prob. of lineage surviving once arising
	- Survival/extinction determined by first few generations after mutation
###### Overall
- Following abrupt change, the effects of variation depend on baseline extinction risk
	- Low-extinction cases see more ext. risk w/ variation
	- High-extinction cases: extinction likely anyway, but variation means runs of good years favorable to adaptation
- Fluctuating phenotypic optimum vs. fluctuating fecundity - similar results
	- variation in pop'n trajectories, not necessarily the lag load, likely are causing changes in extinction
- "our models should be a reasonable approximation of evolutionary rescue dynamcs if there is weak density dependence, but given strong density dependence that acts even at small population sizes the rescue process might differ. In addition, in models without density dependence, once  population is locally adpted (e.g. mean bsolute fitness is greater than 1), the population size can grow without bound. Therefore, most of our results only pertain to extinction over relatiely shorttime-scales due to maladaptation following environmental change and do not pertain to the long-term extinction risk in a population following adaptation. Following adaptation, environmental variation should increase the probability of extinction in a finite population, and the expectd time to extinction is known to scale with carryng capacity"
- Managers: for populations whch are *a priori* at low risk, provide stable environment, but for high-risk populations intervene as soon as possible
- "Effects of environmental variation that we have explored in this study are relatively minor compared with the importance of the magnitude of envronmental change or initial population size."

Environmental assists populations likely to fail (good years), hurts populations likelyt o to survive (bad years). The initial trajectories of populations are very important. (In context of density dependence: this will depress growth rates, making populations inherently higher risk, meaning good years more likely good?)

### Lande, R. 1993. Risks of population extinction from demographic and environmental stochasticity and random catastrophes. American Naturalist.

- Stochastic diffusion approximations to examine expected time to extinction in populations under varying growth rates and levels of stochasticty
- Model: two-piece model with linear/expo growth up to carrying capacity, then stuck there
	- this is the Macarthur and Wilson model (1967)
	- (n.b., our model behaves differently)
- With demographic stochasticity, variance in growth rate is just V1/N (individual variance divided by # individuals)
- Demographic stochasticity: time to extinction is approx. exponential if r > 0, linear if r = 0, sublinear for r < 1
- In long run for sufficently large sze mean extinction time is nearly independent of initial size because population sizes are expected to quickly grow to CC and fluctuate before stochasticity can have large effect

Hmm... not sure how to relate this to ours. Growth form is different (growth rate here is always going to be greater than ours, so upper bound), our variance in growth rate is also slightly higher. Still, we have a carrying capacity that increases - rate of increase in CC has increasingly large effect due to exponential scaling.

### Gomulkiewicz, R., and D. Houle. 2009. Demographic and genetic constraints on evolution. The American Naturalist.

##### Evolutionary constraints
- Evolutionary constraints prevent a population from obtaining a target phenotype
	- **Absolute constraint** on evolution: population lacks genetic variation in some direction of phenotypic space
		- e.g., absence of variation in particular combos of traits
		- produces a singular $G$ matrix of additive genetic variances and covariances
	- **Genetic constraint** - ability to respond to selection is present but limited (??)
- Demographic limitations make quantitative constraints impossible to overcome

(did not fully read - very opaque)

### Whitlock, M.C. 2000. Fixation of new alleles and the extinction of small populations: drift load, benefcal alleles, and sexual selection. Evolution.

- Prob. of fixation of a new mutation: function of its effect and effective population size
- Change in fitness depends a lot on the distribution of allelic effects
- Assume for rest of paper that distribution is exponential
- ISO: critical Ne below which fixation of deleterious alleles overwhelm fixation of beneficial alleles
	- delt W = delt W_b + delt W_d
	- It exists (eq. 9) - proportional to cubic root of ratio of mutation rates of deleterious to beneficial alleles, inv. proportional to mean effect of beneficial and deleterious alleles (product)
	- Higher when:
		- Rate of beneficial mutation is high rel. to deleterious
		- mean effects of mutation (both?) are high
		- Ne is high (compared to N)
		- sexual selection removes deleterious alleles
		- beneficial alleles are more common or have stronger effect when compensating for decline (when does this happen?)
- Below the threshold, loss of fitness becomes rapid with decreasing Ne, but above the threshold, gain in fitness is slow with increasing Ne
	- "short periods at small size will require that the population spend many times longer at large size to compensate for the fitness loss"
- Ne in the hundreds is typically good for indefinite persistence w/o genetic issues
	- with thousands of individuals, prob. of deleterious alleles fixing faster than beneficial ones is low
	- more likelihood of extinction due to catastrophes
- For small Ne, fixation of benefcial alleles can still stem mutation loss
	- but below the critical value, these beneficial fixations have negligible effects

This paper is looking at mutations but it's unclear if all aspects of it relate to mutations. Two interesting ideas here: being very small has disproportionately large effects (comp. to being large) for change in fitness. Also, Ne in the hundreds is good for avoiding potential fitness losses due to genetics. Thinking about how NDD constrains population size here is good and important.

### Frankham, R. 2005. Genetics and extinction. Biological Conservation. 

##### Background
- Inbreeding: H_t / H_0 = (1 - 1/(2Ne))^t = 1 - F (H_t is heterozygosty at gen t)
	- middle term approx. e^(-t/(2Ne)) - exponential decay inv. proportional to population size
- Controversies: does inbreeding impact extinciton risk?
	- Does selection purge deleterious alleles, reducing their impact?
	- Does inbreeding affect wild populations?
	- Do mutations influence extinctions in sexually reproducing species?

##### Inbreeding and extinction
- Inbreeding having deleterious effects (fitness, reproduction, survival) pretty widely demonstrated
- Newman and Pilson, 1997: less inbred populations (F = 4%) had 25% extinction rate but slightly more inbred (8-9%) had ~69% - 26% of variation in extinction risk was attributable to inbreeding
- Simulations of populations with genetic and demographic extinction threats: inbreeding decreases time to extinction
	- esp. in outbreeding taxa
	- pop'ns with already low growth rates were less affected (already headed to extinction)
- If inbreeding/genetic threats do not induce extinctions (see Lande 1988), then genetic variation will not vary by extinction threat levels
	- but Spielman et al., 2004, found 77% of 170 taxa surveyed had reduced genetic diversity (median difference in heterozygosity was 40%)
	- N_e / N is typically around 0.1 rather than 0.25 - 0.5 as previously thouhgt (see Frankham, 1995)
	- Lande 1998 also is overly-optimistic about the degree of purging; see refs herein
- Outbreeding typically removes inbreeding depression
- *Mutational meltdown* - accumulation of negative alleles reducing growth rates and eventually inducing extinction

Good general reference for extinction risks in small populations. Mutational meltdown is an intersting idea outlined here (ths is likely nt the first place). Also some reviews of empirical things in here. 	
### Lynch, M., et al. 1995. Mutational meltdowns in sexual populations. Evolution.

##### Background
- Muller's ratchet (not explained here): accumulaton of deleterious alleles
	- W/o recombination, offspring can't have fewer deleterious mutations than offspring
	- Accumulation of deleterious alleles -> eventual extinction
	- Segregation and recombination lower the probability of fixation
	- Acc'd to Kimura et al. 1963, for high Ne, segregation of deleterious mutations reduces mean fitness by factor of at most exp(-m) for mutation rate m per individual across all loci
		- i.e., large sexual populations that can produce enough offspring to offset this difference will not face extinction risk due to deleterious mutation
- Charlesworth 1993 and Lande 1994 address how much mutation accumulation should affect extinction times in large populations
	- **Hill-Robertson effect** background genetic variation reduces efficiency of slection on indivdual loci
	- HR-effect is a violation of some of the assumptions from Lande 1994
	- Also: Lande 1994 ignores changing population size, reducing arithmetic mean of the Ne
- What is the risk that large sexually reproducing populations face from accumulating deleterious mutations?
- Lynch et al. 1993: mutation accumulaton has three phases
	- (as long as population is initially mutation free and reproductive rate is high)
	- Phase 1: mutations accumulate rapidly because selection is low (due to low genetic variance)
	- Phase 2: Mutation, drift, selection start to balance,
	- Phase 3: WR = 1 (i.e., replacement rate?) increases in the mutation load will reduce population size, increasing the rate of mutaton accumulation as drift overwhelms selection
- 

##### Model
- Simulations similar to Charlesworth 1993
- Monoecious randomly mating populations with discrete generations
	- R is reproductive rate per (surviving) adult
	- N(t) is number of reproducing adults at time t
- Viability selection on number of deleterious mutations
	- W(n1,n2) = 1 for no mutations, W < 1 for mutations
	- Mutations are poisson distributed with rate \mu
- Analytical approximations with diffusion theory
	- assumes: constant population size, gives length of phase 2
	- mean viablity is product of segregting mutations (fluctuating) and fixations (increasing)
- Transition-probability theory approach (see Ewens 1979) 
	- assume constant population size $K$
	- probability of going from i to j loci in a time step is Binomal with size 2K and prob nu

##### Results and discussion
- As found in other places, mean time to extinction scales with log(R)
	- Populations have "steady fractional decline in fitness" during phase 2 but phase 3 doesn't begin until mean viability is 1/R
- Whatever.
- "It is now well established that mutations with small, but intermediate, deleterious effects cause the most cumulative damage to populatons."
	- Lande showed that generally the selection pressure that minimizs time to extinction is inversely proportional to the effective population size
- "Roughly speaking, immunity to the meltdown requires a large enough population size that during the period in which we predict a high probability of extinction (due purely to unconditionally dleterious mutations), there are enough opportunities (individuals x generations) with effectve enough selection to correct the damage.

Jesus, I don't know. Main finding is you need large populations to avoid extinction by the accumulation of mutations. This scales with the log of mean fitness. The idea here of loss of fitness through segregation (?) and fixation is worthwhile.

### Keller, L. and D.M. Waller. 2002. Inbreeding effects in wild populations. TREE.

- Small isolated populations face two threats: fixation/loss of alleles eroding genetic variation (drift) and accumulation of deleterious alleles
	- Leads to mutatonal meltdowns although this happens gradually
	- Inbreeding happens more quickly
- Inbreeding - shared ancestry
	- Wright calculated through correlation, but now calculated by probabilities
	- $F_i$ is prob of individual having two genes being IBD (identical by descent)
	- Inbreeding also measured by F_{is} = 1 - H_o / H_e (heterozygosity observed/expected under HWE) where F_{is} > 0 suggests more inbreeding than expected; this gives deviation from random mating
- Inbreeding would be less prevalent in natural populations if it was purged (inbred individiauls don't survive/reproduce)
	- Purging most effective when
		- Deleterious mutational strength strong rel. to Ne
		- Selection does not have strong cross-loci effects
		- Gradual inbreeding
		- Isolation means purged alleles do not reappear
	- Meta-analysis showed tht purging occurred in ~33% of studies observable in ~10%
- Inbreeding depression - typically within populations measured as slope between inbreeding coefficient of individuls (F) and log fitness
- "Whereas our understanding of how inbreeding affects individuals has improved, we still know little about the effects of inbreeding on population dynamics. In the presence of strong density-dependent regulation, inbreeding depression might eliminate inbred individuals without affecting population dynamics (so-called "soft selecton"). In such circumstances, the magnitude of inbreeding depression would not necessarily directly affect population dnamcs. When selection is 'hard', however, inbreeding could substantially affect population dynamics."
- "It would lso be of interest to know how variance in the degree of inbreeding among individuals within a population affects the expression of inbreeding depression and subsequent population dynamcs."

Good things in here. Drift vs. genetic load and their relationship to inbreding. Measuring inbreeding (may be annoyingly difficult in our populations). Reviews suggesting that inbreeding does in fact exist and influencs individuals, although there is an opening for looking at population dynamics (esp. with NDD - NDD and hard selection should make a difference!)

### Comeron, J.M., Williford, A., and R.M. Kliman. 2008. The Hill-Robertson effect: evolutionary consequences of weak selection and linkage in finite populations. Heredity, 100: 19 - 31.

- Fisher (1930) says that fixation of beneficial mutations depends on presence of other mutations, the presence of which could slow down fixation of beneficial mutation
- Hill and Robertson (1966) - two loci in a finite population (genetic drift), selection at one locus interferes with selection at the other
	- in general multiple loci meas less efficient selection	- degree of interference increases with genetic linkage, such that genomic regions wth reduced recombination will have reduced Ne
- Recombination: Felsenstein (1974) showed that recombination increases the rate of adaptation by reducing linkage disequilibrium 
	- (n.b. apparently this is not the primary feature of HR effect but instead a consequence)
- Other modls of interference between segregating mutations: Muller's ratchet (stochastic accumulation of deleterious mutations through drift and linkage), hitchiking (selective sweeps also bringing bad mutations to fixation), background selection (selection aginst deleterious mutations also removing linked neutral mutations)
- Genomes likely contain a lot of weakly selected sites, which are likely to be clustered in exons and regulatory regions (which would amplify the HR effect)

Not sure... makes more sense to go and read original Hill and Robertson. This review looks mainly at more complicated genomes than we are modeling.

### Orr, H.A., and R.L. Unckless. 2008. Population extincton and the genetics of adaptation. Am. Nat.

##### Framing
- "Nearly all work on the theory of adaptation, however, has been limited in at least one important way: it assumes that populations maintan a stable size following a sudden change in environment, a change to which the population must adapt."
- This study: looking at rescue from new mutations and/or standing variation
	- mutation: "we calculate the probablity that a population substitutes a beneficial allele of sufficently large effect to escape extinction"
- "our results are most directly relevant to organisms such as microbes, in which a small number of genes can respond to n environmental challenge"

##### Model
###### Rescue from mutation
- Population size changes geometrically up until carrying capacity $K$ (the Macarthur and Wilson model)
- Haploid genetics, one locus, two alleles $A$ and initially absent/rare mutant $a$
	- wt has W_A = 1 - r, mutant has W_a = 1 + s - r (s>r)
- Use theory of branching processes
	- so long as a lineage with W_a > 1 becomes large enough it will survive, so we are most interested in probabillity of ultimate survival
- Prob. of arising in generation $t$ and going to fixation: 1 - probabilty that all mutations are lost accidentally
	- under assumptions this is 1 - exp(- (2 N0 u(s-r)) / r)
	- Prob(T = t | surv) approx r exp(-rt) where T is generation where the mutation occurs
###### Rescue from standing variation
- e.g., alleles previously deleterious reaching fixation?
	- if this allele is favored for t > 0 after environmental change, al it has to do is escape accidental loss
	- avoid extinction with approx prob. 1 - exp(- (2N0u(s-r))/sd) for selectve advantage sd
- Mean probablty of survival is 1 - ((s-r+sd)/sd)^(-2N0u)
###### Overall probablty of rescue
- Overall probability of rescue is Pr(standing) + (1-Pr(standing) * Pr(mutation)

###### Discussion
- "our qualitative conclusion agrees with that of Gomulkiewicz and Holt; unless populations are very large or suffer only modest decreases in absolute fitness below 1, adaptation to a sudden environmental change can be difficult"
- As r becomes small, new mutations are more likely to rescue a population because decline is slow
- Evolution from new mutations or from mutation-selection balance are assisted by: s large (large selective advantage for mutant), N0u is large (number of mutations? u not defined), r small (slow population decline), sd small (small initial disadvantage for deleterious)

Use for some similarities to G&H - notably that initial population size is important. Also talk about population size and mutation.

### Orr, H.A. 2005. The genetic theory of adaptation: a brief history. Nature Reviews: Genetics.

##### The "micromutation" paradigm
- "Micromutation" small unspecified mutations of small effect
- The quantitative trait model - from Fisher, regarded today as mostly mathematical convenience
- Fisher's "geometric model" - alternative to the quantitative trait model, features random mutations of unequal size in n-dimensional phentoypic space, pleiotropy, etc.
	- Kimura pointed out in the 1990s that under this model only mutations of intermedate effect are likely to contribute to adaptation
- Quantitative trait loci (QTLs, mid 1980s) - use a large suite of molecular markers (mapped) and associate them with population/species characteristics

##### Sequence-based techniques
- John Maynard Smith: real adaptation occurs in discrete steps, not a continuous phenotypic space (more reflective of molecular mechanism)
	- for a L-length base sequence, there are really only approx 3L possible mutations
	- populations make "adaptive walks" through this space
- The "NK model" features fitness landscapes of varying shape/ruggedness, and probability of arriving at global optima bsed on initial condition, etc.
	- but, rules for moving through the landscape are unrealistic, and this model typically ignores starts from highly fit wild-type phenotypes
- John Gillespie: the distribution of fitnesses is often unknown, but WTs are probably on the right tail and adaptations will be even more extreme
	- use extreme value theorem, which makes few if any assumptions about the shape of the distribution, for this
	- when presented with several mutations, populations will follow mutations proportional to their selective advantage
	- number of steps needed to reach new optimum is very small (2-5)
	- subsequent work demonstrting that adaptation features large jumps in fitness rank (Pareto principle - most of the adaptive benefit can come from few substitutions)
- (thoughts on current hurdles/future challenges/directions on final page)

### Bell, G. 2012. Evolutionary rescue of a green alga kept in the dark. Biol Lett.

##### Experiment
- Prediction: rescue should be most effective in outbred populations, somewhat effective in inbred sexual populations, least effective in asexual populations (i.e., gradient of genetic diversity)
- Algae populations, all popagated for thousands of gens. across many years:
	- Obligate large sexual populations (three replicates)
	- Bottlenecked obligate sexual populations (two lines) propagated each gen. from one individual
	- Mixed mating populatons
	- Asexual populations
- Two experiments involving plates growing in the dark:
	- Long term (50 growth cycles) and short term (20 cycles)
- In long-term experiment, most populations went extinct within 10 cycles
	- surviving populations declined at first but recovered with U-shaped curve
- Short-term sexual lines, survival varied by treatment
	- Outbred > bottlenecked obligate > mixed > asexual
	- diff. between asexual and mixed is n.s.

Sex is important, and within the realm of sexually reproducing populations genetic variation is important. Consider this for contrasting sex vs. asexual and also effects of bottlenecking. 

### Pray, L.A., et al. 1996. The effect of population size on effective population size: an empiritcal study in the red flour beetle *Tribolium castaneum*. Genet. Res. Comb.

##### Background
- Nunney (1993) speculated that Ne/N > 0.5 in most populatons, including small ones (some empirical stuff cited here)
- Frankham (1995) argues that Ne.N approx 0.1 is more likely common

##### Study
- Studied Ne in replicates of populations initiated from size 2, 4, ..., 48, 480, 960 (9 - 25 reps/each)
	- each populaton was initiated with one cSM ++ genotype and remaining cSM bb (geno corresponds to color)
	- Populations run for only one generation? Counting each geno
	- among offspring generations, mean and varance of +b indviduals calculated
- Some statistics and corrections used to estimate Ne
- Ne/N appears to shrink logarithmically with N
	- is approx 1 for N = 2, approx 0.8 for N = 960
	- "a decrease in Ne/N is proportional to an increase in Va"
	- "there is a significant increase in V_k with logN"

##### Discussion
- Decrease along large N-gradient is consistent with Husband and Barrett (1992) as well as recalculations from Nunney (1995). 

### Lande, R., and Shannon, S. 1996. The role of genetic variation in adaptaton and population persistence in a changing environment. Evolution, 50(1): 434 - 437. 

- "In the short term, genetic variability is often less critical than other determinants of population persistence (Lande 1988), but in the long-term, it can play the decisive role in allowing a population to persist and adapt in a changing environment."
	- "However, gnetic variability is thought not to be the rate-limiting factor in long-term evolution... instead, long-term rates of evolution and adaptive radation are constrained by ecological opportunity [refs]"
- Model: rate if increase, r, is the sum
	r_max - (maladaptaton penalty) - density dependent term [(n.b. looks very similar to Lande 76/G&H)]
	- Want: time-averaged fitness (rbar) > 1
	- Note: the maladaptation penalty includes standing load and maladaptation load
	- The relative sizes of maladaptation and standing load depends on the type of selection
- In a constant environment, any amount of phenotypic variance will increase standing load and decrease growth rates
- Under a directional change in fitness (deterministic), mean phenotype tracks optimum with a predictableevolutionary load
	- the evolutionary load is more responsive to genotypic variation than to the rate of adaptation, so more variaton is good (closer tracking)
- In a cyclically changing environment, more variance is good for slow, small-amplitude oscillations
- Table gives evolutionary loads for randomly fluctuating and autocorrelated random environments
	- In totally unpredictable envronment, additive genetic variation increases expected deviation from optimal environment
	- In an autocorrelated environment, load decreases with additive variance as well as autocorrelation time
- So, genetic variance is not always adaptive. t is maladaptive in cases where increasing adaptve capacity will lead it into evolutionary traps
	- autocorrelation is helpful here for ensuring that adaptation doesn't end in traps

Interesting but not necessarily relevant to us due to timescales and type of environmental change (although will give some intuition for more complex models).

### Arnaud, L., Haubruge, E., and Gage, M.J.G. 2001. Morphology of Tribolium castaneum male genitalia and its possible role in sperm competition and cryptic female choice. Belgian Journal of Zoology, 132(2): 111 - 115.

- Sexual selection post-copulation: when mating with multiple males females can select for certain sperm
	- e.g. preferental uptake, directional movement into different storage organs
- In *T. castaneum* the refs. suggesting that the last male to mate achives predcence (fertilizes 2/3 of ova)
- Here, looking at male genetalia in reared *T. castaneum*
	- They have an aedeagus, a scoop with hard spines, which trap and remove rival sperm
	- Aedeagal length is correlated with body size!
		- other studies have shown that rel. body size influences success in sperm competition
	- The bristled hairs on the retractable brush may dislodge rival sperm
		- or it may be used to signal to females

So: the aedeagus (which may vary in size among individuals) is an organ that appears to have a scoop and brushes which can be used to remove rival sperm to increase fertilization potential for the most recent mate.

### Pai, A., and Yan, G. 2020. Long-term study of female multple mating indicates direct benefits in *Triboleim castaneum*. Entomologia Experimentalis et Applicata, doi:10.1111/eea.12906

- Polyandry (multiple-mating) - why could or would it be advantageous for females? 
	- E.g., Multiple mating may allow sperm competition
	- Protects from infertility in case first mate has bad sperm
	- The inter-mating interval is important: short time periods allows sperm competition
- Three hypothses for why periodic multiple mating may be good
	- Fresh sperm means more offspring production
		- either because old sperm is bad or because there's less of it
		- females can have long mating periods (several months) although fertility typically is small after the first month
		- Test: females who can mate at least once a month will have more offspring than less often
	- Ejaculate contents (besides sperm?) may have effects
		- e.g., what if it stimulates more egg production
		- Test: ehanced survivorship and/or offspring for females with more mating opportunites
	- Genetic benefits
		- genetic mate compatbility? lots of matings are unsuccessful
		- polyandry may increase odds of receiving compatible sperm
		- Test: polyandrous offspring production higher than singly-mated female mating at similar rates
- Experiment on *T. castaneum*
	- Nine treatments
		- Mating period (2wk, month, 1/3 months, 1/5 months)
		- Distinct mates (distinct vs. multiple)
	- Recording: n. mating opportunities, n. mates, n. offspring, lifespan
	- ANOVA: sqrt(offspring) as respons
	- Sperm replenishment test: paired t-test for offspring prod.
	- Frequency test: ANOVA comparing offspr ~ mating freq.
	- Ejaculate effect test: product limit survival analysis
- Results
	- Everyone produced offspring
		- so n. mates does not affect offspring yes/no
	- Remating frequency (but not n. distinct partners?) influenced n. offspring
		- More frequent mating, but not more mates per se, increases n. offspring
	- More offspring in younger mothers, and higher n. offspring after mating when mating every 3-5 months but not >=1/month
		- I.e., females produce more offspring if second mating is ~3 months after initial mating
	- Survival was not affected by mating
	- Fig. 2 is good: n. offspring tails off after mating bout in all groups, but increases after extra mating, such that but mating even more frequently has big impact on n. offspring
- Sperm replenishment: more frequent mating is good, because sperm apparnetly doens't stick around for very long (offspring tails off after mating)
- Female longevity: no evidence supporting this in *T. castaneum* although it does in some species
- Genetic compatibility: no evidence supporting this
	- although n.b. it would be helpful for avoiding inbreeding

In the context of our rescue study: Triboleum can mate wth multiple males (won't affect number of offspring). They can mate multple times over the course of their lifetimes.

### Chevin, L, and Lande, R. 2010. When do adaptive plasticity and genetic evolution prevent extinction of a density-regulated population? Evolution, 64(4): 1143 - 1150.

- Single environmental shift in an otherwise-stable environment
- Phenotypic plasticity (see Hendry et al. 2008) which is itself evolvable
- Similar to Lande 2009
	- Lande '09 shows that w/ sufficient variation in plasticity, first step in rescue is rapid increase in plastcty then slow genetic change supplanting it
- trait $z$ with a linear reaction norm, $z = a + beps + e$
	- eps is environmental variable
	- e residual variaton
	- a is breeding value in (arb.) reference environment
	- b is plasticity (slope, how much env. change affects pheno)
	- sig_e, sig_a, sig_b are genetic variances assoc. with each var s/t sig^2_z = sig^2_e + sig^2_a + sig^2-b
		- assume sig_a, sig_b constant over time
- Fitness of individual with trait z, plasticity b, in pop. size N is
	- R(z,b,N) = R0^{1 - f(N)/f(K0)} W(z) C(b)
	- K0 carrying capacity in abs of plasticity and genetic load
	- W(z) is gaussan selection, exp(-(z-Beps)^2 / (2 om_z^2))
		- so fitness changes linearly with environment acc'd to B
	- C(b) is cost of plasticity, exp(-b^2 / (w \om_b^2))
- Selection first acts on b, then on z
- Phase 1 of evolutionary process: rapid increase "mean plasticity" (mean of b)
- Phase 2 of process: loss of plasticity compensated for by trait evoluton (increase in mean of height in reaction norms, a?)
	- trait (z) not affected here
- Density dependence: f(N) = N is Ricker model, f(N) = N^{th} is theta-logistic
	- "By comparison [to GH95], density dependence accelerates the declining phase and slows the increasing phase, reducng population size at all phases"
- More plasticity means less time to reach phenotypic optimum and higher population trough
- Persistence: more plasticity means the initial population size required for persistence can be smaller (rel. to extinction threshold)
	- for alpha > 0.5 persistence seems guaranteed?

Main effects of NDD studied here are simply the obvious ones: steeper decline, slower ascent. Because plasticity is faster and then is replaced at rate by adaptation, more plasticity means more rapid adaptation and lower trough. Without plasticity the initial size needed to avert extinction (under NDD) is much larger.

### Reznick, DN and Ghalambor, CK. 2001. The population ecology of contemporary adaptations: what empirical studies reveal about the conditions that promote adaptive evolution. Genetica, 112-113: 183 - 198.

- "Recent reviews suggest that natural populations are capable of sustained rates of evolution many orders of magnitude higher than those observed in the fossil record [refs]"
- Here: review of studies focusng on rapid adaptive evolution, phenotypic change, genetic basis for observed trait change
##### Studies
- Most empirical studies look at resp. to anthropogenic changes in environment (e.g., heavy metal tolerance)
- 47 studies in all
	- 36 with some genetic basis for change and 11 with only evidence of pheno change
	- 38 with rapid evolution in animals (16 vertebrate species, 19 invert)
	- Nine terrestral isects, five aquatic inverts
- Traits under selection: morphology, phys, lfe history, phen, behavior
- "almost all cases of contemporary adaptation were associated with 'colonization' events." (i.e. into new habitats)
	- all but six studies had colonizaton components
	- coloninies could be isolated, or could result from heterogeneity in environments allowing adaptation (e.g., metapopulations), i.e., is there gene flow possible between populatons
- Adaptations
	- use of novel resource/food source (e.g., soap bugs feeding on novel plants)
	- bio-physical changes (e.g., metal tolerance)
	- changes to mortality (interesting life history), mostly predation
	- character dsplacement? competition?
##### Discussion
- Adaptations are typically polygenic and involve a suite of adaptations (e.g., life history of guppies, size at maturity, devo rate, fecundity, offspring size)
- [...]

So... uh, radid adaptation can occur, although usually in several traits.

### Shaw, RG, and Etterson, JR. 2012. Rapid climate change and the rate of adaptation: insight from experimental quantitative genetics. New Phytologist 199: 752 - 765.

- "a critcal outstanding queston of our time is the extent to which adaptve evoluton will ameliorate the detrimental effects of rapid contemporary climate change"
- Migration: populations could "track" climate by moving polewards or up elevational gradients
	- meta-analyses suggest range contraction in 60-90% of species in USDA Forest Inventory and only ~20% of these sp. show poleward expansion
	- Habitat fragmentation makes dispersal polewards more difficult
- Adaptaton and rescue
	- grass species and metal tolerance (Antonovics refs)
		- although many species by mines dd not adapt
	- widespread angiosperm resistance to herbicide
- Fitness links demography with evolutionary processes
	- What influences the rate of adaptation? 
		- Mutation, recombnation.
		- Drft (random sampling) opposes by reducing VA(w)
	- Lopez et al. 2009 on fragmentation and reduced population size and its effects on drift
	- Gene flow: good by alleviating inbreeding and potentially bringing in useful genetic diversity (esp. for similar environments), bad if incoming genes are maladapted esp. in temporally stable environments
	- See Pease et al. 1989, Polechova et al. 2009 for adaptation + migration + NDD
- Experiments and predictions
	- Traits which are/are not likely to be under selection pressure suggestng populations may have considerable genetic variation with heritability although some trats may have negative correlation (e.g., water-use efficiency and rosette size) (maybe use these refs?)
- "To our knowledge, [VA(w) and mean absolute fitness W] have never been jointly estimated for any population in nature, perhaps largely because of the challenges of doing so"
	- but you can use Aster mods for measuring mean lifetime fitness

So: some stuff in here about how much genetic variation exists in wild populations, but also we lack joint estimation of these and don't have a thorough knowledge of which traits are actually under selection and/or what their variance+heritabilities are. This is all for plants. Good for overviews.

### Godwin, JL, et al. 2020. Mating patterns influence vulnerability to the extinction vortex. Global Change Biology.

- "the extinction vortex has proven to be an informative model for conceptualizing [...]" refs
- "[...] ample evidence across a range of systems showing that individual or combined effects of envronmental stress and genetcs lead to reduced population size or even extinction" refs
	- but, no explicit experimental tests of the vort model
- This study: explicit experimental test of the vortex model, where resp. is "consquences for population viability"
	- Trt 1: dietary restriction (restrain missing protein)
	- Trt 2: environmental stress (heat wave)
	- Trt 3: bottleneck
	- Resp: # offspring produced after ~15 generations
	- Also, trt 4? Mating patterns (monogamy vs. polyandry)
		- sexual selection can create more btwn-individual variance in reproductive success potentially accentuating effects of vortex (due to adaptation loss?)
- "Genic capture theory" - reproductive success in face of competton and choice depends on "wider conditon which is contingent on genome-wide variation" (??)
	- sexual selection can reduce load by (1) making mutations more deleterous in males and (2) load ffects mating and fertilization s/t average father has less fitness load than average male (hmm...)
	- under this theory (?) sexual selection may enhance population fate in the face of vortex by increasing overall condition (??)
	- or, if sex selection operates in opposition to natural selection (e.g., antlers too big) or by creating different reproductive optima for males and females
##### Experiment
- Mating treatment: enforced monogamy vs. polyandrous (five males to one female) 
	- three independent lines within each regme with 20 monogmous or 12 polyandrous groups per line (Ne approx 40)
- Vorte: 95 generations of experimentl evolution wth strong vs. weak selection
	- 15 generations inside the vortex
	- in each line, nine replicate populations initiated with non-sib virgin pairs (9 pop'ns x 3 lines x 2 sexual selection treatments = 54 populations)
	- populations initiated... somehow
	- in each treatment: generation initiated wth 5m + 5f, after mating 10 pupae sexed and used for initiating next gen while remaining pupae are usd to calculate population fitness
	- then run through five three-generation cycles of stress (nutrents then thermal then or genetic)
		- nutrients: replicate small popn's of 5m and 5f in single sex groups w/ nutrient-poor food for 10 days then placed together for 7; some pupae taken for next generation and remainder counted for pop'n fitness
		- thermal: small populations (Ne = 10) exposed to 5 day heat wave during mating stage; heat will damage male fertility
		- bottleneck: inbreeding depresson based on restricting population to size Ne = 2
##### Analysis
- Responses: 
	- time to extinction (accelerated failure time model; ext. time ~ ss)
	- population ftness (neg. bin; offspring by all adults averaged over populations and lines ~ ss + gen)
	- per female fitness (total number of offspring divided by two)
	- relative population fitness (per pop'n fitness per regime-x-gen divided by baseline fitness for pattern
##### Results
- in control condition, total reproductive output of small populations did not vary by ss
	- in pairs, small but sign. difference; monog > polya
	- per female, fewer offspring in small groups than in pairs for both ss treatments (could be due to group oviposition)
- extinctions more likely in monogamous treatments than in polyandrous ones (100% after 15 gens vs 40%)
	- suggests poly survived 2.65x longer than monog
- fitness: more rapid decline in mean fitness for monogamy compared to poly
	- divergence happened after generation four (in second cycle)
	- no strong differences among stressor types on fitn decline
##### Discussion
- prev. work in Triboleum has shown sexual selection -> more male reproductive competitiveness (?), more female abilty to handle sexual conflict (???), increases ability to cope with inbreeding stress, stronger competitors with conspecifics
- could the difference in performance be explained by differences in inbreeding? evidence suggests no (what the heck is going on in this paragraph?)
- support for the genic capture theory, s/t reproductive success in face of competition is achieved by individuals in better condition (what the heck does this mean? genome-wde variation across traits?)
	- sexual selection + genic capture would counter the cost of making males who don't reproduce in polyandry (because the most fit male would reproduce meaning on sum offspring is more fit?)
- para on the vortex: two predictions of the EV model (logarithmic scaling time of pop size to extincton and also increasing rate of decline closer to extinction) not supported (sigmoidal rate of extinction and flattening of fitness declines)
	- could be genetic purging

Okay... this paper is very confusing. Will need a re-read. Extinction vortex plays only a small role it seems; the sexual selection stuff seems to be more important. Were the 95 generations just to establish an evolutionary history? The part about sexual selection and implications of competition among males (this genic capture stuff?) may be related to intraspecific competition. At least there is some stuff about the extinction vortex and predictions of that model (which do not come true here). Does the experimental procedure really fit in with the idea of the vortex? Seems convoluted and confusing, esp. the way individuals are sampled in founding the next generation.

### Yamada, Y. (???) Tribolium as a biological model in quantitative genetics. (???)

###### Trib. mating and genetics, etc.
- 10 (small) chromosomes incl. XY sex chromosome
- "mating is polygamous in nature"
- Quantitative traits include: larval weight, pupal weight, adult weight, development tme, egg number
	- heritability of egg number: "low to medium", high heterosis
	- egg number +correlated with larval weight
- Cite Orozco 1969: heritability for egg laying about 20%, responds to environment

Hmm... weird, where was this even published? Consider the Orozco 1969 piece. Beyond that probably no use.

### Lavie, B, Ritte, U, and Moav, R. 1978. The genetic basis of egg lay response to conditioned medium in the flour betle, Tribolium castaneum. Theoretical and Applied Genetics, 52, 193 - 199.

- Is rate of egg laying in Tc an "independent genetic character"?
- Krause and Bell 1972 say that 75-80% of penotypic variation in fitness is from egg production (so selection should favor higher rates)
	- but, we would expect egg production to be involved in adapting to changing environments and thus it may be good to have fewer eggs in certain cases (e.g., overcrowding)

Is this going in the right direction? This seems to be about selection and treats that differently than the assertion that egg numbers is a quantitative trait... maybe first consider Orazco above. Maybe return to this.

### Fagan, W.F., and Holmes, E.E. 2006. Quantifying the extinction vortex. Ecology letters.

- Observation studies during actual extinction events are rare (but see, Belovsky et al. 1999 and Schoener et al. 2003)
- This study: compiled studies looking at vertebrate extinctions (publshed and unpublished), testing two predctions
	1. In extinct populatons, lifetimes scale logarithmically with populaton size (see Lande et al. 2003)
	2. Extinction vortex characteristics, feedbacks,
		a. proportionately larger declines in abundance as population nears extinction
		b. increase in inter-annual variability in changes in population size for smaller populations (isn't this just stochasticity?)
- Dataset: 10 annual time serieses 12-21 years long of vertebrate population declining to extinction
	- note that the short time-frame may bias extinction times to be sooner compared to pool of all extinct populations
- Metrics from each dataset:
	- smallest size from which an increase occurred
	- population size from which smallest decline occurred
	- population size one year before extinction
	- n.b. when calculating growth rates N_final = 1 (good idea!)
- EV tests:
	- (1) did year-to-year rates of decline increase as extinction neared
		- i.e. dd ln(N_t / N_{t+1}) depend on t?
		- if positive slope with t, supported
		- also does it depend on N_{t+1}
	- (2) did year-to-year rates of decline become more variable closer to extincton
		- variability is square of residuals from regression above
	- (3) did the "value" of a population ??
		- time greater than/less than 10 to see if dynamics are different closer to extinction
##### Results
- In four populatons, no increases once populations went < 50
	- in two populations, increases < 10
- Rate of decline to extinction from population peak did not vary by population size (no dens. dep.?)
- Last peak was 70 (incl. salmon)/ 49 (not incl. salmon)
	- mean last census size before extinction was 8 (2-26)
- Population size: logarithm of pop size prop. to time to extinction (i.e., time to extinction exponential with pop size)
	- initial populaton size also log. relationship w/ time to extinction
- Year to year population decline grow larger closer to extinction
- More variable population size changes closer ot extinction
- The slope in the log curve for tau < 10 steeper than tau > 10 (i.e., steeper decline closer to extinction)

Interesting. EV hypothesis supported here. For our purposes, looking at log(N_t / N_{t+1}) and regressing it against time to extincton seems like a good idea. Some difficulty in here because our model also has rates of adaptation that will change the population growth rate... Can also think about regressing change in genotypes. Extinct/surviving - would this work as a predictor?

### Rogers, A., and Prugel-Bennett, A. 2000. Evolving populations with overlapping generations. Theoretical Population Biology.

- One-locus models of adaptation with overlapping generations: Moran 1958, Ewens 1979
- Multi-locus models of overlapping generations using cumulants or momnts of population efficacy distributon: Woodcock and Higgs 1996, Prugel-Bennett 1997
	- analytical solutions for weak selection
	- effects of genetic drift are twice as large with overlapping generations compared to non-overlapping ("extra randomness with which population members die")
##### Model
- N haploid organisms, L genes with positive/negative alleles s/t efficacy is E = \sum_i^L S_i for S_i = -1 or 1
- Individuals are included in the next generation according to W_j = exp(b E_j) / Z where b is selection pressure (b = log(1+s)), Z is selection pressure-weighted sum of fitnesses
- (not that b' = b / sqrt(L) to scale by number of loci... same scaling as our model)
- Each time step, population member selected to reproduce asexually (offspring identical genotype?) and another population member randomly selected to die (uniformly?)
- Mutations occur - bdirectional at rate g
- Track dynamics using cumulants of the distribution of E
	- Prugel-Bennett and Shapiro 1997 use cumulants to model non-oerlapping generations
- Expressions given for both overlapping and non-overlapping generations, first four cumulants
	- the overlapping generation equations have an extra factor of N in denominator (due to extra reproductive event?)
- Uh... the cumulants are slightly different over time. I guess K1 reaches a smaller steady state and K2 is always slightly smaller...
- Some stuff on mutations...

Okay... the discussion here is two paragraphs. What the heck? I can't rally interpret what's going on here. Very strange. But I like the idea of studying the cumulants. Is this easy to do? Is it amenable to our model? Note the efficacy here is a haploid version of our genetic model although selection is handled slightly differently.

### Klausmeier, C.A., et al. 2020. Ecological limits to evolutionary rescue. Phil Trans R Soc B.

- Review of moving-optimum models: Kopp and Matuszewski and supplements to (this) text
	- Needed for rescue: populaton phenotype closely tracks environment
- Missing from the optimum model: density dependence, species interactions, non-quadratic relationships between phenotype and fitness, limits to the niche
- In continuous models, fitness r is defined by r = r_max - (x-E)^2/(2*s^2) (discrete models have W = e^r)
	- trait value \bar{x} evolves as product (dr/dx)V_g (for V_g genetic variance) from Lande 1976
- for constant rate of change, d, environmnt is E(t) = E0 + dt
	- for x1 = E - \bar{x}, dx1/dt = d - V_g(dr/dx)
	- equibrium lag has V_g(dr/dx) = d i.e. adaptation rate tracking the environmental change
		- can derive conditions for persistence from this
		- although note that a relaxed assumption of constant environmental variance complicates things... see Burger and Lynch 1995
##### Complication 1: non-quadratic fitness functions
- Quadratic phenotypic-fitness relationship comes from a Taylor expansion of fitness around the optimum, which isn't guaranteed to hold up when populations are far from phenotypic optimum
- Osmond and Klausmeier (2017?) show that non-quadratic fitness functions will cause deviaton from expectation (e.g., evolutionary tipping points)
	- (Fig. 2 - having trouble relating panels b, c)
	- Will want to check out Cotto et al., 2019 for some of this stuff with stage-structured populations (GODDAMNIT)
##### Complication 2: population dynamics
- Some models looking at NDD in moving optimum models: B&L 1995, Polchova et al., 2009, Osmond et al. 2017
- "Ecological quantitative genetics": N' = r(x,E,N) * N (see Klausmeier et al., 2020, Geritz t al., 1998)
- Also, trait-dependence on births and deaths (note: not included in our modl)
- oh my god I'm having a breakdown this shit is important but I just can't follow along with it... my career is over. It's over. I'm done. I'm never going anywhere. It's over.
##### Complication 3: Community context
- Prior moving optimum models have considered competition and predaton, no mutualisms yet
- No study of community adaptation yt
##### Complication 4: limits to the fundamental niche
- Central assumption of moving optimum models: the fitness landscape doesn't change shape ever
	- But there may be fundamental limits to how much/far populations can adapt [...]
	- The shape of the fitness landscape likely depends on environmental variables

Okay... I skimmed part (3) and skipped part (4). Good rundown of the moving optimum models. Go back to density dependence sections esp. for trait dependence. Look at reference for other population types. 

### Gilpin, M.E., and Soule, M.E. 1986. Minimum Viable Populations: Processes of Species Extinction. Chapter in Conservaton Biology: The Science of Scarcity and Diversity (ed: Michael Soule).

- Good background reading: Schaffer 1981
	- Differentiates between deterministic and stochastic extinction
		- Deterministic typically caused by something essential being removed/degraded or something lethal being added
		- Stochastic typically due to normal random perturbations or environmental effects (~50)
	- Also has MVP size for organisms with non-trivial life cycles
- Extinction vortices: biological + ecological processes interaction to create feedback loops
	- Can be determnistic extinctions
	- Could also leave populations smaller and more susceptible to stochastic extinctions
- *R vortex*: smaller $N$ means greater variance in r (growth rate)
	- smaller population size means demographic stochasticity
	- also changes in the population distrbution (e.g., age distribution, sex ratio)
	- can cause multi-generational effects of small $N$ (where extinction risk is larger than the sum of per-generation extinction probabilities)
- *D vortex*: increased framgentation means smaller populations $N$ in patches and smaller $N_e$ per patch
	- more of a spatial process?
- *F vortex*: smaller $N_e$ means genetic vortices
	- inbreeding and/or changes in heterozygosity/fitness
	- requires(?) manifestation of low $N_e$ on phenotype
	- e.g., loss of heterozygosity influences fitness, influencing survival or fecundity
- *A vortex*: smaller $N_e$ (and associated inbreeding, drift) means slower adaptation
	- Thus populations are increasingly maladapted to their populations
	- on this point see Frankel and Soule 1981 and Fukuyama 1983
	- has the longest timescale so may be pre-empted by other cycles
- Notably, different life histories are susceptible to different vortices
	- e.g., a quickly, mass-reproducing species may be less susceptible to reductions in $N_e$ because there are so many offspring
	- longer-lived organisms are less susceptible to vortices responsive to large environmental change

Good stuff here about deterministic versus stochastic extinction - fits in with the phased framework I had earlier. Extinction is deterministic without adaptation and/or if adaptation is too slow. But then populations are left at a small size by the time they recover and are thus susceptible to extinction.

Good thought here is that we should look at the minimum size and time of minimum size for surviving populations. It's also interesting to see how often extinctions occur after the expected trough of population size.

These various types of cycles are not mutually exclusive, also are confusingly defined. Need to think about how they manifest themselves here. Off the top of my head, we have potential for *R* vortices based on populations being small/declining and NDD may make that worse. Look at var(r) perhaps? We also have *F* as inbreeding, etc. influence fitness. I love the inclusion of *A* - the rate of adaptation *should* be slowed down by drift. Think about what is needed to quantify each of these.

### Kulbaba, M. W., et al. 2019. Additive genetic variance for lifetime fitness and the capacity for adaptation in an annual plant. Evolution, 73-9: 1746 - 1758. 

- Fisher's FTNS says rate of adaptation (in current environment) is V_A(W) / \bar{W}
	- studies quantifying V_A(W) are rare, although we know rapid adaptation exists in nature
	- rapid adaptation and evolutionary rescue require sufficent V_A(W)
	- what amount of V_A(W) is sufficient for rapid adaptation?
- Why would/would not there be high natural variation in additive variance?
	- Theoretical arguments about "incessant" selection and fixation
	- but, fixation likely takes a very long time in natural populations
	- Estimates of heritability often suggest that heritability is low
	- but, low heritability may be a result of high environmntal variation
	- Difficulty: fully pedigreed individuals are needed to avoid confounding environmental and genetic effects
- Here, measuring V_A(W) in *Chamaecrista fascculata* (annual) growing in native location
- 2-3 fruits from 200 maternal plants (each 10m apart) collected from three populations (two in MN, one in IA)
	- get a pedgree, planting seeds from respective populations into common garden in MN (no mowing, grazing, etc.)
	- five seeds per family, 15-24 replications of half-sib fam per block, [...]
- Fitness is number of seeds collected at end of yar
	- also recording: emergence, flowering, fruit (incl. elongated but immature), fruits on ground in vicinity, pedicels from dehisced fruits on plant
- Aster model: emergence, flowering, fruit production (pois), seed production (pois)
	- going out to seeds but only measuring from certain seeds... scaling up, not quite sure what to make of this, se Stanton-Geddes et al. 2012
	- somehow asters are also estimating V_A(W)...? not quite sure how this works either! is it just variance in the fitness response?
- Significant V_A(W) (as in significantly non-zero?) in all three years of planting, in all three populations
	- significant variation among years suggests parental genetic effects varying among years
- increasing in \bar{W} ranging from 0.83 to 6.12
	- two cases of rescue where \bar{W} went from < 1 to > 1
	- what are they adapting to??? or is this just fluctuation?
- Finding more V_A(W) here than elsewhere - maybe because the varous components of fitness are compounding

[need to finish]

### Hufbauer, R.A. et al. 2015. Three types of rescue can avert extinction in a changing environment. 

- Brown and Kodric-Brown's rescue effect: immigration (esp. human-facilitated) averting extinction in small threatened populations (esp. if inbred or in degraded hab)
- **Demographic rescue** more individuals in population (more resistant to stochastic fluctuations, fewer Allee effects)
	- fitness effects: may allow for more time for adaptation, or may introduce maladaptation
- **Genetic rescue** increased population fitness du to contributions from immigrants (less inbreeding, more genetic variation)
- **Evolutionary rescue** is the result of a population adapting on its own from standing variation
	- however immigrants can enhance the speed of adaptation
- What are relative strengths of these? Demographic vs./in addition to genetic processes
- This study: experimental partition of demographic and genetic effects in a challenging environment
	- Size treatment: 50 vs. 150
	- Treatments:
		- Evolutionary (no migrants)
		- Demographic (population size increased - same source pop)
		- Genetic (genetic variation added - demographic substitution)
##### Results
- Only large populations to go extinct were evolutionary ones (i.e. no immigration)
- More extinctions in small populations (as expected) including in all treatments
- In early generations, slight population growth rate advantage to demographic-treatment populations (avoiding Allee effects)
	-  in subsequent gnerations, demographic advantage fell (compared to evolutionary rescue treatment)
	- in subsequent generations, genetic effects led to substantial increases in growth rate esp. for larger populations
	- notably genetic populations had same growth rates as demographic + genetic effects
##### Discussion
- So both of these contribute to lessened extinction risk (compared to evolutonary rescue alone)
- Demographic treatment did enhance persistence with the most substantial affects immediately after the migration
	- long-term effects on population size for large (but not small) pops
- Genetic effects: lower extinction risk and also higher intrinsic fitness (pop growth) in final generation (less inbreeding depression enhanced adaptive capacity)
- D + G reduced extinction risk rel. to D alone or G alone but did not increase final size or fitness
	- demographic rescue could swamp effects of genetic rescue (weakly supported here - observed in one generation)
	- no effects of size here - similar responses of small and large pops (interpretation is here but I don't buy it)

Here, uses are (1) some empirical demonstration of rescue esp. differences between large and small popultions and (2) showing additional strategies/effects mainly from interactions with other populations.

### Stewart, G.S., et al. 2017. The power of volutionary rescue is constrained by genetic load. Evolutionary applications.

- Here: evolutionary rescue is from standing variation and no migration
- Genetic load: higher in haploids than diploids b/c deleterious alleles are subject to selection (note: not applicable in our model b/c of dominance - deleterious alleles are acted on by selection)
- Genetic load: drift load + inbreeding load
	- inbreeding is estimated through comparing fitness with random mating in population at large
	- drift load is estimated through comparing fitness with outcrossing with (other?) populations
	- inbreeding turns to drift as heterozygosity is lost
##### Methods
- Trib studies: favorable and unfavorabl environments
- Two founding sizes: N_0 = 10, N_0 = 40
	- sex ratio not controlled
- Non-evolving and evolving treatments
- Measuring load: H is that by gen 8 inbreeding and genetic load combat effects of evolution
	- test this by mixing individuals among popultions for mating and comparing with unmixed evolving populations
	- (high difference here is heterosis, evidence of drift load)
	- didn't measure inbreeding load (can't know sibs)
	- able to do this for evolution treatment but not by size
- Extinction analysis: binomial GLM
	- data augmentation for favorable environment (add one extinct trial) because three of four treatment groups did not go extinct
	- analysis: popsize x evolving x environment + fixed block
- Of surviving populations, evaluated population size (repeated measures w/ previous size as predictor) in all time steps?
	- also looked at log growth rates using ...
##### Results
- >50% of poor environment populations went extinct
	- no interactions (power?) but env., fsize, evolution influence extinction risk
	- small populations have 5x greater chanc of extinction compared to 40
	- evolving populations have ~1/5 extinction risk of smaller ones
- Poor environments: large populations quickly returned to smaller size, evolving populations larger than non-evolving ones in most time steps
	- in favorable environmnts, evolving populations were only larger for first few generations before shrinking to approx. same size as large pops
- Admixture: environment mixed populatons performed better than evolving and control populations in both treatments (save for unfavorable evolving population for which there were insufficient pops)
	- suggests heterosis present in all environment types
- What's happening here: in the evolving populations, populations evolve at first but disappear ~7 generations fter founding du to non-adaptive genetic processes e.g. inbreding and drift (esp. due to small size)
	- supported by heterosis study, where the isolated populations likely had different alleles reach fixation
	- similar results by Agashe et al. 2011 (Evolution)
	- may also be due to epigenetics or grand-maternal environmental effects

Okay notably not rescue (and the evidence of adaptation sems kinda weak) but good empirical demonstration of accumulation of load during the adaptation process, limiting adaptive capacity.

### Yun, L., and Agrawal, A.F. 2014. Variation in the strength of inbreeding depression across environments: Effects of stress and density dependence. Evolution.

- Strength of inbreeding depression is variable... why?
	- this can include multiple measurements of the same population producing different values
	- in this case, if degree of inbreeding F and polymorphism degree p are constant, h or s would cause differences in strength of inbreeding (why is h not constant?)
- Is inbreeding selection stronger in harsher environments?
	- Agrawal and Whitlock 2010 suggest that density dependence should influence fitness
	- density dependence and intraspecific compettion has effects that disproportionately affect the weak
		- e.g. with food scarcity stronger genotypes may be better able to find food, weaker genotypes search longer, more predation risk
	- **note: this is density-dependent fitness, no?**
- Here: measure Drosophila melanogaster in lab experiments with varying stress and density dependence
##### Methods 
- 90 lines of *D. melanogaster* with 22 environment types (thermal regime, chemical compositon, food abundance, water availability, substrate texture)
- Four assays:
	- 130 outbred + 130 competitors (DsRed) (normal density + outbred - FO)
	- 130 inbred + 130 competitors (normal density + inbred - FI)
	- 100 outbred + 100 competitors (low density + outbred - LO)
	- 160 outbred + 160 competitors (high densty + outbred)
	- Inbred flies have F = 0.4
- Assays x environment allow quantification of stress, density dependence, inbreeding depression
	- Stress: absolute reduction in absolute fitness below maximum (calculated here for focal value)
	- Density dependence: reduction in fitness with increasing density (calculated here with regression)
	- in all cases fitness is survival
- Differences evaluated by standard parametric testing and also likelihood analyses
##### Results
- NDD and stress non-significantly correlated (p = 0.47)
- Correlation between inbreeding and stress was non-significant (p = 0.18)
- Inbreeding depression was significantly correlated with NDD
- In a joint model, stress x density dependence interaction 
##### Discussion
- Inbreeding depression is stronger where NDD is stronger
	- this may be by increasing the *mean* strength of selection
	- could also be by reducing dominance (?) or covariances between genotypes
- Could be due to NDD caused by increasing competition, where 
genotypic differences are accentuated
	- Weiner 1985 says at high density, competition leaves losers with relatively fewer resources
	- Phenotypic variability (Waller et al. 2008) - inbreding depression greatest where there is more phenotypic variability (??)
	- Long et al. 2013: variation within groups predicts variation between groups (outbred vs. outbred compared to inbred vs. outbred)
- Fox and Reed 2010 found a significant relationship between stress and inbreeding depression, but this study did not. Why?
	- could just be because they were not accounting for NDD
- Stress vs. novelty
	- past selection reducng polymorphism at sites strongly selected for in normal environment
	- past selection will not reduce polymorphism at sites neutral/weakly selected n normal environment means there will be more varaton here
	- this means higher inbreding depression in novel environments? Presumably because bad stuff in the novel environment has not been purged yet and it can be revealed through inbreeding
- But not all stressful environments are novel
	- maybe there is complexity in envronment types requiring either few genes (simple stress) versus more of hte genome (more complex)

Interesting paper. Argues that NDD (via stress increasing competition) increases the severity of inbreeding depression. This is because the good genotypes are good at excluding poorer genotypes (which inbred individuals are more likely to be) under intense competition. Note that inbreeding depression applies to individuals and not necessarily at the population scale.

I don't think this applies to our model because competition affects all individuals equally. I'm not immediately what would be the case if our model did have this. Would inbred individuals just be purged and the population would be fine so long as there were fit indivduals? Ugh - think about this and consult with group.

In general this looks like it's all about density-dependent fitness which is not quite how our model operates. Worth having a paragraph on this?

### Frankham, R., et al. 1999. Do population size bottlenecks reduce evolutionary potential? Animal Conservation. 

- Pop gen theory: a single generation bottleneck reduces heterozygosity by 1/(2Ne) for Ne effective population size
	- single pair bottleneck reduces heterozygosty by 1/4 (also applies to quantitative traits)
	- heterozygosity as measure of genetic variation for evolutionary potential
- Can there be cases of higher genetc variation in bottlenecked vs. non-bottlenecked populations?
	- non-additive genetic variation
	- although, bottlenecks could expose rare deleterious alleles that are partially recessive so this variation exposed may be non-adaptive
- Do population bottlenecks influence evolutionary potential (in a changing and challenging environment?) and if so how?
##### Methods
- *D. melanogaster*: in single pair bottlenecks for one (1B treatment) and three generations (3B)
	- 1B: 15 reps, increased from single pair to size 500 for 18 gens
	- 3B: 5 reps increased to 500 for 18 gens
	- Inbred: 5 reps with 35 generations (!) of full-sib mating and inbreeding coefficient of 0.999
- Evolvability: 250 pairs of flies (w/in a line) in a standard medium cage for two weeks then several weeks of successivel saltier media
	- continually changing environment
	- extinction: two successive weeks with no larvae in plots
	- salt concentration in first of two weeks w/ no larvae is concentration of extinction
	- terminal test to evaluate if adaptation occurred by taking some of the experimental stocks (?) in medium with no stocks and compared to an unevolved stock (counting egg laying)
- Statistical comparisons
	- outbred base pops versus single pair bottlenecks
	- outbred base pops against eachother
	- single pair bottleneck versus three-generation bottleneck
	- between bottlnecked populations after size restriction versus 18 gens after
##### Results
- 1B had more extinctions than outbred populations
- No differences among wild populations
- No differences between 1B3 and 1B18 (after/late after size reduction)
- Relationship between extinction concentrations and proportion of genetic variation remainng (how was gnetic variation determined??): increasing genetic variation means increased concentration of extinction
- Longer bottlneck/more inbreeding means quicker extinction (and lower extinction concentration)
	- i.e. no bottleneck means can withstand lower levels of change
- Unevolved populations unable to produce viable offspring at 6% salt while the evolved populations could (at 10% of the salt-free concentrations)
##### Discussion
- The 1B treatment did have the 25% expected reducton in heterozygosity (again - how measured?) but the 3B treatmnt only had 29% reduction - far less than the expected 50%
	- could just be that bottlenecking is noisy and there were few reps
	- selection opposing fixation (interesting... how? at what stage?)
	- also happened in another study cited here
- "Population bottlenecks should be avoided whenever possble in managing populations"
	- "populations shown to have suffred bottlenecks, such as island populations, should be recognized as having a reduced ability to response to novel environmental threats"

Bottlenecking hastens extinction and lowers the amount of tolerable environmental change. Note the environmental change regime here is similar to that of Burger and Lynch, not G&H.

Uses to us: helpful in thinking about the extinction vortex. Populations shrinking too much are likely to enter the vortex, which has inbreeding/bottlenecks. This reduces their adaptive capacity.

Worth also thinking about how the vortex may affect initially small populations and/or low diversity more. Do they face additonal additive effects of bottlenecking, compounded effects, or have they already suffered the worst genetic effects simply by virtue of their size/diversity? I think this would be a good angle to use in Discussion that would link together ideas of the vortex and our treatments. 

### Hansson, B., Morales, H., and van Oosterhout, C. Comment on "Individual heterozygosity predicts..." Science.

- Comment on Scott et al. 2020 about translocation of desert tortoises over 18 years 
	- 10K tortoises released, of survivors there is much higher genome-wide heterozygosity
	- Scott et al. suggest that we should prioritize heterozygosity in managemnt/cons
- Hansson et al. are arguing here:
	- This is oversimplifying
	- Advice may not be best way to establish viable pops
	- Translocation should also include assessment of genetic load
- Positive correlation btwn heterozygosty and fitness are found in a number of taxa
	- Pop viability is influenced by load but also by expression of recessive deleterious alleles in homozygosity ("exposure of genetic loads")
	- study cited looking at genetic load vs. demographic history and extinction risk in IUCN
- Metapopulation simulations under the advice to maimize genomic variability and heterozygosity	
	- translocating individuals with high diversity does increase fitness
	- but, benefits are transient, increasing load
	- translocating those with the smallest load produces both initial and sustained high fitness because deleterious mutations have been purged
	- can strike a balance by combining individuals with lowest load and highest "genetic uniqueness"
- Also the heterozygosity advice ignores cological heterogeneity and local adaptation
	- here rather than prioritzing heterozygosity it would be better to consider closely-located donors

Sims... not really sure how to interpret although admittedly I didn't look at very hard. Seems like a lot of this is due to good performance of the "random" selection of individuals for translocation. Load + drift seems optimal. Two-step strategy.

Translocation is something akin to rescue (colonization of a novel habitat)? Suppose how this would happen in our model is heterozygosity/diversity increases speed, but reducing load reduces distance so speed is less important. The heterozygosity without reduced drift increases the odds of maladaptive alleles reaching fixation - during rescue this is especially important due to the chance of populations becoming small.

Fascinating stuff!

### Pauls, S.U., et al. 2013. The impact of global climate change on genetic diversity within populations of species. Molecular Ecology.

- How climate change will/can influence genetic diversity:
	- change the distributions of genetic variants (due to species range changes)
	- change levels of phenotypic plasticity
	- influence adaptation
	- (in most cases we expect reduced variation)
- Species avert climate change impacts through plasticity or adaptation
- In context of range shifts
	- Leading edge tends to have low neutral diversity (recurring founder effects)
	- Leading edge also though has most persistent alleles (trailing edge has alleles more likely to go extinct)
	- Leading edge: can have selective sweeps for newly advantageous alleles, including sweeps on diversity
	- Gene surfing can happen: drift means that certain alleles are over-represented
	- Range contractions should decrease effective population size and introduce bottlenecks
		- speed of contractions: Arenas et al. 2012 shows slower range contractions may have lower diversity (not sure why)
	- Climate change could cause hybridization events (possibly through moving ranges)
- Micro-evolutionary responses in stationary populations (e.g., local adaptation)
	- see Riddle et al. 2008, Hoffman and Sgro 2011, Peischel and Kirkpatrick 2012
	- not only could reduced genetic variation at certain loci occur, but reduced variation at other hitchhiking parts of the genome (see Charlesworth et al. 1997)
	- see Vogt et al. 2007 showing exposure to stress can reduce neutral diversity through adaptation and drift (see also Nowak et al. 2009, Vogt et al. 2010)
- Speed of adaptation and climate change
	- Global climate change is rapid and requires strong selection pressures
	- Not only mean trends in climate but changes in variance will have influences on variation
	- See Drake 2005 for something about climate change and population growth rates
	- Avolio et al. 2012: genotypic diversity in a prairie grass decreased appreciably with a decade of increased variation in precipitation
		- but also possibly an increase in plot-level genotypic diversity? possibly due to increased spatial and temporal variation in moisture condition in plots
		- one conclusion is that with sufficiently high initial variation, more climate variation can cause genetic differentiation through spatial and temporal variation
	- Another possibility is that climate change selects for dispersal which ultimately homogenizes populations and reduces drift, increasing genetic diversity
- Feedbacks and multiple stressors
	- Reductions in effective population size through diminished census population size, changing sex ratio, reduced fecundity, increased mortality (Lande 1988, Reed 2005)
		- This will further reduce diversity through drift, inbreeding (Nowak et al. 2009)
	- Yes - reference to the extinction vortex here
[Some applied stuff I skimmed]
- Neutral and non-neutral diversity and adaptive capacity
	- "There is an inherent link between the standing genetic variation in a population and their potential for adaptation, but this link is not well documented in the [climate change] literature"
	- Jensen et al. 2008 (trout) and Fischer and Karl 2010 (butterflies) showing life history trait responses on an environmental gradient...
- Phenotypic plasticity as a pathway for "adapting" to changing conditions
	- Can be done by maximizing an individual's ability to respond to changing environments, even if coming at a cost to some other function such as reproduction
	- Phenotypic plasticity should be highest in longest-lived organisms with long generations (Canale and Henry 2010) and low annual reproduction
		- (but see trees as a counter example? Petit and Hampe 2006; argument is I guess that population sizes are typically very large)
	-  Canale and Henry find that a lot of adaptive capacity is due to plasticity rather than genomic adaptation??
- Gene expression: some empirical examples (re-check?)
- Next steps/improvements
	- Neutral variation is good, but it would be better to get a more functional genomic understanding, get candidate loci and assess their genetic diversity
	- Space-for-time substitutions are currently widely used, but maybe we can use natural archives w/ DNA, museum collections, cryobanking etc. for climatic variation

Lots of good stuff in here.

Most important is some empirical refs looking at loss of genetic variation in certain contexts. Nowak et al. 2007, and the Vogt. et al. papers look useful here. Drake 2005, Riddle et al. 2008, Hoffman and Sgro 2011, Peischel and Kirkpatrick 2012 also are useful here.

The range shift and life history pieces were interesting for other projects. General rule seems to be that diversity will shift throughout a range when there are range shifts. Also the idea that slower-growth, longer-lived species are more likely to rely on plasticity and have less genomic capacity for adaptation to global change. 

Very intrigued by this idea of selection for more dispersal (to escape bad habitats) and ultimate homogenization. Drawing this out slightly, if there's selection for dispersal, then presumably that will happen in sites that are low-quality? Would one expect to see selection for dispersal in high-quality sites? Dynamics here? Interesting to think about in context of the assisted migratino study.

Would also be good to cite in introduction, making case for why changing genetic diversity is important to incorporate into models.


### Nowak, C., et al. 2009. Rapid genetic erosion in pollutant-exposed experimental chironomid populations. Environmental Pollution.

- Predictions that chemical exposure would reduce genetic variation - is it so?
	- Some refs of "reported" genetic erosion, but these don't show strong evidence of it as a result of stress
	- Genetic diversity is influenced not only by past and present selection but also by migration and chance and field studies do not necessarily control for these
- Here: does exposing midges (*Chironomus riparius*) to a biocide (TBT = tributylin) at an environmentally relevant concentration cause reduced genetic variation at neutral loci in a short timespan? Are there changes to life history?
##### Stats
- Larval mortality ~ TBT and Mean emergence time ~ TBT tested with nested two-way ANOVA (repeated measures)
- Lienar trends in life history (mortality, emergence time, fertile egg ropes (?) per female, population growth rate) and pop gen parameters assessed with a Mann-Kendall test
- Heterozygosity, alleles per locus, polymorphic loci, and deviations from Hardy-Weinberg tested (microsatellite markres involved here)
##### Results
- Overall mortality w/ TBT: 44%, w/o TBT: 28%
	- mortality was only higher in TBT groups in 4/10 generations
- Higher emergence in TBT groups in all but two generations
- Fewer clutches per female in TBT groups
- Allelic richness did decline in control populations over time
- No temporal trends in any of the life history groups (i.e. no selection)
- Number of alleles decreased from 3 to 2 in both stressed treatments (trend looks similar in controls though?)
	- monomorphy (fixation?) at one locus in one control group, while TBT groups had monomorphy at 3/5
- Observed heterozygosity decreased in both stressed groups, no corresponding trend in controls
	- loss of (expected) heterozygosity was slower in controls, ending with higher degree of heterozygosity
- Because there was no evidence of a trend in any LH trait and therefore no selection, most likely culprit for loss of diversity is reduced effective population size (because fewer offspring in each generation)
	- Gene flow would have increased N_e

Use for this paper: good controlled evidence that loss of diversity happens in populations exposed to a stressor in a setting where all else (incl. especially gene flow) is controlled. This appears to be due to the lowered effective population size associated with the reductions in fecundity (it looks like actual census population size was held constant over time though). It's possible that selection is occurring but there isn't enough statistical power to detect it. Good empirical demonstration of link between stress/conditions and genetic variation!

### Rego, A., Messina, F.J., and Gompert, Z. 2019. Dynamics of genomic change during evolutionary rescue in the seed beetle *Callosobruchus maculatus*. Molecular Ecology.

- Probability of evolutionary rescue under varying conditions has been well-established, but evolutionary dynamics and genomic consequences are less well-known
	- Some previous work suggests that mutations of large effect (could) have a large influence on rescue during early stages (b/c of distance to phenotypic optimum)
	- Some work showing that standing variation is helpful but also soft selective sweeps could be occurring during rescue
- Some outstanding questions:
	- Does rescue depend on a few, or many loci?
	- Is selection on loci consistent throughout the process? (not quite sure what this means)
	- Does selection cause alleles to fix? Or are the changes more subtle?
	- Answer these questions with a legume-eating beetle moved to a marginal host
##### Study
- Focus on one line; other attempts to establish on lentil hosts failed
	- but three independently derived long-established lines (L1-L3) and another failed line (L11) can be used
	- L14 is a new line (?)
	- All lines established from the same ancestral population, reared at first on mung bean for >300 generations
	- L1-3 had rapidly increasing survival over a few generations along with other life history changes
	- Genomic analysis of L1-3 in Gompert and Messina 2016) but only after ~80 generations on new host
	- L14 is one extinct population on lentils that did experience a transient rapid increase in population size
	- L15 placed in each of 40 petri dishes with 100 lentil seeds, then seeds with a hatched egg were isolated for emergence over F5, F10, F20 generations
		- Survival analysis?
- Genetic data from above: 48 beetles per generation in L14 founders and F1-4
	- then L14 line split in to sublines (A and B) and generations subsampled from each...
	- Look at >21K SNPs - not necessarily finding loci responsible for adaptation, but yes finding SNPs affected by selection through linkage disequilibrium in causal variants (how??)
	- SNPs in L1-3 described in Gompert and Messina 2016
- Estimate N_e for L14, use that to parameterize and test a null model of evolution solely by drift
	- here: identifies SNPs most likely to be affected even if indirectly through LD
	- This yielded 198 focal SNPs (of which 188 were available for L1-3)
- Selection strength estimated for a few generations for L14 (didn't read carefully)
- Also estimating genomic change on L1-3 and L11 to see if these patterns of initial change (L11) and final outcome (L1-3) look similar to those of L14's rescue process
	- if yes: repeatability of dynamics in genome
##### Results
- (Skipped for time being)
##### Discussion
- In 198 focal loci, strong selection occurring
	- different loci associated with early versus later stages of adaptation (i.e. selection acts differently during different periods of rescue process)
- Strong selection in 10 generations causing huge increase in survival (1% to 90%)
	- during this time population went through bottleneck where N_e < 10
	- strong selection in >100 SNP markers in 6-16 high LD clusters
	- results are consistent with some previous results (Orr 2005, Bell and Gonzalez 2009) as well as theory suggesting mutations of large effect have disproportionate impact
- For >100 SNPs, a minor allele reached frequency > 90% within 16 generations
	- results here differ from other eukaryotic studies, mostly on *Drosophila*, showing subtle shifts in allele frequencies and not selective sweeps
	- may be result of the very strong selection pressure and the fact that many of the other studies did not feature strong enough selection for there to be decline
- Varying selection strength on loci... could be due to
	- rapid fixation occurring, after which selection on these loci may have ceased
	- epistasis
	- direct selection on causal variants could still be constant but this study wasn't looking at causal variants
- The patterns in allele frequency and selection in the 198 focal SNPs in L14 looked different than those of L1-3, L11
	- suggests that evolution is not "repeated" (same/similar path) during adaptation to new plant
	- high sensitivity to initial conditions? chaotic dynamics...

Very cool. Rapid adaptation and fixation occurring on some loci, esp. during bottleneck. Supports some of the genomic results of our model. 

### Bijlsma, R., and Loeschcke, V. 2011. Genetic erosion empedes adaptive responses to stressful environments. Evolutionary Applications.

- Ability to cope with changing and stressful environments relies on individual phenotypic change (e.g. plasticity) and also genetic diversity available for adaptation
	- Possible empirical examples of rescue? See Roberts and Andre 1994 (DDT) and Dahlgaard et al. 2001 (temperature stress)
	- migration of course also a response to environmental change
	- plasticity also a response to rescue (although not always adaptive?)
	- of course this response relies on having sufficient necessary genetic variation in a population - often not present
		- Orr and Unckless 2008 and Teotonio et al. 2009 - adaptation relying on standing beneficial variants and not necessarily production of new variants?
##### Habitat fragmentation and reduced population size
- Population losses (and isolation) increase drift and inbreeding
	- "because the force of genetic drift increases with decreasing population size, the potential to respond to natural selection will, in general, decrease with decreasing population size"
	- Selection coefficients below (1/2N_e) drift is stronger than selection and genetic dynamics determined by genetic variation rather than the alleles themselves (Kimura...)
	- but also purging will happen for deleterious alleles with large effects on fitness
	- because alleles becoming fixed through drift is equal to its initial frequency, stress-avoidance/coping alleles that are rare under beneficial conditions are unlikely to reach fixation
##### Inbreeding, plasticity, and stress tolerance
- "The magnitude of inbreeding depression generally increases under adverse environmental conditions"
	- increasing inbreeding depression in Drosophila under different types of stress (Kristensen et al. 2008, Joubert and Bijlsma 2010)
	- Fox and Reed 2011 has a metaanalysis on this
	- Genotype x environment interactions proposed: could be due to fixation of maladaptive alleles when the environment is favorable that gets worse when the environment is unfavorable (would have likely been purged otherwise)
	- Also, potentially more variation in fitness among inbred individuals, meaning a greater proportion of them would not be viable under hard selection where minimum viability is needed to survive
	- Biljsma et al. 2000 showed experimentally that extinction was more likely in small populations under stress conditions
- Plasticity has a genetic basis and there is genetic variation for plasticity (i.e. plasticity as a trait)
	- Individuals in small populations may have less plastic response
	- Mixed results of inbreeding on plasticity (negative effects but also little effects)
	- Genetic erosion may hamper plasticity, then
##### Genetic effects
- Some empirical stuff about loss of neutral diversity in experimental and natural populations: Frankham et al. 2002, Johnson et al. 204, DiBattista 2008, Hoeck et al. 2010
- Genetic erosion can be countered by increasing gene flow (esp. in fragmented populations) and mutations in local population
	- new mutations require large population sizes (remember N_e is smaller than N typically)
- Adaptive response in smaller populations is slower due to drift and inbreeding
- Rescue process means smaller populations which contributes to genetic erosion... (our argument!)

Good review that should absolutely be used to highlight the reduction in genetic diversity that occurs in (1) stressful environments and (2) small populations. Would go well even in introduction.

The stress-inbreeding argument is interesting and new (GxE interactions). Different way of thinking about it than what I had been thinking previously. Here and throughout hte paper the inbreeding here is more mechanistic than our model. 

### Hao, Y., Brockhurst, M.A., Petchey, O.L., and Zhang, Q. 2015. Evolutionary rescue can be impeded by temporary environmental amelioration. Ecology Letters.

- Most rescue studies look at monotonic environmental deterioration, but fluctuations are important too
	- Could act as "temporal refuge" slowing down environmental decline
	- Or, could relax selection for stress tolerance, slowing the fixation of beneficial mutations
	- Note that "later steps of adaptation are contingent on the success of earlier ones" so genetic effects would be compounding (both would be)
- Here: looking at rescue in a lytic bacteriophage previously demonstrated to have adaptive capaity in face of high temperatures
	- here: experimental evolution in gradually warming environment with/without "amelioration" via lowered temperatures
	- also looking for size interaction here
		- hypothesis: positive effects of amelioration would be better for small populations avoiding demographic stochasticity
##### Methods
- Phage, typically failes to reproduce above 30C but can adapt, grown in microcosms that form plaques w/ 10^7 bacterial cells (in an ancestral state)
- Study: phage populations either have 
	- no amelioration: temperature increasing in stepwise fashion from 29 to 31.4C in six steps
		- each temperature has three transfers/passages
	- periodic amelioration: second transfer in each temperature block is 0.4C lower
	- Four incubators for each regime, 12 lines in each incubator
		- six have 1% of phage cultures transferred in each step, six with 0.1% transferred in each step (large vs. small bottlenecks)
- Analysis of time to extinction analyzed with a survival regression model
	- demographic recovery: look at change in population size for each chapter
##### Results
- More extinctions, and quicker extinctions, in the amelioration treatment
	- In early stages of experiment ameliorations led to increased population sizes, but had worse decline in later stages
- Progressive/multi-step adaptation study... not sure I 100% follow but basically populations tolerant to a mildly stressful temperature to a more stressful environment didn't grow
	- so... multiple beneficial mutations required in this case?
	- conclusion here is that the amelioration means less opportunity for some of the mutations to reach fixation, slowing adaptation
- No evidence that the beneficial mutations to stress incur a cost in the benign environment
	- so amelioration relaxes but does not reverse adaptation
##### Discussion
- In this study, early demographic recovery due to amelioration is outweighed by later evolutionary shortcomings
	- (wait... demographic benefits would transfer over during the transfers of 0.1-1%, right?)
- Amelioration could be an evolutionary trap akin to Ferriere and Legendre 2013
- Temporal genetic rescue... large population size in ameliorated periods provide a demographic source (relatively at least) but also impede fixation through maladaptation

Really like the framing here - great interesting contrast between demography and evolution. Also thinking about the early vs. later stages in adaptation... does it seem similar to other patterns that early stages are disproportionately important (does it? geometric/multiplicative process?).

The successive evolution is interesting. Does it matter in our model? I guess not because loci are iid and not sequential? Does this sytem have sequential evolution...? Hard to think about in terms of our model because we have a single shift - here the effect of the amelioration would depend on its direction/magnitude in relation to the current state of the population.

I wonder if here the fact that the amelioration led to two successive increases rather than one had any impact here? 

The demographic recovery - would provide more benefit through increased mutations, but also, does it depend on the population size? The recovery for populations that are currently very small would be very helpful, one would thing. Come to mention it - I'm not really sure what the result involving size is/was here? I guess the difference between proportions of ameliorated/static populations is smaller for small populations than large ones... but the figure shows no survival past generation 15 for small ameliorated populations in a 21 (?) generation experiment?

NDD here - not sure what effect it has here. Maybe makes populations smaller such that the effects of amelioration would be (better? worse?) than DI populations... populations may be smaller during the ameliorated time periods. Hmm...

The Willi et al. 2006 paper looks useful. What's up with the Kirkpatrick and Peischl paper? Mutations, space and time, interesting but perhaps less relevant than some?

### Schluter, D., Price, T.D., and Grant, P.R. 1985. Ecological character displacement in Darwin's finches. Science.

- "Ecological character displacyement occurs when morphological differences between coexisting species are enhanced as a result of competitive interactions between them"
	- presumably phenotypic differences within a population as well?
- Darwin's Galapagos finches (Geospiza fortis, G. fuliginosa) are sympatric
	- varying beak size on shared islands, beak sizes intermediate where alone
- (1) Assess character displacement by comparing morphology in sympatry and allopatry, while controling for differences in food supply)
	- quantitative description of food characteristics to create an adaptive landscape (for beak morphology phenotypes)
	- using this landscape one can predict what food supply would be on various islands based on food supply alone; deviation from this (greater morphological differences in sympatry) is evidence as character displacement
	- adaptive landscape: estimate expected density of a population as a function of beak depth
		- identify seeds in a sp's diet given beak depth, get density of seeds within these limits on island, estimate finches that can be supported by that density of seeds (how is explained)
- (2) Assess character displacement by observing selection pressures on species in allopatry
	- differential survival of phenotypes (w/in sp?) affected by use of foods by other species... (?)
- Look at eight sympatric islands and a singleton island for each species
	- population density curves show distinct peak on several islands, but often distributions are polymodal
	- (wait... why is the expected density in Fig. 1 not necessarily the expected one for G. fuliginosa... last para on first page probably explains)
	- Fig. 1 - on sympatric island there is a greater difference between two species beak depths with at least G. fortis further off its expected peak... not as convinced for G. fuliginosa
	- Fig. 2 - single (no competition) G. fortis population is closest to its expected peak but for sympatric populations are displaced by more
	- Fig. 2 - G. fuliginosa has much larger beak than expected in single population and other populations are closer to expectation? but somehow this is also evidence of displacement...?
	- I'm a little confused by the directionality here... or at least, how is the food supply and predicted density factoring in here? This only makes sense to me by thinking about the sympatric populations in relationship to the single allopatric ones, but one sympatric is basically where it is expected based on the food supply (and the other populations aren't!) while the other species's sympatric is much larger than expected
		- oh - G. fuliginosa actually is observed at the larger peak when alone but I guess when paired with G. fortis it is shunted in to the smaller part... why is the expected allopatric density at the lower peak though??
- Analysis under drought: severe drought in 1977 with mass mortality; regression shows that the survivors were able to eat large seeds (larger beaks... second density peak?)
	- but a small number with small beaks did survive by exploiting what few small seeds were present
- G. fortis size arises from trade-off on seed size; when G. fuliginosa is present, then small seeds are not available in the dry season (due to competition)
	- so, G. fortis individuals experience selection for larger beaks (larger seeds) where the competitor is present
	- likewise G. fortis on Daphne (allopatric island) are more likely to eat the smaller seeds in the absence of competition

(read last two paras)
also think about this in the context of our model. The fitness landscape is monotonic... but NDD changes it I guess? There isn't a case where the less fit individuals are favored - it's just a matter of them not being disfavored I suppose?

### Willi, Y., Van Buskirk, J., and Hoffman, A.A. 2006. Limits to the adaptive potential of small populations. ARES.

- "How good is the evidence that small populations have reduced evolvability?"
	- (Houle 1992 is good source for evolvability, it seems)
- Frankham 1996: small populations have low neutral genetic variation and low heterozygosity
	- assumes that heterozygosity at markers are proportional to heterzygosity underlying phenotypically expressed traits (and therefore traits experiencing selection)
	- but, genetic variation depends on more than population size (e.g., number of loci, dominance/epistasis, mutation rates, selection intensity, mode of selection)
	- other factors influence adaptation: enviromental stress, genetic stress

##### Theory of genetic variation and population size
###### Single locus
- allelic diversity may be a better measure of genetic variation than henterozygosity
- with only drift and mutation, heterozygosity increaes monotonically with population size
	- both because of the relationship to effective population size and the fact that more individuals means more mutations
- single-locus models with selection and additive gene effects (single-locus but additive effects?): heterozygosity declines because the favored alleles will reach fixation
	- selection is more efficient in large populations than small ones (because of the lack of drift?)
	- so, elimination of deleterious mutations increases with N
	- "A locus under selection begins to behave like a neutral locus when s <= 1/(2N_e) (Wright 1931)"
		- i.e. with weak selection strength (or small population size) loci are more subject to drift?
		- "Therefore as population size declines, stronger selection is required to prevent loss of rare beneficial mutations due to drift, and to ensure elimination of deleterious mutations. Consequently, genetic load becomes a problem"
- with enough dispersal, isolated populations behave as if panmictic
	- for isolated populations with low dispersal, loss of heterzygosity depends on number of demes/populations (probability of fixation happening at all demes?)
##### Polygenic models
- Assuming neutrality and additive effects, equilibrium V_A of a quantitative trait increases linearly with N_e (see Chakraborty and Nei 1982, Lynch and Hill 1986)
	- because variation is eroded in small populations
	- but not all small populations have low variation; drift at low N_e can produce a lot of different outcomes
- Under selection, equilibrium V_A depends not just on population size but type/intensity of selection, rate and effect of mutations, and number of loci involved (Houle 1989, Keightley and Hill 1987)
	- with directional selection V_A can increase unbounded but stabilizing selection has a limit on the possible V_A for large populations (because directional selection favors alleles in one direction... but this would reduce variation wouldn't it?)
- Maximum heritability under stabilizing selection depends on population size and number of loci
	- selection eliminates genetic variation (but effect is small when many loci are involved - Robertson?)
	- selection also favors genotypes with beneficial alleles at many loci -> linkage equilibrium (which reduces variation by reducing the independence of loci)
	- selection also decreases the number of individuals that reproduce and increases relatedness above them and produces variance in reproductive success, reducing N_e
- Dispersal: movement increases V_A because immigrants counteract losses from drift more rapidly than mutation does
- For longer-term directional selection, drift decreases chances of fixation of rare beneficial alleles and counteracts selection
	- the selection differential also itself decreases under longer-term selection (?)
	- Hill 1982 and Nomura 1997 show that there's more response to selection in populations that are small - variation of response is 2FV_A for inbreeding coefficient F
		- staying smaller for longer increases F and thus increases variation in response
	- think I should also look at Verrier et al. 1991?
- Bottlenecks: reduce V_A due to drift and also influence individual genes (does this depend on frequency of dominant alleles...)

##### Empirical studies
- Swindell and Bouzat 2005: Drosophila heritability of bristle number declines with increasing time at small size, consistent with models of drift+mutation or drift+mutation+stabilizing selection
- Cross-study comparison looking at F (can be derived from generalized study design) and and V_A
	- F gives the expected decline in heterozygosity at neutral loci
	- V_A declines linearly with F as with a neutral model with additive effects
	- Fig. 2... should digest what's going on with fitness traits although I suppose that the number of studies with fitness is very small; but what does it mean if y axis is greater than 1??
- In natural populations, often N_e does not predict quantitative genetic variation
	- Willi et al. 2007 did find smaller populations had lower heritability when averaging over eight traits in 13 natural plant populations
	- Waldmann 2001 and Widen and Andersson 1993 did not find relationships between heritability and population size although they were underpowered
	- a combined meta-analysis of five studies here found no combined evidence of an effect of heritability on population size
	- Reed and Frankham 2001: meta-analysis shows no evidence that quantitative variation and population size (15 studies?)
		- heritability and neutral variation also appear unrelated
	- Could be that small populations only declined only recently and still have variance-increasing effects of dominance or epistatic variance
		- if V_A drops under drift-mutation balance at 1/(2N_e) per generation, substantial time may be needed to see effects of drift
	- traits may be experiencing weak stabilizing selection such that V_A is only sensitive to population sizes at small size (at larger size the variance would asymptote)
	- Also could be that the smallest populations where the effect would be strongest could already have gone extinct
	- Or, gene flow may decouple V_A from population size
- Abundant experimental tests showing populations held at small size for a long time have reduced short- and long-term adaptive response
	- refs herein
	- likely explanation here is drift counteracting selection
		- Weber & Diggins (1990) looked at nine long-term experiments and found reduced selection capacity was stronger than expected due to a drift/selection model
		- could be few loci or linkage disequilibrium (Bulmer effect)
- Bottlenecks: original work by Bryant (1986 and 1993 papers) found no evidence of V_A reductions in bottlenecked populations
	- subsequent studies have both found increases in V_A but also decreases consistent with additive expectations
	- life history traits tend to show the increases in additive variance with inbreeding - possibly due to inbreeding depression and strong selection (purging?)
	- hypothesis: bottlenecks increase V_A by changing nonadditive variance to additive variance (traits associated with fitness, experiencing inbreeding depression, have high dominance variance)
	- then, can bottlenecks increase adaptive potential?
		- bottlenecks typically reduce fitness (inbreeding depression) and surfaces deleterious recessive alleles - increased adaptive ability probably can't overcome this
		- increases in V_A for some traits may be accompanied by decreases in V_A for others

##### Other things reducing population adaptability

- Small populations could be marginalized to less favorable environments, exposing them to less favorable conditions and lowering fitness
	- what is needed to adapt in these environments may be different than in more suitable habitat, e.g., adaptation to environmental variation (plasticity?), adaptation to stronger stress
- Lower individual fitness due to inbreeding depression, genetic load, nongenetic Allee effects
	- Crnokrak and Roff 1999 review inbreeding depression in natural populations
	- Lynch and Lande 1993 look at critical rates of change of a trait conferring maximal fitness
	- Fischer and Mathies 1998 is a good empirical look with plant populations of various sizes - large populations can adapt to rapid environmental change but smaller populations can not
	- Armbruster and Reed 2005 found inbreeding depression on average 69% more severe in stressful environments

Hmm... okay good. Theory shows that most of the time but with some exceptions (not relevant to our model) small populations will lose genetic variation and heritability will fall. Selection seems to also reduce V_A (especially stabilizing selection, which isn't in our model).

Empirical results are mixed. May be due to the fact that multiple traits are involved and/or more complicated genetic architectures as well as the presence of dominance or epistasis. The additive/non-additive variances were slightly confusing to me - I guess this is dominance variance? Means that our results depend on the trait (surprise!)

The other points about marginal habitat and different types of selection are very interesting but not as relevant to us. The point about inbreeding depression is good - populations aren't just threatened by lack of adaptability, but by the reduced fitness! In fact, they reinforce. Which we see in our model. 

Some empirical stuff in here would be good to look at. Willi's other paper look cool. Reviews: Reed and Frankham, Armbruster and Reed. Fischer and Mathies looks super interesting but perhaps not for this project. Theory by Lande as well as Lynch and Lande would be good to read. Need to read S. Wright. 

### Bell, G., and Gonzalez, A. 2009. Evolutionary rescue can prevent extinction following environmental change. Ecology Letters.

- Rescue in response to a monogenic (?) threat - yeast exposed to salt
	- populations will be likely although not certain to adapt if Nf > 1 where f is mutations capable of rescuing populations from stress
- Yeast: haploid, asexually reproducing, susceptible to osmotic and ionic stress reducing growth (individual or population?) in presence of salt
	- 10g/L inhibits growth, 150 g/L is lethal to wild types, experiment run at 125 g/L
	- Expmt 1: can rescue occur in a large population (initial size 6 x 10^5)
	- Expmt 2: can rescue occur in any size?
		- 10 levels of pouplation size, ranging form 10^7 to ~10, with two different growth media
	- Expmt 3: what is the genetic basis of adaptation?
		- look at 14 ORFs sensitive to cation regulation and halosensitivity
		- repeat study with one ORF deleted at a time?
		- (if a strain with ann ORF deleted survives, gene is unlikely to be effective)
- Expmt 1: geometric decline over first 24 hours to ~10% of original size, but returns to original size geometrically (G&H's U shape confirmed)
	- estimating the frequency of resistant cells in the original population: extrapolating the ascending u-shaped recovery curve backwards to the original size (note: must be asexual reproduction for this to work) and also just spread samples of OG population on a salt plate
		- initial frequency fo 2-2.5%
- Using data here to get the initial population size needed to get 50% probability of rescue and also the beneficial mutation rate
- Most of the deletion strains had growth similar to wild type, but one had highly reduced growth and another had eessentially zero growth
- "50% probability of adaptation may require many 100s of individuals even only in the presence of demographic stochasticity"
	- this means more individuals are likely needed in variable environments, or in populations with reduced standing variation
- Also: empirical evidence that size matters (see Fig. 3)

Uses to us: size is important (empirical demonstration), rescue can happen of course, simple study system, contrast the comment about 100s of individuals being needed for rescue to occur

### Vinton, A.C., and Vasseur, D.A. 2020. Evolutionary tracking is determined by differential selection on demographic rates and density dependence. Ecology and Evolution.

- For fast birth, fast death populations, response to selection will be faster (than a slow birth, slow death population with the same intrinsic growth rate) due to more population turnover
- Likewise, stochastic processes can cause deviation from expected outcomes (in - evolution?) esp. for small populations (so this should be modeled explicitly)
- So with adaptation, in addition to knowing $r$, we need to understand both the birth/death parameters and also how r depends on density dependence and environmental change
- Here: how does evolutionary tracking (and thereby rescue?) depend on environmental change and its effects on demographic rates
##### Model
- Continuous time, individual-based logistic growth model
	- environment alters birth or death in a density-independent or density-dependent way (four cases)
	- N' = g(N) * N
		- g(N) is difference between per cap birth rate B and death rate D
		- N'/N = b0 - bI - Nbd - (d0 + dI + NdD)
		- b0, d0 are background birth/death, bI/dI are density independent modifications to them, bD and dD are density-dependent modifications
	- optimum trait value varies sinusoidally, and abs of environmental deviation from it is incorporated into demographic rates
		- equilibrium population size (w/o environmental change) is same in all cases
	- using these, one can solve for parameters that will give equilibrium population size as a function of (environmental change?)
- Stochasticity added through the Gillespie algorithm
	- heritability variation in individual traits incorporated (somehow)
	- random selection on birth and death, also random mutations during reproduction
##### Results
- Persistence was greater in populations with birth rates environmentally influenced and when environment affected density dependence
	- n.b. according to Figure 1 the density dependence means at small population size all genotypes have the same birth/death (think about this...)
	- similarly for death rates persistence is more likely when environment influences density dependence
	- (supposedly has something to do with behavior when the population is small... I suppose a lot of births when small is more helpful than fewer deaths? - interesting question)
	- but also this relationship exists in the no-evolution cases (why does extinction not happen here?)
	- note: Nisbet & Gurney 2003 and Palamara et al. 2016 show that stochasticity is determined by the sum of birth and death rates; somehow the death models have higher stochasticity at equilibrium points (is this due to parameterization only?)
	- ohhhh I see Fig. 1 is important - the vertical line is the other (non-density dependent vital rate) and CC is where they intersect ofc...  
- In all four scenarios, with a greater frequency of environmental variation,   persistence time increases
	- "although the environmental density effect increases variation at high population sizes, it is favorable when populations are small as they are better able to rebound" - counterintuitive... or actually does this match some of the amelioration terms? there's some stuff in here e.g. Ferguson and Ponciano 2015 that look at environmental terms that influence density and show a "multiplicative effect"
	- ecological tracking: ability to track a changing environment diminishes with increasing frequency? stronger tracking at low frequencies means greater variation in population density which means more time spent below mean?
- U-shaped relationship between frequency of env. change and mean persistence times
	- not monotonic; declining at first when going from small to medium rates of change but then increasing persistence going from medium to high (true in all cases)
	- trade-off between ecological and evolutionary tracking; with stronger evolutionary (phenotypic?) tracking, there is less ecological tracking because there's less variation in vital rates (wait... why? is it because the tracking means that populations are less adapted, and therefore large, so there's less demographic stochasticity?)
	- ecological tracking is bad (??) (also what is it?)

Okay... I still don't totally get what's going on but I understand it more now than before. Important things: (1) density dependent *environmental effects* aid persistence rather than hurting and (2) the demographic stochasticity is helpful... somehow

The stochasticity - at some point they discussed it being helpful for populations at low density. Was the idea here that smaller populations can be aided by the environmental fluctuations similarly to the Peniston et al. paper? Not sure. Actually I guess it's not super important!

The use here is also writing about this in the context of the "shape" of density dependence. A density dependent demographic effect that influences all individuals equally vs. one favoring more fit individuals over less fit ones. I do think a lot of the result here is just that at low densities the NDD-effect populations have the same (high) growth rate whereas in the independent case the populations all have lower fitness regardless... ah I wonder if this means that a population that's small can rebound even if it's mostly composed of fit individuals? Actually couldn't this impede selection by stopping purging? Ahhh!!... but also cool!

### Houle, D. 1992. Comparing evolvability and variability of quantitative traits. Genetics.

- "we would like to be able to predict the ability of a population to respond to natural or artificial selection; this I will term 'evolvability'"
- Partitioning trait variation into additive variance and remaining variance
	- heritability h^2 is the ratio of V_A to total phenotypic variance
	- response to selection is R = S h^2 i.e. response is selection differential times heritability (S is difference in population mean before and after selection)
	- life history traits tend to have lower heritabilities (likely because there is such little genetic variation in them)
- "A good measure of evolvability is the expected response to selection."
	- h2 is not a sufficient descriptor (h2 dimensionless, R is not, so info about R is missing)
		- because we are "most often interested in" proportional change so standardizing the response is appropriate and useful - use the population mean (why)
		- also S (strength of selection) is not just a function of the fitness landscape/function but also the selected population (and thus V_P the pheno variance)
	- when X (trait) is fitness, then Fisher's fundamental theorem says (with algebra and substitution) that the response of fitness is V_A / \bar{W}
	- Also: under weak Gaussian selection, there's a formula for the response to selection (Bulmer 1980 has this)
		- with algebra and subs., V_A / \bar{X} is a factor in the proportional response and it is the only factor that is knowable (??)
	- "In both of these important cases, directional and weak stabilizing selection, a variance to mean ratio is the most appropriate measure of evolvability of a trait"
- Some stuff looking at data (not clear... what on earth are these regressions supposed to be telling us)
- "it is not possible to define a single measure appropriate for all circumstances
	- e.g. if w2 is on the order of phenotypic variance or less,  then the analytic comparison fails because V_p won't drop out of some expression
	- is this true for us? we have w2 approx 3.57 and V_p approx between 0.5 and 1, but does it really matter if V_p doesn't drop out of our expression? The argument seems hand-waivy.

Hmm, not what I was expecting. The bulk of this paper is looking at variation (in phenotypes, heritability, and V_a) across traits. Only a small part is about this concept of "evolvability". The thrust of the part on "evolvability" appears to be that V_A is a better measure of "evolvability" (response to selection I guess) than heritability although this wasn't in question for us (I guess a historical dispute). 

Some stuff about Gaussian optimizing selection which is/was good to think about. But also there was some stuff about the strength of selection and the phenotypic variance. We may violate this but it wasn't explained what harm this does (if any!). 

The standardization - does this matter to us? Seems like the standardization matters for comparison across traits - does it matter across populations? Or across time, when the optimum is changing? I really don't have a good sense because (see below)...

Is the concept even really thoroughly explained out here? Seems to just be "response to selection" but maybe broader I guess? No aspect of temporality here. Is that what our "adaptive capacity" entails?



