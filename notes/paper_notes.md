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

