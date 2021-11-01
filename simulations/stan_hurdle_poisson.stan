// Courtesy of: https://mc-stan.org/docs/2_20/stan-users-guide/zero-inflated-section.html

data {
  // Number of data points
  int<lower=0> N;
  // Response (generation of rescue)
  int<lower=0> y[N];
  // Input data
  //  - Four columns: (1) size (2) diversity (3) NDD (4) genotype
  // real x[N,4];
}
parameters {
  // Probability of a "trivial" immediate rescue
  real<lower=0, upper=1> theta; // [N];
  // Mean log-time to rescue conditioned on not having a trivial rescue
  real<lower=0> lambda; // [N];
  // Parameters for time to rescue (for non-trivial rescue)
  // intercept plus 16 coefficients
  // real beta[N,16];
  // Parameters for probability of trivial rescue
  // real gamma[N,16];
}
model {
  
  for (n in 1:N) {
    //     // Probability of trivial immediate rescue
    // theta[n] = gamma[n,1] + 
    //             gamma[n,2] * x[n,1] + // size
    //             gamma[n,3] * x[n,2] + // diversity
    //             gamma[n,4] * x[n,3] + // NDD
    //             gamma[n,5] * x[n,4] + // genotype;
    //             gamma[n,6] * x[n,1] * x[n,2] + // size x diversity
    //             gamma[n,7] * x[n,1] * x[n,3] + // size x NDD
    //             gamma[n,8] * x[n,1] * x[n,4] + // size x genotype
    //             gamma[n,9] * x[n,2] * x[n,3] + // diversity x NDD
    //             gamma[n,10] * x[n,2] * x[n,4] + // diversity x genotype
    //             gamma[n,11] * x[n,3] * x[n,4] + // NDD x genotype
    //             gamma[n,12] * x[n,1] * x[n,2] + x[n,3] + // size x diversity x NDD
    //             gamma[n,13] * x[n,1] * x[n,2] + x[n,4] + // size x diversity x genotype
    //             gamma[n,14] * x[n,1] * x[n,3] + x[n,4] + // size x NDD x genotype
    //             gamma[n,15] * x[n,2] * x[n,3] + x[n,4] + // diversity x NDD x genotype
    //             gamma[n,16] * x[n,1] * x[n,2] + x[n,3] + x[n,4]; // size x diversity x NDD x genotype 
                
  }
  
  // Loop over data
  for (n in 1:N) {
                
    // If trivial immediate rescue
    if (y[n] == 0)
      // Likelihood function on log scale is probability of trivial rescue
      // plus probability of not immediate rescue times the Poisson pdf at zero
      target += log_sum_exp(bernoulli_lpmf(1 | theta),
                            bernoulli_lpmf(0 | theta)
                              + poisson_lpmf(y[n] | lambda));
    else
    // Probability of rescue time is poisson distributed
      target += bernoulli_lpmf(0 | theta)
                  + poisson_lpmf(y[n] | lambda);
  }
  
}