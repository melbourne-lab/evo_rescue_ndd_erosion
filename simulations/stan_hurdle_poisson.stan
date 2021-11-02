// Courtesy of: https://mc-stan.org/docs/2_20/stan-users-guide/zero-inflated-section.html

data {
  // Number of data points
  int<lower=0> N;
  // Response (generation of rescue)
  int<lower=0> y[N];
  // Input data
  //  - Four columns: (1) size (2) diversity (3) NDD (4) genotype
  // real x[N,4];
  int<lower=0, upper=1> x[N];
}
parameters {
  // Probability of a "trivial" immediate rescue
  // real<lower=0, upper=1> theta; // [N];
  // Mean log-time to rescue conditioned on not having a trivial rescue
  // real<lower=0> lambda; // [N];
  
  // Parameters for time to rescue (for non-trivial rescue)
  // intercept plus 16 coefficients
  // real beta[16]; 
  real beta[2];
  // Parameters for probability of trivial rescue
  // real gamma[16];
  real gamma[2];
  
  // real<lower=0> sigma_theta[2];
  // real<lower=0> sigma_lambda[2];
}
model {
  
  real theta[N];
  
  real lambda[N];
  
  for (n in 1:N) {
    
    // (logit-)Probability of trivial immediate rescue
    // theta[n] =  gamma[1] +
    //             gamma[2] * x[n,1] + // size
    //             gamma[3] * x[n,2] + // diversity
    //             gamma[4] * x[n,3] + // NDD
    //             gamma[5] * x[n,4] + // genotype;
    //             gamma[6] * x[n,1] * x[n,2] + // size x diversity
    //             gamma[7] * x[n,1] * x[n,3] + // size x NDD
    //             gamma[8] * x[n,1] * x[n,4] + // size x genotype
    //             gamma[9] * x[n,2] * x[n,3] + // diversity x NDD
    //             gamma[10] * x[n,2] * x[n,4] + // diversity x genotype
    //             gamma[11] * x[n,3] * x[n,4] + // NDD x genotype
    //             gamma[12] * x[n,1] * x[n,2] * x[n,3] + // size x diversity x NDD
    //             gamma[13] * x[n,1] * x[n,2] * x[n,4] + // size x diversity x genotype
    //             gamma[14] * x[n,1] * x[n,3] * x[n,4] + // size x NDD x genotype
    //             gamma[15] * x[n,2] * x[n,3] * x[n,4] + // diversity x NDD x genotype
    //             gamma[16] * x[n,1] * x[n,2] * x[n,3] * x[n,4]; // size x diversity x NDD x genotype
    
    theta[n] = gamma[1] + gamma[2]*x[n];
                // sigma_theta);
   
   // (log-)Rate for rescue time given non-trivial rescue time
   // lambda[n] =  beta[1] +
   //              beta[2] * x[n,1] + // size
   //              beta[3] * x[n,2] + // diversity
   //              beta[4] * x[n,3] + // NDD
   //              beta[5] * x[n,4] + // genotype;
   //              beta[6] * x[n,1] * x[n,2] + // size x diversity
   //              beta[7] * x[n,1] * x[n,3] + // size x NDD
   //              beta[8] * x[n,1] * x[n,4] + // size x genotype
   //              beta[9] * x[n,2] * x[n,3] + // diversity x NDD
   //              beta[10] * x[n,2] * x[n,4] + // diversity x genotype
   //              beta[11] * x[n,3] * x[n,4] + // NDD x genotype
   //              beta[12] * x[n,1] * x[n,2] * x[n,3] + // size x diversity x NDD
   //              beta[13] * x[n,1] * x[n,2] * x[n,4] + // size x diversity x genotype
   //              beta[14] * x[n,1] * x[n,3] * x[n,4] + // size x NDD x genotype
   //              beta[15] * x[n,2] * x[n,3] * x[n,4] + // diversity x NDD x genotype
   //              beta[16] * x[n,1] * x[n,2] * x[n,3] * x[n,4]; // size x diversity x NDD x genotype
   
   lambda[n] = beta[1] + beta[2]*x[n];
                // sigma_lambda);
                
  }
  
  // Loop over data
  for (n in 1:N) {
                
    // If trivial immediate rescue
    if (y[n] == 0)
      // Likelihood function on log scale is probability of trivial rescue
      // plus probability of not immediate reåscue times the Poisson pdf at zero
      target += log_sum_exp(bernoulli_logit_lpmf(1 | theta[n]),
                            bernoulli_logit_lpmf(0 | theta[n])
                              + poisson_log_lpmf(y[n] | lambda[n]));
    else
    // Probability of rescue time is poisson distributed
      target += bernoulli_logit_lpmf(0 | theta[n])
                  + poisson_log_lpmf(y[n] | lambda[n]);
  }
  
}
