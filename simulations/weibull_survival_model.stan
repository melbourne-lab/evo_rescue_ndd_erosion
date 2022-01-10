/*  Variable naming:
 obs       = observed
 cen       = (right) censored
 N         = number of samples
 M         = number of covariates
 bg        = established risk (or protective) factors
 tau       = scale parameter
*/
// Tomi Peltola, tomi.peltola@aalto.fi
// Downloaded by SN - 10 Jan 2022
// biostan repo: https://github.com/jburos/biostan/blob/master/inst/stan/weibull_survival_model.stan

functions {
  vector sqrt_vec(vector x) {
    vector[dims(x)[1]] res;

    for (m in 1:dims(x)[1]){
      res[m] = sqrt(x[m]);
    }

    return res;
  }

  vector bg_prior_lp(real r_global, vector r_local) {
    r_global ~ normal(0.0, 10.0);
    r_local ~ inv_chi_square(1.0);

    return r_global * sqrt_vec(r_local);
  }
}

data {
  int<lower=0> Nobs;
  int<lower=0> Ncen;
  int<lower=0> M_bg;
  vector[Nobs] yobs;
  vector[Ncen] ycen;
  matrix[Nobs, M_bg] Xobs_bg;
  matrix[Ncen, M_bg] Xcen_bg;
}

transformed data {
  real<lower=0> tau_mu;
  real<lower=0> tau_al;

  tau_mu = 10.0;
  tau_al = 10.0;
}

parameters {
  real<lower=0> tau_s_bg_raw;
  vector<lower=0>[M_bg] tau_bg_raw;

  real alpha_raw;
  vector[M_bg] beta_bg_raw;

  real mu;
}

transformed parameters {
  vector[M_bg] beta_bg;
  real alpha;

  beta_bg = bg_prior_lp(tau_s_bg_raw, tau_bg_raw) .* beta_bg_raw;
  alpha = exp(tau_al * alpha_raw);
}

model {
  yobs ~ weibull(alpha, exp(-(mu + Xobs_bg * beta_bg)/alpha));
  target += weibull_lccdf(ycen | alpha, exp(-(mu + Xcen_bg * beta_bg)/alpha));

  beta_bg_raw ~ normal(0.0, 1.0);
  alpha_raw ~ normal(0.0, 1.0);

  mu ~ normal(0.0, tau_mu);
}

/*
generated quantities {
    real yhat_uncens[Nobs + Ncen];
    real log_lik[Nobs + Ncen];
    real lp[Nobs + Ncen];

    for (i in 1:Nobs) {
        lp[i] = mu + Xobs_bg[i,] * beta_bg;
        yhat_uncens[i] = weibull_rng(alpha, exp(-(mu + Xobs_bg[i,] * beta_bg)/alpha));
        log_lik[i] = weibull_lpdf(yobs[i] | alpha, exp(-(mu + Xobs_bg[i,] * beta_bg)/alpha));
    }
    for (i in 1:Ncen) {
        lp[Nobs + i] = mu + Xcen_bg[i,] * beta_bg;
        yhat_uncens[Nobs + i] = weibull_rng(alpha, exp(-(mu + Xcen_bg[i,] * beta_bg)/alpha));
        log_lik[Nobs + i] = weibull_lccdf(ycen[i] | alpha, exp(-(mu + Xcen_bg[i,] * beta_bg)/alpha));
    }
}
*/
