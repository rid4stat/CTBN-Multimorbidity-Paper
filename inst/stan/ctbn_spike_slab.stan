// =============================================================================
// ctbn_spike_slab.stan
// Continuous spike-and-slab CTBN model for multimorbidity.
//
// Embeds the Bayesian CTBN log-linear Poisson likelihood with the
// order-dependent continuous spike-and-slab prior of Section 3.3 of the paper.
// Posterior inclusion probabilities (PIPs) are returned as generated quantities.
//
// Inputs (data block):
//   N_rows    : number of at-risk intervals (panel rows)
//   Q         : number of network coefficients (main + interaction)
//   P         : number of covariate coefficients
//   N_obs[i]  : event count in interval i (0 or 1 for absorbing onset)
//   X         : (N_rows x Q) design matrix of indicator products
//   Z         : (N_rows x P) covariate design matrix
//   T         : (N_rows)     exposure widths Δ_{i,s}
//   beta_order[j] in {0,1,...,max_order} : interaction order of β_j
//   pi0       : base inclusion probability (slab)
//   theta     : order penalty
//   a0, b0    : Inv-Gamma hyperparameters for the slab variance
//   spike_var : variance of the spike component (e.g. 1e-4)
//   max_order : largest interaction order P
// =============================================================================
data {
  int<lower=1> N_rows;
  int<lower=1> Q;
  int<lower=1> P;
  array[N_rows] int<lower=0> N_obs;
  matrix[N_rows, Q] X;
  matrix[N_rows, P] Z;
  vector<lower=0>[N_rows] T;
  array[Q] int<lower=0> beta_order;
  real<lower=0, upper=1> pi0;
  real<lower=0>          theta;
  real<lower=0>          a0;
  real<lower=0>          b0;
  real<lower=0>          spike_var;
  int<lower=1>           max_order;
}
transformed data {
  vector[Q] pi_beta;
  for (j in 1:Q)
    pi_beta[j] = pi0 * exp(-theta * beta_order[j]);
  real pi_gamma = pi0;
}
parameters {
  vector[Q]                       beta;
  vector[P]                       gamma;
  vector<lower=0>[max_order + 1]  sigma2;
  real<lower=0>                   sigma2_gamma;
}
model {
  for (p in 0:max_order)
    sigma2[p + 1] ~ inv_gamma(a0, b0);
  sigma2_gamma ~ inv_gamma(a0, b0);

  for (j in 1:Q)
    target += log_mix(pi_beta[j],
                      normal_lpdf(beta[j]  | 0, sqrt(sigma2[beta_order[j] + 1])),
                      normal_lpdf(beta[j]  | 0, sqrt(spike_var)));

  for (k in 1:P)
    target += log_mix(pi_gamma,
                      normal_lpdf(gamma[k] | 0, sqrt(sigma2_gamma)),
                      normal_lpdf(gamma[k] | 0, sqrt(spike_var)));

  vector[N_rows] log_mu = X * beta + Z * gamma + log(T);
  target += poisson_log_lpmf(N_obs | log_mu);
}
generated quantities {
  vector[Q] pip_beta;
  vector[P] pip_gamma;
  for (j in 1:Q) {
    real ls = log(pi_beta[j])   + normal_lpdf(beta[j]  | 0, sqrt(sigma2[beta_order[j] + 1]));
    real lp = log1m(pi_beta[j]) + normal_lpdf(beta[j]  | 0, sqrt(spike_var));
    pip_beta[j] = exp(ls - log_sum_exp(ls, lp));
  }
  for (k in 1:P) {
    real ls = log(pi_gamma)   + normal_lpdf(gamma[k] | 0, sqrt(sigma2_gamma));
    real lp = log1m(pi_gamma) + normal_lpdf(gamma[k] | 0, sqrt(spike_var));
    pip_gamma[k] = exp(ls - log_sum_exp(ls, lp));
  }
}
