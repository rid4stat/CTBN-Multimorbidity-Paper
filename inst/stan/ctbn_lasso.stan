// =============================================================================
// ctbn_lasso.stan
// Bayesian LASSO CTBN model (Park & Casella 2008) with order-penalised scale.
//
//   β_j | λ²_{p_j} ~ DoubleExp(0, √λ²_{p_j} · exp(-θ p_j / 2))
//   λ²_p          ~ Inv-Gamma(a0, b0)
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
  real<lower=0> theta;
  real<lower=0> a0;
  real<lower=0> b0;
  int<lower=1>  max_order;
}
parameters {
  vector[Q]                      beta;
  vector[P]                      gamma;
  vector<lower=0>[max_order + 1] lambda2;
  real<lower=0>                  lambda2_gamma;
}
model {
  for (p in 0:max_order)
    lambda2[p + 1] ~ inv_gamma(a0, b0);
  lambda2_gamma ~ inv_gamma(a0, b0);

  for (j in 1:Q) {
    real scale_j = sqrt(lambda2[beta_order[j] + 1]) * exp(-theta * beta_order[j] / 2.0);
    beta[j] ~ double_exponential(0, scale_j);
  }

  gamma ~ double_exponential(0, sqrt(lambda2_gamma));

  vector[N_rows] log_mu = X * beta + Z * gamma + log(T);
  target += poisson_log_lpmf(N_obs | log_mu);
}
generated quantities {
  vector[Q] kappa;
  for (j in 1:Q) {
    real lam2_eff = lambda2[beta_order[j] + 1] * exp(-theta * beta_order[j]);
    kappa[j] = 1.0 / (1.0 + 2.0 * lam2_eff);
  }
}
