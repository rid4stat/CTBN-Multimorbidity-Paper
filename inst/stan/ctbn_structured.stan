// =============================================================================
// ctbn_structured.stan
// Hierarchical Structured (HS) prior CTBN model.
//
// Order-penalised normal slab variance:
//   β_j | σ²_{p_j} ~ N(0, σ²_{p_j} · exp(-θ · p_j))
//   σ²_p          ~ Inv-Gamma(a0, b0)
//
// This is a smooth continuous shrinkage prior — no hard selection.
// Use the posterior shrinkage factor as a pseudo-PIP for downstream selection.
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
  vector<lower=0>[max_order + 1] sigma2;
  real<lower=0>                  sigma2_gamma;
}
model {
  for (p in 0:max_order)
    sigma2[p + 1] ~ inv_gamma(a0, b0);
  sigma2_gamma ~ inv_gamma(a0, b0);

  for (j in 1:Q) {
    real eff_var = sigma2[beta_order[j] + 1] * exp(-theta * beta_order[j]);
    beta[j] ~ normal(0, sqrt(eff_var));
  }

  gamma ~ normal(0, sqrt(sigma2_gamma));

  vector[N_rows] log_mu = X * beta + Z * gamma + log(T);
  target += poisson_log_lpmf(N_obs | log_mu);
}
generated quantities {
  vector[Q] post_sd_beta;
  for (j in 1:Q)
    post_sd_beta[j] = sqrt(sigma2[beta_order[j] + 1] * exp(-theta * beta_order[j]));
}
