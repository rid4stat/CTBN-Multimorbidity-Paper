// =============================================================================
// ctbn_horseshoe.stan
// Regularised Horseshoe CTBN model (Piironen & Vehtari 2017) with
// order-penalised global scale.
//
//   β_j = z_j · √(λ̃_j²) · τ_j,   τ_j = τ · exp(-θ · p_j)
//   λ̃_j² = c² λ_j² / (c² + τ_j² λ_j²)        (slab-regularised)
//   z_j ~ N(0,1),  λ_j ~ HalfCauchy(0,1),  τ ~ HalfCauchy(0, τ_0)
//   c² ~ InvGamma(slab_df/2, slab_df · slab_scale²/2)
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
  real<lower=0> tau0;
  real<lower=0> slab_df;
  real<lower=0> slab_scale;
  int<lower=1>  max_order;
}
parameters {
  vector[Q]          z_beta;
  vector[P]          gamma;
  vector<lower=0>[Q] lambda;
  real<lower=0>      tau;
  real<lower=0>      c2;
}
transformed parameters {
  vector[Q] beta;
  for (j in 1:Q) {
    real tau_j     = tau * exp(-theta * beta_order[j]);
    real lambda2_j = square(lambda[j]);
    real lambda2_t = c2 * lambda2_j / (c2 + square(tau_j) * lambda2_j);
    beta[j]        = z_beta[j] * sqrt(lambda2_t) * tau_j;
  }
}
model {
  tau    ~ cauchy(0, tau0);
  c2     ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df * square(slab_scale));
  lambda ~ cauchy(0, 1);
  z_beta ~ std_normal();
  gamma  ~ normal(0, 1);

  vector[N_rows] log_mu = X * beta + Z * gamma + log(T);
  target += poisson_log_lpmf(N_obs | log_mu);
}
generated quantities {
  vector[Q] kappa;
  for (j in 1:Q) {
    real tau_j     = tau * exp(-theta * beta_order[j]);
    real lambda2_j = square(lambda[j]);
    real lambda2_t = c2 * lambda2_j / (c2 + square(tau_j) * lambda2_j);
    kappa[j] = 1.0 / (1.0 + N_rows * square(tau_j) * lambda2_t);
  }
}
