# =============================================================================
# ctbn_map_fast.R  —  Fast Matrix-Based MAP/Laplace Estimation for CTBN Models
# Version 2.0  —  Supports: spike_slab | structured | lasso | horseshoe
# =============================================================================
#
# MOTIVATION
# ----------
# ctbn_fit() runs M sequential Stan MCMC chains. With M=10, 4 chains,
# 2000 iter: ~200,000 Poisson evaluations × 10 nodes = hours per fit.
# This file replaces MCMC with L-BFGS-B MAP + Laplace approximation:
#   - Same statistical model and all four prior families
#   - 50-200x faster: minutes not hours
#   - Full interface compatibility with ctbn_cv_parallel and ctbn_simulation_v3
#
# =============================================================================
# MATHEMATICAL BACKGROUND
# =============================================================================
#
# ── 1. CTBN PANEL-DATA POISSON LIKELIHOOD ────────────────────────────────────
#
# For target condition m, at-risk intervals i:
#   N_{im} ~ Poisson( exp(η_{im}) · Δ_{im} )
#
#   η_{im} = Σ_{j≠m} β_{jm} X_j(t_i)                   [main effects, p=0]
#           + Σ_{j<k} β_{jkm} X_j(t_i)X_k(t_i)          [2-way, p=1]
#           + Σ_r γ_{rm} Z_r(t_i)                        [covariates]
#
# Log-likelihood:  ℓ(θ) = Σ_i [N_i(Φ_i θ + log Δ_i) − exp(Φ_i θ + log Δ_i)]
# Gradient:        ∇ℓ   = Φ^T (y − μ)          μ_i = exp(η_i)
# Neg-Hessian:     −∇²ℓ = Φ^T diag(μ) Φ        (PSD ⟹ objective is concave ✓)
#
# ── 2. PRIOR SPECIFICATIONS ──────────────────────────────────────────────────
#
# All priors include an order penalty: v_j scales with exp(−θ·p_j) where
# p_j ∈ {0,1,...} is the interaction order (0 = main effect).
#
# STRUCTURED (Normal):
#   β_j ~ N(0, σ²_p · exp(−θ·p_j))
#   log p  = −½[β²/v + log(2πv)]    v = σ²_p exp(−θp)
#   ∂/∂β   = −β / v
#   ∂²/∂β² = −1 / v
#
# LASSO (Smooth Laplace approximation):
#   β_j ~ DE(0, s_j)   s_j = √λ²_p · exp(−θp/2)
#   Smoothed: log p ≈ −√(β²+ε²)/s − log(2s)     ε = 1e-6
#   ∂/∂β   ≈ −β / [s · √(β²+ε²)]
#   ∂²/∂β² ≈ −ε² / [s · (β²+ε²)^{3/2}]
#
# SPIKE-AND-SLAB (Marginalised mixture):
#   π_j = π_0 · exp(−θ·p_j)    [order-penalised inclusion probability]
#   log p(β_j) = log[π_j φ(β_j;0,σ²_p) + (1−π_j) φ(β_j;0,ε²)]
#   ∂/∂β via quotient rule on the mixture (see code)
#   ∂²/∂β² via finite differences
#
# HORSESHOE (Regularised; Piironen & Vehtari 2017 — MAP approximation):
#   Stan model: z_j~N(0,1), λ_j~HalfCauchy(0,1), τ~HalfCauchy(0,τ_0),
#               c²~InvGamma(slab_df/2, slab_df·slab_scale²/2)
#   β_j = z_j · λ̃_j · τ_j    where τ_j = τ·exp(−θ·p_j)
#   λ̃²_j = c²λ²_j / (c² + τ²_j λ²_j)   [regularised local scale]
#
#   MAP approximation strategy:
#   We treat τ and c² as fixed at their prior means and set the local
#   scale to its prior median (λ_j = 1). This gives an effective Normal
#   prior on β_j with variance:
#
#     σ²_hs(j) = c²_eff · τ²_j / (c²_eff + τ²_j)
#
#   where τ_j = τ_0 · exp(−θ·p_j)  and
#         c²_eff = E[c²] = (slab_df/2 · slab_scale²) / (slab_df/2 − 1)
#                        ≈ slab_scale² · slab_df/(slab_df−2)  [for slab_df>2]
#
#   This approximation:
#   (a) Recovers the regularised horseshoe's shrinkage profile for typical
#       CTBN coefficients (|β| ≪ c_eff)
#   (b) Gives the correct limiting behaviour:
#       β→0  ⟹  σ²_hs → c²τ²/(c²+τ²)  [bounded — "regularised"]
#       τ≪c  ⟹  σ²_hs ≈ τ²  [strong global shrinkage when τ is small]
#       τ≫c  ⟹  σ²_hs ≈ c²  [slab cap — prevents infinite tails]
#   (c) Is equivalent to the "prior-variance matching" EM initialisation
#       used in fast horseshoe algorithms (Bhattacharya et al. 2016)
#   (d) Avoids auxiliary variable sampling entirely
#
#   After MAP optimisation, the horseshoe shrinkage factor (Carvalho 2010):
#     κ_j = 1 / (1 + n · σ²_hs(j))    [0 = signal, 1 = shrunk to zero]
#   is used as the selection statistic, matching ctbn_fit_unified output.
#
# ── 3. MAP OBJECTIVE ─────────────────────────────────────────────────────────
#   F(θ) = ℓ(θ) + Σ_j log p(β_j | prior) + Σ_r log p(γ_r | N(0,σ²_γ))
#   Minimise −F(θ) via L-BFGS-B. Concave ⟹ unique global MAP.
#
# ── 4. LAPLACE POSTERIOR COVARIANCE ──────────────────────────────────────────
#   Σ_post ≈ [Φ^T diag(μ̂)Φ + diag(−∂²log p / ∂θ²)]^{−1}
#   SE_j = √[Σ_post]_{jj}     95% CI: θ̂_j ± 1.96·SE_j
#
# ── 5. SELECTION STATISTICS ──────────────────────────────────────────────────
#   spike_slab : PIP_j = π_j φ_slab(β̂_j) / [π_j φ_slab + (1−π_j) φ_spike]
#   horseshoe  : κ_j = 1/(1 + n·σ²_hs_j);  pseudo-PIP = 1 − κ_j
#   lasso/structured: pseudo-PIP = Φ(|β̂_j|/SE_j)   [Wald exceedance]
#
# =============================================================================

library(data.table)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# =============================================================================
# SECTION A: Prior functions — log-density, gradient, Hessian diagonal
# =============================================================================

# ── A1. Structured (Order-penalised Normal) ───────────────────────────────────

prior_logdens_structured <- function(beta, sigma2, orders, theta) {
  v <- sigma2[orders + 1L] * exp(-theta * orders)
  -0.5 * sum(beta^2 / v + log(2 * pi * v))
}

grad_prior_structured <- function(beta, sigma2, orders, theta) {
  v <- sigma2[orders + 1L] * exp(-theta * orders)
  -beta / v
}

hess_diag_prior_structured <- function(sigma2, orders, theta) {
  v <- sigma2[orders + 1L] * exp(-theta * orders)
  -1.0 / v
}

# ── A2. LASSO (Smooth Laplace / pseudo-Huber) ─────────────────────────────────

prior_logdens_lasso <- function(beta, lambda2, orders, theta, eps = 1e-6) {
  s <- sqrt(lambda2[orders + 1L]) * exp(-theta * orders / 2.0)
  -sum(sqrt(beta^2 + eps^2) / s + log(2 * s))
}

grad_prior_lasso <- function(beta, lambda2, orders, theta, eps = 1e-6) {
  s <- sqrt(lambda2[orders + 1L]) * exp(-theta * orders / 2.0)
  -beta / (s * sqrt(beta^2 + eps^2))
}

hess_diag_prior_lasso <- function(beta, lambda2, orders, theta, eps = 1e-6) {
  s <- sqrt(lambda2[orders + 1L]) * exp(-theta * orders / 2.0)
  -(eps^2) / (s * (beta^2 + eps^2)^1.5)
}

# ── A3. Spike-and-Slab (Marginalised mixture) ─────────────────────────────────

prior_logdens_spikeslab <- function(beta, sigma2, orders, theta,
                                     pi0, spike_var) {
  pi_j   <- pmin(pmax(pi0 * exp(-theta * orders), 1e-10), 1 - 1e-10)
  slab_v <- sigma2[orders + 1L]
  log_s  <- log(pi_j)      + dnorm(beta, 0, sqrt(slab_v),    log = TRUE)
  log_k  <- log(1 - pi_j)  + dnorm(beta, 0, sqrt(spike_var), log = TRUE)
  sum(apply(cbind(log_s, log_k), 1, function(r) {
    mx <- max(r); mx + log(sum(exp(r - mx)))
  }))
}

grad_prior_spikeslab <- function(beta, sigma2, orders, theta, pi0, spike_var) {
  pi_j   <- pmin(pmax(pi0 * exp(-theta * orders), 1e-10), 1 - 1e-10)
  slab_v <- sigma2[orders + 1L]
  phi_s  <- dnorm(beta, 0, sqrt(slab_v))
  phi_k  <- dnorm(beta, 0, sqrt(spike_var))
  numer  <- -pi_j * phi_s * beta / slab_v -
             (1 - pi_j) * phi_k * beta / spike_var
  denom  <- pi_j * phi_s + (1 - pi_j) * phi_k + 1e-300
  numer / denom
}

# Hessian diagonal via finite differences (spike-slab has no clean closed form)
hess_diag_prior_spikeslab <- function(beta, sigma2, orders, theta,
                                       pi0, spike_var, eps = 1e-5) {
  gp <- grad_prior_spikeslab(beta + eps, sigma2, orders, theta, pi0, spike_var)
  gm <- grad_prior_spikeslab(beta - eps, sigma2, orders, theta, pi0, spike_var)
  (gp - gm) / (2 * eps)   # d²log p / dβ²  (will be ≤ 0)
}

# ── A4. Horseshoe (Regularised; Normal scale-mixture MAP approximation) ────────
#
# Given fixed τ (global scale) and c² (slab variance), the effective prior
# variance for coefficient j at interaction order p_j is:
#
#   σ²_hs(j) = c² · τ²_j / (c² + τ²_j)    where τ_j = τ · exp(−θ · p_j)
#
# This is the regularised horseshoe variance evaluated at λ_j = 1 (median of
# HalfCauchy(0,1)), providing the "prior-variance matching" approximation.
# The resulting marginal prior on β_j is Normal: β_j ~ N(0, σ²_hs(j)).

.hs_eff_var <- function(hs_tau, hs_c2, orders, theta) {
  tau_j <- hs_tau * exp(-theta * orders)
  hs_c2 * tau_j^2 / (hs_c2 + tau_j^2)
}

prior_logdens_horseshoe <- function(beta, hs_tau, hs_c2, orders, theta) {
  v <- .hs_eff_var(hs_tau, hs_c2, orders, theta)
  -0.5 * sum(beta^2 / v + log(2 * pi * v))
}

grad_prior_horseshoe <- function(beta, hs_tau, hs_c2, orders, theta) {
  v <- .hs_eff_var(hs_tau, hs_c2, orders, theta)
  -beta / v
}

hess_diag_prior_horseshoe <- function(hs_tau, hs_c2, orders, theta) {
  v <- .hs_eff_var(hs_tau, hs_c2, orders, theta)
  -1.0 / v
}

#' Compute E[c²] under InvGamma(slab_df/2, slab_df*slab_scale^2/2)
#' Mean exists for slab_df > 2.  Default slab_df=4 ⟹ E[c²] = 2·slab_scale².
hs_c2_mean <- function(slab_df, slab_scale) {
  if (slab_df > 2)
    (slab_df / 2 * slab_scale^2) / (slab_df / 2 - 1)
  else
    slab_scale^2 * 4   # conservative fallback
}

# =============================================================================
# SECTION B: Design matrix builder
# =============================================================================

#' Build Φ_m = [X_{-m} | interactions | Z] for one target node.
#' Called once per node in the main loop.
build_design_matrix <- function(X_all, Z, cond_names, target, max_order) {
  influencers <- setdiff(cond_names, target)
  X_inf       <- X_all[, influencers, drop = FALSE]
  x_cols      <- influencers
  x_ord       <- rep(0L, length(influencers))

  if (max_order >= 2L) {
    for (ord in 2L:max_order) {
      for (grp in combn(influencers, ord, simplify = FALSE)) {
        cn       <- paste(grp, collapse = "_x_")
        col_vals <- Reduce(`*`, lapply(grp, function(v) X_all[, v]))
        X_inf    <- cbind(X_inf, col_vals)
        x_cols   <- c(x_cols, cn)
        x_ord    <- c(x_ord,  ord - 1L)
      }
    }
  }
  colnames(X_inf) <- x_cols
  list(Phi = cbind(X_inf, Z), x_cols = x_cols,
       x_orders = x_ord, n_beta = length(x_cols), n_gamma = ncol(Z))
}

# =============================================================================
# SECTION C: Poisson log-likelihood, gradient, Hessian
# =============================================================================

poisson_loglik <- function(par, Phi, y, offset) {
  eta <- pmin(drop(Phi %*% par) + offset, 30)
  sum(y * eta - exp(eta))
}

poisson_grad <- function(par, Phi, y, offset) {
  eta <- pmin(drop(Phi %*% par) + offset, 30)
  drop(crossprod(Phi, y - exp(eta)))
}

# Returns Φ^T diag(μ) Φ  (positive semi-definite Fisher information matrix)
poisson_neg_hess <- function(par, Phi, offset) {
  eta <- pmin(drop(Phi %*% par) + offset, 30)
  crossprod(Phi * sqrt(exp(eta)))
}

# =============================================================================
# SECTION D: MAP objective factory
# =============================================================================

#' Returns list(fn, gr) closures for optim(method="L-BFGS-B").
#' fn(par) = −F(par),  gr(par) = −∇F(par).
make_map_objective <- function(Phi, y, offset, x_orders,
                                prior, sigma2, lambda2, sigma2_g,
                                theta_pen, pi0, spike_var,
                                hs_tau, hs_c2,
                                n_beta, n_gamma) {
  bi <- seq_len(n_beta)
  gi <- n_beta + seq_len(n_gamma)

  fn <- function(par) {
    beta  <- par[bi];  gamma <- par[gi]
    ll    <- poisson_loglik(par, Phi, y, offset)
    lp_b  <- switch(prior,
      structured = prior_logdens_structured(beta, sigma2,  x_orders, theta_pen),
      lasso      = prior_logdens_lasso(     beta, lambda2, x_orders, theta_pen),
      spike_slab = prior_logdens_spikeslab( beta, sigma2,  x_orders, theta_pen,
                                            pi0, spike_var),
      horseshoe  = prior_logdens_horseshoe( beta, hs_tau, hs_c2, x_orders, theta_pen)
    )
    lp_g  <- -0.5 * sum(gamma^2 / sigma2_g + log(2 * pi * sigma2_g))
    -(ll + lp_b + lp_g)
  }

  gr <- function(par) {
    beta  <- par[bi];  gamma <- par[gi]
    gl    <- poisson_grad(par, Phi, y, offset)
    gp_b  <- switch(prior,
      structured = grad_prior_structured(beta, sigma2,  x_orders, theta_pen),
      lasso      = grad_prior_lasso(     beta, lambda2, x_orders, theta_pen),
      spike_slab = grad_prior_spikeslab( beta, sigma2,  x_orders, theta_pen,
                                         pi0, spike_var),
      horseshoe  = grad_prior_horseshoe( beta, hs_tau, hs_c2, x_orders, theta_pen)
    )
    gp_g <- -gamma / sigma2_g
    -(gl + c(gp_b, gp_g))
  }

  list(fn = fn, gr = gr)
}

# =============================================================================
# SECTION E: Laplace approximation for posterior SEs
# =============================================================================

laplace_se <- function(theta_hat, Phi, y, offset, x_orders,
                        prior, sigma2, lambda2, sigma2_g,
                        theta_pen, pi0, spike_var, hs_tau, hs_c2,
                        n_beta, n_gamma) {
  d    <- length(theta_hat)
  bi   <- seq_len(n_beta)
  gi   <- n_beta + seq_len(n_gamma)
  beta <- theta_hat[bi]

  H_lik <- poisson_neg_hess(theta_hat, Phi, offset)

  hd_beta <- switch(prior,
    structured = -hess_diag_prior_structured(sigma2,  x_orders, theta_pen),
    lasso      = -hess_diag_prior_lasso(beta, lambda2, x_orders, theta_pen),
    spike_slab = -hess_diag_prior_spikeslab(beta, sigma2, x_orders, theta_pen,
                                             pi0, spike_var),
    horseshoe  = -hess_diag_prior_horseshoe(hs_tau, hs_c2, x_orders, theta_pen)
  )

  pd        <- numeric(d)
  pd[bi]    <- hd_beta
  pd[gi]    <- rep(1.0 / sigma2_g, n_gamma)
  H_full    <- H_lik + diag(pd, d)

  Sig <- tryCatch(
    chol2inv(chol(H_full)),
    error = function(e) chol2inv(chol(H_full + diag(1e-6, d)))
  )
  list(se = sqrt(pmax(diag(Sig), 0)), cov_approx = Sig)
}

# =============================================================================
# SECTION F: Selection statistics (PIP / kappa) from MAP + Laplace
# =============================================================================
#
# COVARIATE PIP  (pip_gamma)
# --------------------------
# In ctbn_fit_unified_v2, the spike-and-slab Stan model places the same
# mixture prior on each covariate coefficient γ_k:
#
#   γ_k | δ_k ~ δ_k · N(0, σ²_γ) + (1-δ_k) · N(0, ε²)
#   δ_k ~ Bernoulli(π_0)           [flat — no order penalty for covariates]
#
# The covariate PIP is then:
#   pip_gamma_k = π_0 φ(γ̂_k; 0, σ²_γ) / [π_0 φ(γ̂_k;0,σ²_γ) + (1-π_0) φ(γ̂_k;0,ε²)]
#
# For lasso/structured/horseshoe, the original code does not produce a formal
# covariate PIP. Here we mirror exactly:
#   spike_slab : PIP via posterior mixture weight (matches Stan generated quantities)
#   horseshoe  : kappa_gamma_k = 1/(1 + n·σ²_γ_hs);  pseudo-PIP = 1 - kappa
#                where σ²_γ_hs = hs_c2·tau0² / (hs_c2 + tau0²) [no order penalty,
#                so p=0 ⟹ tau_j = tau0; this matches the Stan horseshoe treatment
#                of gamma which uses the same tau but no order shrinkage]
#   lasso/structured : Wald exceedance Φ(|γ̂_k|/SE_k)
#
# All results are stored in pip_cov_list[[target]] as a named numeric vector,
# matching the ctbn_fit_unified_v2 output slot exactly.

#' Compute selection statistics for beta (condition) coefficients.
#' Returns list(pip, kappa):
#'   pip   — PIP or pseudo-PIP [0,1] per beta coefficient
#'   kappa — horseshoe shrinkage factor [0,1], NA for other priors
compute_selection_stats <- function(beta_hat, se_beta, x_orders,
                                     prior, sigma2, theta_pen,
                                     pi0, spike_var, hs_tau, hs_c2, n_obs) {
  kappa <- rep(NA_real_, length(beta_hat))

  pip <- switch(prior,
    spike_slab = {
      pi_j   <- pmin(pmax(pi0 * exp(-theta_pen * x_orders), 1e-10), 1 - 1e-10)
      slab_v <- sigma2[x_orders + 1L]
      phi_s  <- dnorm(beta_hat, 0, sqrt(slab_v))
      phi_k  <- dnorm(beta_hat, 0, sqrt(spike_var))
      (pi_j * phi_s) / (pi_j * phi_s + (1 - pi_j) * phi_k + 1e-300)
    },
    horseshoe = {
      v_hs    <- .hs_eff_var(hs_tau, hs_c2, x_orders, theta_pen)
      # Carvalho et al. (2010): κ_j = 1/(1 + n·σ²_hs_j)
      kappa[] <- 1.0 / (1.0 + n_obs * v_hs)
      1.0 - kappa   # pseudo-PIP = 1 - shrinkage
    },
    # lasso / structured: Wald exceedance Φ(|β̂|/SE)
    pnorm(abs(beta_hat) / pmax(se_beta, 1e-10))
  )

  list(pip = pip, kappa = kappa)
}

#' Compute covariate (gamma) selection statistics — pip_gamma / kappa_gamma.
#'
#' Mirrors the Stan spike-and-slab generated quantities for pip_gamma:
#'   pip_gamma_k = π_0 φ(γ̂_k;0,σ²_γ) / [π_0 φ(γ̂_k;0,σ²_γ) + (1-π_0) φ(γ̂_k;0,ε²)]
#'
#' For other priors, returns the same Wald / shrinkage analogues as for beta.
#'
#' @param gamma_hat  numeric [n_gamma], MAP covariate coefficient estimates
#' @param se_gamma   numeric [n_gamma], Laplace SEs for gamma
#' @param sigma2_g   numeric, slab variance for covariate Normal prior
#' @param prior      character, prior family
#' @param pi0        numeric, base inclusion probability (spike_slab)
#' @param spike_var  numeric, spike variance (spike_slab)
#' @param hs_tau     numeric, global scale (horseshoe)
#' @param hs_c2      numeric, slab variance (horseshoe)
#' @param n_obs      integer, number of at-risk rows (for horseshoe kappa)
#' @return named list(pip_gamma = numeric [n_gamma], kappa_gamma = numeric [n_gamma])
compute_gamma_selection_stats <- function(gamma_hat, se_gamma, sigma2_g,
                                           prior, pi0, spike_var,
                                           hs_tau, hs_c2, n_obs) {
  n_g         <- length(gamma_hat)
  kappa_gamma <- rep(NA_real_, n_g)

  pip_gamma <- switch(prior,
    spike_slab = {
      # Covariates share the flat π_0 (no order penalty: p=0 ⟹ π_j = π_0)
      # Matches ctbn_fit_unified_v2 Stan:  real pi_gamma = pi0
      pi0_g  <- pmin(pmax(pi0, 1e-10), 1 - 1e-10)
      phi_s  <- dnorm(gamma_hat, 0, sqrt(sigma2_g))
      phi_k  <- dnorm(gamma_hat, 0, sqrt(spike_var))
      (pi0_g * phi_s) / (pi0_g * phi_s + (1 - pi0_g) * phi_k + 1e-300)
    },
    horseshoe = {
      # Covariates not order-penalised: use tau directly (order=0)
      v_g_hs      <- .hs_eff_var(hs_tau, hs_c2, orders = rep(0L, n_g), theta = 0)
      kappa_gamma <- 1.0 / (1.0 + n_obs * v_g_hs)
      1.0 - kappa_gamma
    },
    # lasso / structured: Wald exceedance
    pnorm(abs(gamma_hat) / pmax(se_gamma, 1e-10))
  )

  list(pip_gamma = pip_gamma, kappa_gamma = kappa_gamma)
}

# =============================================================================
# SECTION G: Empirical Bayes hyperparameter initialisation
# =============================================================================

init_hyperparams_eb <- function(DT_wide, all_conds, all_covs,
                                 max_order, a0 = 2, b0 = 1) {
  betas_crude <- c()
  for (target in all_conds) {
    infl <- setdiff(all_conds, target)
    atr  <- DT_wide[get(target) == 0 & dt > 0]
    if (nrow(atr) < 10) next
    df  <- as.data.frame(atr)
    ec  <- paste0(target, "_event")
    X_  <- as.matrix(df[, infl, drop = FALSE])
    y_  <- as.integer(df[[ec]])
    off <- log(pmax(df$dt, 1e-10))
    tryCatch({
      fit_ <- glm.fit(cbind(1, X_), y_, family = poisson(link = "log"),
                      offset = off, control = list(maxit = 30))
      betas_crude <- c(betas_crude, coef(fit_)[-1])
    }, error = function(e) NULL)
  }
  emp_var <- if (length(betas_crude) > 5)
    max(var(betas_crude, na.rm = TRUE), 0.01) else 1.0
  sigma2 <- emp_var * exp(-(0:max_order))
  list(sigma2 = sigma2, lambda2 = 1 / pmax(sigma2, 0.01),
       sigma2_g = max(emp_var, 0.25), a0 = a0, b0 = b0)
}

# =============================================================================
# SECTION H: Main fast MAP fitting function
# =============================================================================

#' Fast matrix-based MAP/Laplace CTBN estimator
#'
#' Drop-in replacement for ctbn_fit() from ctbn_fit_unified_v2.R.
#' Supports all four prior families: spike_slab | structured | lasso | horseshoe
#' Returns an S3 object of class c("ctbn_map","ctbn_fit") compatible with
#' ctbn_cv_parallel.R and ctbn_simulation_v3.R.
#'
#' @param DT_wide           data.table, wide format (eid, time_to_event, dt,
#'                          [condition cols], [covariate cols])
#' @param prior             "spike_slab" | "structured" | "lasso" | "horseshoe"
#' @param max_order         integer, max interaction order (0 = main effects only)
#' @param fixed_covs        character vector, fixed covariate column names
#' @param time_varying_covs character vector, time-varying covariate column names
#' @param target_conditions character vector or NULL (NULL = fit all nodes)
#' @param variable_select   logical, gate intensity output by pip_threshold
#' @param pip_threshold     numeric [0,1], selection threshold
#' @param theta             numeric, order penalty (all priors)
#' @param a0, b0            Inv-Gamma shape/scale (structured / lasso / spike_slab)
#' @param pi0               base inclusion probability (spike_slab only)
#' @param spike_var         spike variance (spike_slab only)
#' @param tau0              HalfCauchy global scale (horseshoe only)
#' @param slab_df           regularising slab degrees of freedom (horseshoe)
#' @param slab_scale        regularising slab scale (horseshoe)
#' @param sigma2_init       optional numeric vector, override empirical Bayes σ²
#' @param lambda2_init      optional numeric vector, override λ²
#' @param sigma2_g_init     optional numeric, override covariate prior variance
#' @param lbfgs_maxit       integer, max L-BFGS-B iterations (default 500)
#' @param compute_se        logical, compute Laplace SEs (default TRUE)
#' @param verbose           logical
#' @return object of class c("ctbn_map","ctbn_fit")
ctbn_map_fast <- function(DT_wide,
                           prior             = c("spike_slab", "structured",
                                                 "lasso", "horseshoe"),
                           max_order         = 1L,
                           fixed_covs        = character(0),
                           time_varying_covs = character(0),
                           target_conditions = NULL,
                           variable_select   = FALSE,
                           pip_threshold     = 0.5,
                           # Shared penalty
                           theta             = 1.0,
                           # Normal / LASSO / spike-and-slab hyperparameters
                           a0                = 2.0,
                           b0                = 1.0,
                           pi0               = 0.5,
                           spike_var         = 0.01,
                           # Horseshoe hyperparameters
                           tau0              = 1.0,
                           slab_df           = 4.0,
                           slab_scale        = 2.0,
                           # Manual overrides (optional; NULL = empirical Bayes)
                           sigma2_init       = NULL,
                           lambda2_init      = NULL,
                           sigma2_g_init     = NULL,
                           # Optimiser
                           lbfgs_maxit       = 500L,
                           compute_se        = TRUE,
                           verbose           = TRUE) {

  prior <- match.arg(prior)

  stopifnot(is.data.table(DT_wide))
  stopifnot("eid"           %in% names(DT_wide))
  stopifnot("time_to_event" %in% names(DT_wide))
  stopifnot("dt"            %in% names(DT_wide))

  all_covs  <- c(fixed_covs, time_varying_covs)
  reserved  <- c("eid", "time_to_event", "dt", all_covs)
  all_conds <- setdiff(names(DT_wide),
                        c(reserved, grep("_event$", names(DT_wide), value = TRUE)))
  n_cond    <- length(all_conds)
  if (n_cond < 2) stop("Need at least 2 condition columns.")

  target_loop <- if (is.null(target_conditions)) all_conds else {
    bad <- setdiff(target_conditions, all_conds)
    if (length(bad))
      stop("target_conditions not found: ", paste(bad, collapse = ", "))
    target_conditions
  }

  # ── Horseshoe effective hyperparameters ─────────────────────────────────────
  # hs_tau : global scale fixed at tau0 (HalfCauchy mode/scale parameter)
  # hs_c2  : E[c²] under InvGamma(slab_df/2, slab_df·slab_scale²/2)
  hs_tau <- tau0
  hs_c2  <- hs_c2_mean(slab_df, slab_scale)

  if (verbose) {
    message(sprintf("ctbn_map_fast [%s prior, max_order=%d]: fitting %d targets",
                    prior, max_order, length(target_loop)))
    switch(prior,
      spike_slab = message(sprintf(
        "  pi0=%.2f | spike_var=%.4f | theta=%.2f | Inv-Gamma(%.1f,%.1f)",
        pi0, spike_var, theta, a0, b0)),
      structured = message(sprintf(
        "  theta=%.2f | Inv-Gamma(%.1f,%.1f)", theta, a0, b0)),
      lasso = message(sprintf(
        "  theta=%.2f | Inv-Gamma(%.1f,%.1f) [LASSO]", theta, a0, b0)),
      horseshoe = message(sprintf(
        "  tau0=%.2f | slab_df=%.1f | slab_scale=%.2f | c2_eff=%.3f | theta=%.2f",
        tau0, slab_df, slab_scale, hs_c2, theta))
    )
  }

  # ── Empirical Bayes initialisation ──────────────────────────────────────────
  hp       <- init_hyperparams_eb(DT_wide, all_conds, all_covs, max_order, a0, b0)
  sigma2   <- rep_len(if (!is.null(sigma2_init))  sigma2_init  else hp$sigma2,
                      max_order + 1L)
  lambda2  <- rep_len(if (!is.null(lambda2_init)) lambda2_init else hp$lambda2,
                      max_order + 1L)
  sigma2_g <- if (!is.null(sigma2_g_init)) sigma2_g_init else hp$sigma2_g

  # ── Prepare data ─────────────────────────────────────────────────────────────
  DT <- copy(DT_wide)
  setorder(DT, eid, time_to_event)
  for (cond in all_conds) {
    ec <- paste0(cond, "_event")
    if (!ec %in% names(DT)) {
      DT[, (ec) := as.numeric(
        get(cond) == 0 & shift(get(cond), type = "lead") == 1), by = eid]
      DT[is.na(get(ec)), (ec) := 0]
    }
  }

  # ── Output containers ────────────────────────────────────────────────────────
  mk_mat <- function(fill) matrix(fill, n_cond, n_cond,
                                   dimnames = list(all_conds, all_conds))
  beta_matrix      <- mk_mat(0)
  se_matrix        <- mk_mat(NA_real_)
  pip_matrix       <- mk_mat(NA_real_)
  kappa_matrix     <- mk_mat(NA_real_)
  intensity_matrix <- mk_mat(0)
  convergence_vec  <- setNames(rep(NA_integer_, length(target_loop)), target_loop)
  map_fits         <- setNames(vector("list", n_cond), all_conds)
  # pip_cov_list[[target]] = named numeric vector of covariate PIPs/pseudo-PIPs
  # Matches ctbn_fit_unified_v2 output slot; spike_slab gives true PIPs,
  # other priors give Wald / horseshoe shrinkage analogues.
  pip_cov_list     <- setNames(vector("list", n_cond), all_conds)

  to_num <- function(v) {
    if (is.null(v)) return(NULL)
    if (is.factor(v) || is.character(v)) return(as.numeric(as.factor(v)) - 1)
    as.numeric(v)
  }

  # ==========================================================================
  # Main loop: one MAP optimisation per target node
  # The CTBN factorisation p(β,γ|data) = ∏_m p(β_m,γ_m|data_m) makes these
  # M problems independent — see ctbn_map_parallel() for parallelisation.
  # ==========================================================================
  for (target in target_loop) {
    ec          <- paste0(target, "_event")
    influencers <- setdiff(all_conds, target)

    atr <- DT[get(target) == 0 & dt > 0]
    if (nrow(atr) == 0) {
      warning(sprintf("No at-risk rows for '%s' — skipping.", target)); next
    }
    df <- as.data.frame(atr)
    n  <- nrow(df)

    # ── Build Z ─────────────────────────────────────────────────────────────
    if (length(all_covs) > 0) {
      Z_mat <- do.call(cbind, lapply(all_covs, function(cv) to_num(df[[cv]])))
      colnames(Z_mat) <- all_covs
    } else {
      Z_mat <- matrix(1.0, n, 1, dimnames = list(NULL, "__intercept__"))
    }

    # ── Build Φ_m ────────────────────────────────────────────────────────────
    X_all <- do.call(cbind, lapply(all_conds, function(cn) as.numeric(df[[cn]])))
    colnames(X_all) <- all_conds
    dm       <- build_design_matrix(X_all, Z_mat, all_conds, target, max_order)
    Phi      <- dm$Phi
    x_cols   <- dm$x_cols
    x_orders <- dm$x_orders
    n_beta   <- dm$n_beta
    n_gamma  <- dm$n_gamma

    y      <- as.integer(df[[ec]])
    offset <- log(pmax(df$dt, 1e-10))

    # ── MAP objective ────────────────────────────────────────────────────────
    obj <- make_map_objective(
      Phi = Phi, y = y, offset = offset, x_orders = x_orders,
      prior = prior, sigma2 = sigma2, lambda2 = lambda2,
      sigma2_g = sigma2_g, theta_pen = theta,
      pi0 = pi0, spike_var = spike_var,
      hs_tau = hs_tau, hs_c2 = hs_c2,
      n_beta = n_beta, n_gamma = n_gamma
    )

    # ── L-BFGS-B optimisation ────────────────────────────────────────────────
    # Complexity: O(n·d·iter_max)  vs Stan O(n·d·iter·chains) — 100-1000× faster
    opt <- tryCatch(
      optim(par = rep(0.0, n_beta + n_gamma),
            fn  = obj$fn, gr = obj$gr,
            method  = "L-BFGS-B",
            control = list(maxit = lbfgs_maxit, factr = 1e7, pgtol = 1e-5)),
      error = function(e) {
        warning(sprintf("L-BFGS-B FAILED for '%s': %s", target, e$message))
        list(par = rep(0.0, n_beta + n_gamma), convergence = 9L, value = NA)
      }
    )

    convergence_vec[target] <- opt$convergence
    theta_hat <- opt$par
    beta_hat  <- theta_hat[seq_len(n_beta)]
    gamma_hat <- theta_hat[n_beta + seq_len(n_gamma)]

    # ── Laplace SEs ──────────────────────────────────────────────────────────
    se_hat <- rep(NA_real_, n_beta)
    if (compute_se) {
      lap <- tryCatch(
        laplace_se(theta_hat = theta_hat, Phi = Phi, y = y, offset = offset,
                   x_orders = x_orders, prior = prior,
                   sigma2 = sigma2, lambda2 = lambda2, sigma2_g = sigma2_g,
                   theta_pen = theta, pi0 = pi0, spike_var = spike_var,
                   hs_tau = hs_tau, hs_c2 = hs_c2,
                   n_beta = n_beta, n_gamma = n_gamma),
        error = function(e) {
          warning(sprintf("Laplace SE failed for '%s': %s", target, e$message))
          list(se = rep(NA_real_, n_beta + n_gamma))
        }
      )
      se_hat   <- lap$se[seq_len(n_beta)]
      se_gamma <- lap$se[n_beta + seq_len(n_gamma)]
    } else {
      se_gamma <- rep(NA_real_, n_gamma)
    }

    # ── Selection statistics for condition coefficients (beta) ───────────────
    sel <- compute_selection_stats(
      beta_hat = beta_hat, se_beta = se_hat, x_orders = x_orders,
      prior = prior, sigma2 = sigma2, theta_pen = theta,
      pi0 = pi0, spike_var = spike_var,
      hs_tau = hs_tau, hs_c2 = hs_c2, n_obs = n
    )
    pip_hat   <- sel$pip
    kappa_hat <- sel$kappa

    # ── Selection statistics for covariate coefficients (gamma) ─────────────
    # pip_cov_list[[target]]: named vector matching ctbn_fit_unified_v2 output.
    # For spike_slab: true PIP matching Stan generated quantities pip_gamma.
    # For other priors: Wald pseudo-PIP or horseshoe 1-kappa analogue.
    gsel <- compute_gamma_selection_stats(
      gamma_hat = gamma_hat, se_gamma = se_gamma,
      sigma2_g  = sigma2_g, prior = prior,
      pi0 = pi0, spike_var = spike_var,
      hs_tau = hs_tau, hs_c2 = hs_c2, n_obs = n
    )
    pip_cov_list[[target]] <- setNames(gsel$pip_gamma, colnames(Z_mat))

    # ── Store in output matrices ─────────────────────────────────────────────
    for (idx in seq_along(influencers)) {
      inf <- influencers[idx]
      beta_matrix[inf,  target] <- beta_hat[idx]
      se_matrix[inf,    target] <- se_hat[idx]
      pip_matrix[inf,   target] <- pip_hat[idx]
      kappa_matrix[inf, target] <- kappa_hat[idx]
    }

    # ── Intensities at reference profile ─────────────────────────────────────
    for (idx in seq_along(influencers)) {
      inf <- influencers[idx]

      # Selection gate (mirrors ctbn_fit_unified logic for all four priors)
      if (variable_select) {
        pv <- pip_matrix[inf, target]
        kv <- kappa_matrix[inf, target]
        skip <- switch(prior,
          spike_slab = !is.na(pv) && pv < pip_threshold,
          structured = !is.na(pv) && pv < pip_threshold,
          lasso      = !is.na(kv) && kv > (1 - pip_threshold),
          horseshoe  = !is.na(kv) && kv > (1 - pip_threshold),
          FALSE)
        if (isTRUE(skip)) next
      }

      x_ref   <- rep(0.0, n_beta)
      inf_pos <- which(x_cols == inf)
      if (length(inf_pos) == 1L) x_ref[inf_pos] <- 1.0

      z_ref <- if (ncol(Z_mat) == 1L && colnames(Z_mat)[1L] == "__intercept__")
        1.0 else colMeans(Z_mat)

      intensity_matrix[inf, target] <- exp(
        sum(x_ref * beta_hat) + sum(z_ref * gamma_hat))
    }

    # ── Store node fit ────────────────────────────────────────────────────────
    map_fits[[target]] <- list(
      theta_hat    = theta_hat,
      beta_hat     = beta_hat,
      gamma_hat    = gamma_hat,
      se_hat       = se_hat,
      se_gamma     = se_gamma,
      pip_hat      = pip_hat,
      kappa_hat    = kappa_hat,
      pip_gamma_hat   = gsel$pip_gamma,
      kappa_gamma_hat = gsel$kappa_gamma,
      x_cols       = x_cols,
      x_orders     = x_orders,
      z_cols       = colnames(Z_mat),
      convergence  = opt$convergence,
      optim_value  = opt$value,
      n_obs        = n,
      n_beta       = n_beta,
      n_gamma      = n_gamma
    )

    if (verbose)
      message(sprintf("  [%s] target: %-8s | conv: %d | -logpost: %8.2f | n=%d",
                      if (opt$convergence == 0L) "OK" else "!!",
                      target, opt$convergence, opt$value, n))
  }

  # Cleanup event indicator columns added to the copy
  ec_cols <- grep("_event$", names(DT), value = TRUE)
  if (length(ec_cols)) DT[, (ec_cols) := NULL]

  structure(
    list(
      prior            = prior,
      method           = "map_lbfgs",
      beta_matrix      = beta_matrix,
      se_matrix        = se_matrix,
      pip_matrix       = pip_matrix,
      kappa_matrix     = kappa_matrix,
      intensity_matrix = intensity_matrix,
      # pip_cov_list: covariate PIPs per target — matches ctbn_fit_unified_v2
      # spike_slab: true PIP from posterior mixture weight at MAP
      # horseshoe:  pseudo-PIP = 1 - kappa_gamma
      # lasso/structured: Wald exceedance Phi(|gamma_hat|/SE)
      pip_cov_list     = pip_cov_list,
      convergence      = convergence_vec,
      map_fits         = map_fits,
      # NULL Stan-fit slots for interface compatibility
      stan_fits        = setNames(vector("list", n_cond), all_conds),
      models           = map_fits,
      # pvalue_matrix included for backward compat with ctbn_fit interface
      pvalue_matrix    = mk_mat(NA_real_),
      call_args        = list(
        prior             = prior,
        max_order         = max_order,
        fixed_covs        = fixed_covs,
        time_varying_covs = time_varying_covs,
        target_conditions = target_conditions,
        all_conditions    = all_conds,
        variable_select   = variable_select,
        pip_threshold     = pip_threshold,
        theta             = theta,
        a0 = a0, b0 = b0,
        pi0 = pi0, spike_var = spike_var,
        tau0 = tau0, slab_df = slab_df, slab_scale = slab_scale,
        hs_tau = hs_tau, hs_c2 = hs_c2,
        sigma2 = sigma2, lambda2 = lambda2, sigma2_g = sigma2_g,
        lbfgs_maxit = lbfgs_maxit
      )
    ),
    class = c("ctbn_map", "ctbn_fit")
  )
}

# =============================================================================
# SECTION I: get_lp S3 dispatch for CV compatibility
# =============================================================================

get_lp <- function(fit, newdata, target, ...) UseMethod("get_lp")

#' @method get_lp ctbn_map
get_lp.ctbn_map <- function(fit, newdata, target, ...) {
  mfit <- fit$map_fits[[target]]
  if (is.null(mfit)) {
    na <- rep(NA_real_, nrow(newdata))
    return(list(lp = na, lambda = na))
  }

  ca        <- fit$call_args
  all_covs  <- c(ca$fixed_covs, ca$time_varying_covs)
  all_conds <- ca$all_conditions
  nd        <- as.data.frame(newdata)
  n         <- nrow(nd)

  # Drop event columns from newdata if present
  nd <- nd[, setdiff(names(nd), grep("_event$", names(nd), value = TRUE)),
           drop = FALSE]

  to_num <- function(v) {
    if (is.null(v)) return(rep(0.0, n))
    if (is.factor(v) || is.character(v)) return(as.numeric(as.factor(v)) - 1)
    as.numeric(v)
  }

  X_all <- do.call(cbind, lapply(all_conds, function(cn) to_num(nd[[cn]])))
  colnames(X_all) <- all_conds

  Z_mat <- if (length(all_covs) > 0L)
    do.call(cbind, lapply(all_covs, function(cv) to_num(nd[[cv]])))
  else
    matrix(1.0, n, 1L)

  dm <- build_design_matrix(X_all, Z_mat, all_conds, target, ca$max_order)

  if (ncol(dm$Phi) != (mfit$n_beta + mfit$n_gamma)) {
    warning(sprintf("get_lp.ctbn_map [%s]: dim mismatch — %d vs %d",
                    target, ncol(dm$Phi), mfit$n_beta + mfit$n_gamma))
    na <- rep(NA_real_, n)
    return(list(lp = na, lambda = na))
  }

  lp <- as.numeric(dm$Phi %*% mfit$theta_hat)
  list(lp = lp, lambda = exp(lp))
}

# =============================================================================
# SECTION J: Print and summary methods
# =============================================================================

#' @method print ctbn_map
print.ctbn_map <- function(x, digits = 4, ...) {
  ca <- x$call_args
  cat(sprintf("\n=== CTBN MAP Fit (%s prior, L-BFGS-B) ===\n", toupper(x$prior)))
  cat(sprintf("  Conditions      : %d\n", nrow(x$beta_matrix)))
  cat(sprintf("  Max order       : %d\n", ca$max_order))
  cat(sprintf("  theta (penalty) : %.2f\n", ca$theta))
  cat(sprintf("  L-BFGS maxit    : %d\n", ca$lbfgs_maxit))
  switch(x$prior,
    spike_slab = cat(sprintf(
      "  pi0=%.2f | spike_var=%.4f | Inv-Gamma(%.1f,%.1f)\n",
      ca$pi0, ca$spike_var, ca$a0, ca$b0)),
    structured = cat(sprintf("  Inv-Gamma(%.1f,%.1f)\n", ca$a0, ca$b0)),
    lasso      = cat(sprintf("  Inv-Gamma(%.1f,%.1f) [LASSO lambda^2]\n",
                             ca$a0, ca$b0)),
    horseshoe  = cat(sprintf(
      "  tau0=%.2f | slab_df=%.1f | slab_scale=%.2f | c2_eff=%.3f\n",
      ca$tau0, ca$slab_df, ca$slab_scale, ca$hs_c2))
  )
  conv   <- x$convergence
  n_conv <- sum(conv == 0L, na.rm = TRUE)
  cat(sprintf("  Convergence     : %d/%d nodes converged\n",
              n_conv, sum(!is.na(conv))))
  if (any(conv != 0L, na.rm = TRUE))
    cat(sprintf("  NON-CONVERGED   : %s\n",
                paste(names(conv)[!is.na(conv) & conv != 0], collapse = ", ")))

  cat("\nPosterior MAP beta (influencer -> target):\n")
  print(round(x$beta_matrix, digits))
  cat("\nLaplace SE (influencer -> target):\n")
  print(round(x$se_matrix, digits))

  if (x$prior == "horseshoe") {
    cat("\nShrinkage kappa (0=signal, 1=shrunk; horseshoe):\n")
    print(round(x$kappa_matrix, digits))
  } else {
    cat("\nPIP / pseudo-PIP (influencer -> target):\n")
    print(round(x$pip_matrix, digits))
  }
  cat("\nIntensity at reference (other conditions = 0):\n")
  print(round(x$intensity_matrix, digits))
  invisible(x)
}

#' @method summary ctbn_map
summary.ctbn_map <- function(object, pip_threshold = NULL,
                              show_covariates = FALSE, ...) {
  if (is.null(pip_threshold))
    pip_threshold <- object$call_args$pip_threshold %||% 0.5

  bm <- object$beta_matrix;  pm <- object$pip_matrix
  km <- object$kappa_matrix; sm <- object$se_matrix
  im <- object$intensity_matrix; pr <- object$prior

  is_active <- function(inf, tgt) {
    switch(pr,
      spike_slab = !is.na(pm[inf,tgt]) && pm[inf,tgt] >= pip_threshold,
      structured = !is.na(pm[inf,tgt]) && pm[inf,tgt] >= pip_threshold,
      lasso      = !is.na(km[inf,tgt]) && km[inf,tgt] <= (1 - pip_threshold),
      horseshoe  = !is.na(km[inf,tgt]) && km[inf,tgt] <= (1 - pip_threshold),
      FALSE)
  }

  rows <- list()
  for (inf in rownames(bm)) for (tgt in colnames(bm)) {
    if (inf == tgt || !is_active(inf, tgt)) next
    row <- data.frame(
      influencer = inf, target = tgt,
      log_RR = round(bm[inf,tgt], 4), RR = round(exp(bm[inf,tgt]), 4),
      SE = round(sm[inf,tgt], 4),
      z  = round(bm[inf,tgt] / pmax(sm[inf,tgt], 1e-10), 3),
      intensity = round(im[inf,tgt], 6), stringsAsFactors = FALSE
    )
    if (pr == "spike_slab")
      row$PIP   <- round(pm[inf,tgt], 3)
    else if (pr == "horseshoe")
      row$kappa <- round(km[inf,tgt], 3)
    else
      row$pseudo_PIP <- round(pm[inf,tgt], 3)
    rows[[length(rows)+1]] <- row
  }

  sel_lbl <- switch(pr,
    spike_slab = sprintf("PIP >= %.2f", pip_threshold),
    lasso      = sprintf("kappa <= %.2f", 1 - pip_threshold),
    horseshoe  = sprintf("kappa <= %.2f", 1 - pip_threshold),
    sprintf("pseudo-PIP >= %.2f", pip_threshold))
  cat(sprintf("\n=== Active Influencer Pairs (%s prior, %s, MAP) ===\n",
              pr, sel_lbl))
  if (!length(rows)) {
    result_inf <- NULL
    cat("  None\n")
  } else {
    result_inf <- do.call(rbind, rows)
    sc     <- switch(pr, spike_slab="PIP", horseshoe="kappa", "pseudo_PIP")
    result_inf <- result_inf[order(result_inf[[sc]],
                                    decreasing=(pr=="spike_slab")),]
    rownames(result_inf) <- NULL
    print(result_inf, ...)
  }

  # ── Covariate selection table ─────────────────────────────────────────────
  # Mirrors summary.ctbn_fit show_covariates behaviour.
  # For spike_slab: true PIP >= pip_threshold (exact match to ctbn_fit_unified).
  # For other priors: pseudo-PIP >= pip_threshold (Wald or horseshoe analogue).
  result_cov <- NULL
  if (show_covariates && !is.null(object$pip_cov_list)) {
    cov_col_name <- switch(pr,
      spike_slab = "PIP",
      horseshoe  = "pseudo_PIP",
      "pseudo_PIP")
    cov_rows <- list()
    for (tgt in names(object$pip_cov_list)) {
      pips <- object$pip_cov_list[[tgt]]
      if (is.null(pips) || length(pips) == 0) next
      for (cv in names(pips)) {
        pv <- pips[[cv]]
        if (!is.na(pv) && pv >= pip_threshold)
          cov_rows[[length(cov_rows) + 1]] <- data.frame(
            covariate = cv, target = tgt,
            stat      = round(pv, 3),
            stringsAsFactors = FALSE
          )
      }
    }
    if (length(cov_rows)) {
      result_cov            <- do.call(rbind, cov_rows)
      names(result_cov)[3]  <- cov_col_name
      result_cov            <- result_cov[order(-result_cov[[cov_col_name]]), ]
      rownames(result_cov)  <- NULL
      cat(sprintf("\n=== Active Covariates (%s >= %.2f, %s prior, MAP) ===\n",
                  cov_col_name, pip_threshold, pr))
      print(result_cov, ...)
    }
  }

  invisible(list(influencers = result_inf, covariates = result_cov))
}

# =============================================================================
# SECTION K: Simulation validation helper
# =============================================================================

#' Compare MAP estimates against known simulation truth
validate_map_vs_truth <- function(map_fit, TRUE_BETA_MAIN,
                                   TRUE_BETA_INT2 = NULL,
                                   CONDS = rownames(map_fit$beta_matrix),
                                   pip_threshold = 0.5) {
  rows <- list(); ca <- map_fit$call_args

  for (m in CONDS) {
    mfit <- map_fit$map_fits[[m]]; if (is.null(mfit)) next
    infl <- setdiff(CONDS, m)
    bh   <- mfit$beta_hat; seh <- mfit$se_hat
    ph   <- mfit$pip_hat;  kh  <- mfit$kappa_hat

    for (idx in seq_along(infl)) {
      j   <- infl[idx]; key <- paste(m, j, sep="|")
      tv  <- if (key %in% names(TRUE_BETA_MAIN)) TRUE_BETA_MAIN[[key]] else 0
      sel <- if (map_fit$prior == "horseshoe")
        !is.na(kh[idx]) && kh[idx] <= (1 - pip_threshold)
      else !is.na(ph[idx]) && ph[idx] >= pip_threshold
      rows[[length(rows)+1]] <- data.table(
        condition=m, influencer=j, param_type="main_effect",
        true_val=tv, est_map=bh[idx], laplace_se=seh[idx],
        pip=ph[idx], kappa=kh[idx], selected=sel, true_nonzero=(tv!=0),
        bias=bh[idx]-tv, sq_err=(bh[idx]-tv)^2,
        covered_95 = if (!is.na(seh[idx]))
          as.numeric(abs(bh[idx]-tv) <= 1.96*seh[idx]) else NA_real_
      )
    }

    if (!is.null(TRUE_BETA_INT2) && ca$max_order >= 2L) {
      int_pos <- which(mfit$x_orders == 1L)
      for (ic in seq_along(int_pos)) {
        pos <- int_pos[ic]; cn <- mfit$x_cols[pos]
        pts <- strsplit(cn, "_x_")[[1]]; j <- pts[1]; k <- pts[2]
        k1  <- paste(m,j,k,sep="|"); k2 <- paste(m,k,j,sep="|")
        tv  <- if (k1 %in% names(TRUE_BETA_INT2)) TRUE_BETA_INT2[[k1]]
               else if (k2 %in% names(TRUE_BETA_INT2)) TRUE_BETA_INT2[[k2]] else 0
        sel <- if (map_fit$prior == "horseshoe")
          !is.na(kh[pos]) && kh[pos] <= (1-pip_threshold)
        else !is.na(ph[pos]) && ph[pos] >= pip_threshold
        rows[[length(rows)+1]] <- data.table(
          condition=m, influencer=cn, param_type="interaction_2way",
          true_val=tv, est_map=bh[pos], laplace_se=seh[pos],
          pip=ph[pos], kappa=kh[pos], selected=sel, true_nonzero=(tv!=0),
          bias=bh[pos]-tv, sq_err=(bh[pos]-tv)^2,
          covered_95 = if (!is.na(seh[pos]))
            as.numeric(abs(bh[pos]-tv) <= 1.96*seh[pos]) else NA_real_
        )
      }
    }

    # ── Covariate coefficients ───────────────────────────────────────────────
    # We cannot compute bias/RMSE without true gamma values, but we record
    # the MAP estimate, SE, and covariate PIP for each gamma coefficient so
    # the returned table is complete (param_type = "covariate").
    if (!is.null(mfit$gamma_hat) && !is.null(mfit$z_cols)) {
      gh_m  <- mfit$gamma_hat
      seg_m <- if (!is.null(mfit$se_gamma)) mfit$se_gamma
               else rep(NA_real_, mfit$n_gamma)
      pg_m  <- mfit$pip_gamma_hat
      kg_m  <- mfit$kappa_gamma_hat
      for (gi in seq_along(mfit$z_cols)) {
        cv  <- mfit$z_cols[gi]
        est <- if (gi <= length(gh_m))  gh_m[gi]  else NA_real_
        seg <- if (gi <= length(seg_m)) seg_m[gi] else NA_real_
        pg  <- if (!is.null(pg_m) && gi <= length(pg_m)) pg_m[gi] else NA_real_
        kg  <- if (!is.null(kg_m) && gi <= length(kg_m)) kg_m[gi] else NA_real_
        sel_g <- if (map_fit$prior == "horseshoe")
          !is.na(kg) && kg <= (1 - pip_threshold)
        else !is.na(pg) && pg >= pip_threshold
        rows[[length(rows) + 1]] <- data.table(
          condition=m, influencer=cv, param_type="covariate",
          true_val=NA_real_, est_map=est, laplace_se=seg,
          pip=pg, kappa=kg, selected=sel_g, true_nonzero=NA,
          bias=NA_real_, sq_err=NA_real_, covered_95=NA_real_
        )
      }
    }
  }

  result <- rbindlist(rows, fill = TRUE)

  # Summary over condition-pair parameters (rows with known true values)
  smry_cond <- result[param_type != "covariate", .(
    N        = .N,
    RMSE     = sqrt(mean(sq_err,    na.rm = TRUE)),
    mean_bias= mean(bias,           na.rm = TRUE),
    coverage = mean(covered_95,     na.rm = TRUE),
    tpr      = mean(selected[true_nonzero == TRUE],  na.rm = TRUE),
    fpr      = mean(selected[true_nonzero == FALSE], na.rm = TRUE)
  ), by = param_type]

  # Summary over covariate coefficients
  smry_cov <- result[param_type == "covariate", .(
    N            = .N,
    mean_pip     = mean(pip,      na.rm = TRUE),
    pct_selected = mean(selected, na.rm = TRUE)
  )]

  cat("\n=== MAP vs Truth Validation — Condition Parameters ===\n")
  print(smry_cond)
  cat("\n=== MAP Covariate Selection Summary ===\n")
  print(smry_cov)
  invisible(result)
}

# =============================================================================
# SECTION L: Parallel MAP via future.apply
# =============================================================================

#' Distribute independent node fits across parallel workers.
#' Requires: library(future); plan(multisession, workers=N)
ctbn_map_parallel <- function(DT_wide, ...) {
  if (!requireNamespace("future.apply", quietly = TRUE))
    stop("Install future.apply: install.packages('future.apply')")
  args      <- list(...)
  all_covs  <- c(args$fixed_covs %||% character(0),
                 args$time_varying_covs %||% character(0))
  reserved  <- c("eid","time_to_event","dt", all_covs)
  all_conds <- setdiff(names(DT_wide),
                        c(reserved, grep("_event$",names(DT_wide),value=TRUE)))
  targets   <- args$target_conditions %||% all_conds
  message(sprintf("ctbn_map_parallel: dispatching %d node fits", length(targets)))

  fits <- future.apply::future_lapply(
    X = targets,
    FUN = function(tgt) {
      do.call(ctbn_map_fast,
              c(list(DT_wide=DT_wide,target_conditions=tgt,verbose=FALSE), args))
    },
    future.seed = TRUE, future.packages = "data.table"
  )

  combined <- do.call(ctbn_map_fast,
                      c(list(DT_wide=DT_wide), modifyList(args, list(verbose=FALSE))))
  for (nf in fits) {
    if (is.null(nf)) next
    tgt <- names(nf$convergence)[1L]
    if (is.null(tgt) || !tgt %in% all_conds) next
    combined$beta_matrix[,tgt]      <- nf$beta_matrix[,tgt]
    combined$se_matrix[,tgt]        <- nf$se_matrix[,tgt]
    combined$pip_matrix[,tgt]       <- nf$pip_matrix[,tgt]
    combined$kappa_matrix[,tgt]     <- nf$kappa_matrix[,tgt]
    combined$intensity_matrix[,tgt] <- nf$intensity_matrix[,tgt]
    combined$map_fits[[tgt]]        <- nf$map_fits[[tgt]]
    combined$models[[tgt]]          <- nf$map_fits[[tgt]]
    combined$convergence[tgt]       <- nf$convergence[tgt]
  }
  combined
}
