% Generated manually for the CTBN-Multimorbidity package.
% Replace with roxygen2-generated file by running devtools::document().
\name{ctbn_map_fast}
\alias{ctbn_map_fast}
\title{Fast MAP + Laplace estimation for Continuous-Time Bayesian Networks}

\description{
Fits a CTBN log-linear Poisson model with one of four order-dependent
structured shrinkage priors (\code{spike_slab}, \code{structured},
\code{lasso}, \code{horseshoe}) by L-BFGS-B optimisation of the posterior
mode, with Laplace-approximated standard errors. Drop-in replacement for
the full Stan / NUTS implementation in \code{\link{ctbn_fit}}, typically
50-200x faster and with identical output slots.
}

\usage{
ctbn_map_fast(DT_wide,
              max_order   = 2L,
              prior       = c("spike_slab", "structured", "lasso", "horseshoe"),
              theta       = 1.0,
              pi0         = 0.10,
              a0          = 2.0,
              b0          = 1.0,
              spike_var   = 1e-4,
              tau0        = NULL,
              slab_df     = 4,
              slab_scale  = 2,
              sigma2_g    = 4.0,
              targets     = NULL,
              covariates  = NULL,
              parallel    = FALSE,
              n_cores     = NULL,
              eb_init     = TRUE,
              control     = list(),
              verbose     = TRUE)
}

\arguments{
  \item{DT_wide}{A \code{\link[data.table]{data.table}} in long-by-target /
    wide-by-parent format. Each row is one at-risk interval; columns
    \code{X_*} are the binary states of candidate parent conditions at the
    start of the interval; \code{target} identifies the target condition;
    \code{N_obs} is the event indicator; \code{T_is} is the exposure width.}
  \item{max_order}{Maximum interaction order \eqn{P}; main effects are
    \eqn{p = 0}, pairwise interactions \eqn{p = 1}, three-way \eqn{p = 2},
    etc. Use \code{max_order = 1L} for main effects only.}
  \item{prior}{Name of the shrinkage prior. See Details.}
  \item{theta}{Order-penalty hyperparameter \eqn{\theta > 0}; the
    effective prior variance contracts by \eqn{e^{-\theta p_j}}.}
  \item{pi0}{Base inclusion probability for the spike-and-slab prior;
    ignored otherwise.}
  \item{a0, b0}{Inv-Gamma hyperparameters for the slab variance / LASSO
    rate parameter.}
  \item{spike_var}{Variance of the spike component (spike-and-slab only).}
  \item{tau0, slab_df, slab_scale}{Regularised-horseshoe hyperparameters
    (see Piironen & Vehtari, 2017). \code{tau0} defaults to
    \eqn{p_0 / Q \cdot 1/\sqrt{n_m}}.}
  \item{sigma2_g}{Variance of the weakly informative N(0, \code{sigma2_g})
    prior on the covariate coefficients \eqn{\gamma_m}.}
  \item{targets}{Optional character vector restricting the fit to a subset
    of target conditions. Defaults to all.}
  \item{covariates}{Optional character vector of covariate column names in
    \code{DT_wide}. Defaults to all columns not matching \code{X_*},
    \code{target}, \code{N_obs}, \code{T_is}, or \code{id}.}
  \item{parallel}{If \code{TRUE}, parallelise across target conditions using
    \code{\link[doParallel]{registerDoParallel}}.}
  \item{n_cores}{Number of worker cores when \code{parallel = TRUE}.
    Defaults to \code{getOption("CTBN.cores")}.}
  \item{eb_init}{Logical; if \code{TRUE} (default), initialise
    hyperparameters by an empirical-Bayes pass before MAP.}
  \item{control}{Optional list passed to \code{\link[stats]{optim}} via
    \code{control}.}
  \item{verbose}{Logical; print per-target progress to the console.}
}

\details{
The CTBN log-linear Poisson likelihood for target condition \eqn{m} is
  \deqn{N_{m,i,s} \sim \mathrm{Poisson}\!\left(\exp(\eta_{m,i,s}) \cdot T_{i,s}\right)}
with linear predictor
  \deqn{\eta_{m,i,s} = \beta^0_m + \sum_{\ell \in \mathrm{Pa}(m)} \beta^{(1)}_{m,\{\ell\}} X_\ell(t^-_s) + \sum_{L \subseteq \mathrm{Pa}(m), 2 \le |L| \le P} \beta^{(|L|)}_{m,L} \prod_{\ell \in L} X_\ell(t^-_s) + \mathbf{Z}^\top \boldsymbol{\gamma}_m.}

The four priors apply an order-dependent contraction \eqn{e^{-\theta p_j}}
to the effective prior variance of \eqn{\beta_j}, encoding the hierarchical
principle that higher-order interactions are penalised more strongly than
main effects:
\itemize{
  \item \strong{spike_slab}: two-component normal mixture with PIPs
        recovered analytically per draw via Bayes' theorem.
  \item \strong{structured}: normal scale mixture with Inv-Gamma slab variance.
  \item \strong{lasso}: Bayesian LASSO (Park & Casella 2008) with
        order-penalised Laplace scale.
  \item \strong{horseshoe}: regularised horseshoe (Piironen & Vehtari 2017)
        with slab cap on the local scales.
}

The posterior mode is computed by L-BFGS-B on the concave negative
log-posterior, and Laplace standard errors by inverting the negative
Hessian. The output structure is identical to that of
\code{\link{ctbn_fit}} so that downstream metric and plotting routines
work unchanged.
}

\value{
A named list with one element per target condition \eqn{m}, each a list
containing:
\describe{
  \item{beta_hat, gamma_hat, beta0_hat}{Posterior MAP estimates.}
  \item{se_beta, se_gamma, se_beta0}{Laplace standard errors.}
  \item{pip_hat}{Exact PIP (spike-and-slab) or pseudo-PIP
                 (continuous priors).}
  \item{kappa_hat}{Posterior shrinkage factor (continuous priors only).}
  \item{x_orders}{Interaction order of each \eqn{\beta_j}.}
  \item{loglik, aic, bic}{Diagnostics evaluated at the MAP.}
  \item{convergence}{The \code{convergence} code returned by \code{optim}.}
}
The returned object has class \code{"ctbn_map"} and supports
\code{\link{print}}, \code{\link{summary}}, and the S3 generic
\code{\link{get_lp}}.
}

\section{References}{
Olaniran, O. R., Paria, S. S., Khondoker, M., MacGregor, A., and Lewin, A.
(2026). Continuous-Time Bayesian Networks with Structured Shrinkage Priors
for Modelling Multimorbidity Trajectories in Large-Scale Electronic Health
Records. \emph{Manuscript}.

Park, T. and Casella, G. (2008). The Bayesian LASSO.
\emph{Journal of the American Statistical Association} \bold{103}, 681-686.

Piironen, J. and Vehtari, A. (2017). Sparsity information and regularization
in the horseshoe and other shrinkage priors.
\emph{Electronic Journal of Statistics} \bold{11}, 5018-5051.
}

\seealso{
\code{\link{ctbn_fit}} for the full Stan / NUTS implementation;
\code{\link{ctbn_cv_compare}} for cross-validated comparison across priors;
\code{\link{compute_F_m}} for survival / cumulative incidence functions.
}

\examples{
\dontrun{
# Small synthetic example shipped with the package
data("toy_DT_wide", package = "CTBNMultimorbidity")

fit <- ctbn_map_fast(toy_DT_wide,
                     max_order = 2L,
                     prior     = "spike_slab",
                     theta     = 1.0,
                     pi0       = 0.10,
                     verbose   = FALSE)

summary(fit, pip_threshold = 0.5)
}
}

\keyword{models}
\keyword{multivariate}
\keyword{survival}
