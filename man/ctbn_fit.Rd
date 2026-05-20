% Generated manually for the CTBN-Multimorbidity package.
\name{ctbn_fit}
\alias{ctbn_fit}
\title{Full Bayesian CTBN fit via Stan / NUTS}

\description{
Fits a Continuous-Time Bayesian Network log-linear Poisson model with
one of four order-dependent structured shrinkage priors using the
No-U-Turn Sampler (NUTS) in Stan. This is the inferential workhorse of
the manuscript and produces the full joint posterior on
\eqn{(\boldsymbol{\beta}^m, \boldsymbol{\gamma}_m)} for each target
condition, together with PIPs (spike-and-slab) or posterior shrinkage
factors (continuous priors).
}

\usage{
ctbn_fit(DT_wide,
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
         targets     = NULL,
         covariates  = NULL,
         chains      = 4L,
         iter        = 2000L,
         warmup      = 1000L,
         adapt_delta = 0.95,
         max_treedepth = 12L,
         cores       = NULL,
         seed        = 2026L,
         verbose     = TRUE)
}

\arguments{
  \item{DT_wide, max_order, prior, theta, pi0, a0, b0, spike_var, tau0,
        slab_df, slab_scale, targets, covariates}{See
        \code{\link{ctbn_map_fast}}; identical interface.}
  \item{chains, iter, warmup}{Number of Stan chains, total iterations
        (including warm-up), and warm-up iterations per chain.}
  \item{adapt_delta, max_treedepth}{Stan NUTS control parameters.}
  \item{cores}{Number of cores for parallel chains.}
  \item{seed}{Integer random seed.}
  \item{verbose}{Logical; print Stan progress.}
}

\value{
A named list with one \code{stanfit}-like object per target condition,
augmented with the same downstream summary slots used by
\code{\link{ctbn_map_fast}}.
}

\seealso{
\code{\link{ctbn_map_fast}} for the fast MAP alternative;
\code{\link{ctbn_cv_compare}} for cross-validation;
\code{\link[rstan]{sampling}} for the underlying Stan interface.
}

\examples{
\dontrun{
data("toy_DT_wide", package = "CTBNMultimorbidity")
fit_stan <- ctbn_fit(toy_DT_wide,
                     max_order = 2L,
                     prior     = "spike_slab",
                     theta     = 1.0,
                     chains    = 2,
                     iter      = 1000,
                     warmup    = 500,
                     seed      = 2026)
summary(fit_stan)
}
}

\keyword{models}
\keyword{multivariate}
\keyword{survival}
