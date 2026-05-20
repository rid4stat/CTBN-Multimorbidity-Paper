% Generated manually for the CTBN-Multimorbidity package.
\name{compute_F_m}
\alias{compute_F_m}
\title{Posterior-mean cumulative incidence function for a target condition}

\description{
Computes the cumulative incidence function \eqn{F_m(\tau \mid \mathbf{x}_i,
\mathbf{z}_i)} for target condition \eqn{m} at horizon \eqn{\tau}, given
the patient covariate vector and parent-state trajectory. Both the fully
Bayesian posterior-mean estimator (averaged over MCMC draws or MAP +
Laplace samples) and the plug-in estimator (evaluated at posterior mean
coefficients) are supported.
}

\usage{
compute_F_m(fit,
            target,
            patient_dt,
            tau,
            zbar       = NULL,
            type       = c("plugin", "posterior_mean"),
            n_draws    = 1000L)
}

\arguments{
  \item{fit}{A fitted CTBN object from \code{\link{ctbn_fit}} or
        \code{\link{ctbn_map_fast}}.}
  \item{target}{Character name of the target condition.}
  \item{patient_dt}{One-row \code{data.table} with the patient's parent-state
        covariates encoded in \code{X_*} columns.}
  \item{tau}{Numeric scalar or vector of evaluation horizons in years.}
  \item{zbar}{Optional vector of covariate values; if \code{NULL}, the
        population-average covariate vector is used.}
  \item{type}{\code{"plugin"} (default) evaluates at the posterior mean of
        the coefficients; \code{"posterior_mean"} averages \eqn{F_m} over
        \code{n_draws} posterior draws.}
  \item{n_draws}{Number of posterior draws when \code{type =
        "posterior_mean"}.}
}

\details{
Because the intensity is piecewise constant between observations, the
cumulative incidence simplifies to
  \deqn{F_m(\tau \mid \mathbf{x}_i, \mathbf{z}_i) = 1 - \exp\!\left(-\sum_{s:\,t_{i,s} < \tau} q_{m,i,s} \cdot T_{i,s}^*\right),}
where \eqn{T_{i,s}^* = \min(t_{i,s+1}, \tau) - t_{i,s}}.

The paper shows that across the UK Biobank cohort the plug-in and fully
Bayesian estimators agree to within \eqn{5 \times 10^{-4}} on the
cumulative incidence scale and yield identical patient-level rankings;
\code{"plugin"} is therefore used by default.
}

\value{
A numeric vector of length \code{length(tau)} giving the cumulative
incidence at each requested horizon.
}

\seealso{
\code{\link{compute_interaction_effects}} for the §5.4 synergy decomposition.
}

\examples{
\dontrun{
data("toy_DT_wide", package = "CTBNMultimorbidity")
fit <- ctbn_map_fast(toy_DT_wide, max_order = 2L, prior = "spike_slab",
                     theta = 1.0, verbose = FALSE)
F5 <- compute_F_m(fit, target = "IHD",
                  patient_dt = toy_DT_wide[1, ],
                  tau = 5)
F5
}
}

\keyword{models}
\keyword{survival}
