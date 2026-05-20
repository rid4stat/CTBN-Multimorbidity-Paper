# =============================================================================
# zzz.R — Package load helpers and global option defaults
# =============================================================================

.onLoad <- function(libname, pkgname) {
  # Sensible default options for the CTBN routines. Each can be overridden
  # by the user with options(<name> = ...).
  op <- options()
  op_ctbn <- list(
    CTBN.cores              = max(1L, parallel::detectCores() - 1L),
    CTBN.stan_chains        = 4L,
    CTBN.stan_iter          = 2000L,
    CTBN.stan_warmup        = 1000L,
    CTBN.stan_adapt_delta   = 0.95,
    CTBN.stan_max_treedepth = 12L,
    CTBN.pip_threshold      = 0.50,
    CTBN.eval_times         = c(1, 2, 5, 7, 10),
    CTBN.verbose            = TRUE
  )
  toset <- !(names(op_ctbn) %in% names(op))
  if (any(toset)) options(op_ctbn[toset])
  invisible()
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "CTBN-Multimorbidity (v", utils::packageVersion(pkgname), ")\n",
    "  Use ctbn_map_fast() for the recommended fast MAP fit, or\n",
    "  ctbn_fit() for full Stan / NUTS MCMC. See ?ctbn_map_fast for help.\n",
    "  Repository: https://github.com/InflAim/CTBN-Multimorbidity")
}
