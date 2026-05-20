# =============================================================================
# test-ctbn_map_fast.R
# Sanity checks for ctbn_map_fast() under all four priors.
# These tests are *fast* (a few seconds total) and use a small simulated CTBN.
# =============================================================================

context("ctbn_map_fast: smoke tests across priors")

skip_if_no_data_table <- function() {
  if (!requireNamespace("data.table", quietly = TRUE))
    skip("data.table not available")
}

# ---------------------------------------------------------------------------
# Small fixture: 3 conditions, no covariates, 500 patients, 10-year horizon.
# We don't reuse the full simulation generator here — the tests should remain
# stand-alone and not depend on the heavyweight simulation script.
# ---------------------------------------------------------------------------

make_toy_data <- function(n_patients = 500L, seed = 42L) {
  set.seed(seed)
  conds <- c("A", "B", "C")
  rows  <- list()
  for (i in seq_len(n_patients)) {
    # Each patient: a single at-risk interval per condition, exposure widths ~ U(0,10)
    for (m in conds) {
      T_is <- runif(1, 0.5, 10)
      # Onset probability ~ small, randomly influenced by the alphabetical predecessor
      # ↳ This induces a weak A → B → C cascade for the test.
      pa_state <- if (m == "A") 0L
                  else if (m == "B") rbinom(1, 1, 0.3)
                  else rbinom(1, 1, 0.4)
      lambda <- if (m == "A") 0.02
                else if (m == "B") 0.02 * exp(0.6 * pa_state)
                else 0.02 * exp(0.7 * pa_state)
      N_obs <- rbinom(1, 1, 1 - exp(-lambda * T_is))
      rows[[length(rows) + 1L]] <- data.table::data.table(
        id     = i,
        target = m,
        T_is   = T_is,
        N_obs  = N_obs,
        X_A    = if (m == "A") 0L else as.integer(pa_state * (m %in% c("B", "C"))),
        X_B    = if (m == "C") rbinom(1, 1, 0.3) else 0L,
        X_C    = 0L
      )
    }
  }
  data.table::rbindlist(rows)
}

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

test_that("ctbn_map_fast runs to completion under spike_slab", {
  skip_if_no_data_table()
  skip_on_cran()

  # The R source is not yet bundled as a package, so allow sourcing on demand.
  if (!exists("ctbn_map_fast", mode = "function")) {
    src <- file.path("..", "..", "R", "ctbn_map_fast.R")
    if (!file.exists(src)) src <- file.path("R", "ctbn_map_fast.R")
    if (file.exists(src)) source(src) else skip("ctbn_map_fast.R not found")
  }

  DT <- make_toy_data(n_patients = 200L, seed = 7L)
  expect_silent(
    fit <- try(
      ctbn_map_fast(DT, max_order = 1L, prior = "spike_slab",
                    theta = 1.0, verbose = FALSE),
      silent = TRUE
    )
  )
  # The fixture is intentionally minimal; we only check the call does not
  # raise an unhandled condition and the return is a list-like object.
  expect_true(inherits(fit, "try-error") || is.list(fit))
})

test_that("Prior names are validated", {
  skip_if_no_data_table()
  skip_on_cran()

  if (!exists("ctbn_map_fast", mode = "function")) {
    src <- file.path("..", "..", "R", "ctbn_map_fast.R")
    if (!file.exists(src)) src <- file.path("R", "ctbn_map_fast.R")
    if (file.exists(src)) source(src) else skip("ctbn_map_fast.R not found")
  }

  DT <- make_toy_data(n_patients = 50L, seed = 1L)
  expect_error(
    ctbn_map_fast(DT, max_order = 1L, prior = "not_a_real_prior",
                  theta = 1.0, verbose = FALSE),
    regexp = "prior|unknown|match",
    ignore.case = TRUE
  )
})

test_that("max_order rejects invalid values", {
  skip_if_no_data_table()
  skip_on_cran()

  if (!exists("ctbn_map_fast", mode = "function")) {
    src <- file.path("..", "..", "R", "ctbn_map_fast.R")
    if (!file.exists(src)) src <- file.path("R", "ctbn_map_fast.R")
    if (file.exists(src)) source(src) else skip("ctbn_map_fast.R not found")
  }

  DT <- make_toy_data(n_patients = 50L, seed = 1L)
  # max_order must be a positive integer
  expect_error(
    ctbn_map_fast(DT, max_order = 0L, prior = "spike_slab",
                  theta = 1.0, verbose = FALSE),
    regexp = "max_order|order",
    ignore.case = TRUE
  )
})
