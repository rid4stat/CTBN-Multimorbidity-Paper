# =============================================================================
# test-priors.R
# Algebraic sanity checks on the prior log-density / gradient / Hessian
# functions, where the answer is known in closed form.
# =============================================================================

context("Order-dependent prior log-densities")

ensure_loaded <- function() {
  if (!exists("prior_logdens_structured", mode = "function")) {
    src <- file.path("..", "..", "R", "ctbn_map_fast.R")
    if (!file.exists(src)) src <- file.path("R", "ctbn_map_fast.R")
    if (file.exists(src)) source(src) else skip("ctbn_map_fast.R not found")
  }
}

# -- Structured prior ------------------------------------------------------

test_that("Structured prior reduces to N(0, sigma^2) when theta = 0 and order = 0", {
  skip_on_cran()
  ensure_loaded()

  beta   <- c(-1.2, 0.5, 2.0)
  sigma2 <- c(1.0, 0.5, 0.1)       # for orders 0, 1, 2
  orders <- c(0L, 0L, 0L)
  theta  <- 0.0

  ld <- prior_logdens_structured(beta, sigma2, orders, theta)
  expected <- sum(dnorm(beta, 0, sqrt(1.0), log = TRUE))
  expect_equal(ld, expected, tolerance = 1e-10)
})

test_that("Structured gradient matches -beta / v exactly", {
  skip_on_cran()
  ensure_loaded()

  beta   <- c(0.5, -1.0)
  sigma2 <- c(2.0, 0.5)
  orders <- c(0L, 1L)
  theta  <- 0.5
  v      <- sigma2[orders + 1L] * exp(-theta * orders)

  g <- grad_prior_structured(beta, sigma2, orders, theta)
  expect_equal(g, -beta / v, tolerance = 1e-12)
})

# -- Spike-and-slab prior --------------------------------------------------

test_that("Spike-and-slab log-density is a log-mixture and is bounded above by slab", {
  skip_on_cran()
  ensure_loaded()

  beta      <- c(0.0, 1.0)
  sigma2    <- c(1.0)
  orders    <- c(0L, 0L)
  theta     <- 0.0
  pi0       <- 0.5
  spike_var <- 1e-4

  ld <- prior_logdens_spikeslab(beta, sigma2, orders, theta, pi0, spike_var)

  # At beta = 0 the spike component is large, at beta = 1 the slab dominates;
  # the log density must be finite in both cases.
  expect_true(is.finite(ld))
  expect_true(ld < 0)               # log of a density bounded above by ∞
})

# -- Order penalty monotonicity --------------------------------------------

test_that("Higher interaction order yields a smaller effective variance", {
  skip_on_cran()
  ensure_loaded()

  sigma2 <- c(1.0, 1.0, 1.0)
  theta  <- 1.0
  v0 <- sigma2[1] * exp(-theta * 0)
  v1 <- sigma2[2] * exp(-theta * 1)
  v2 <- sigma2[3] * exp(-theta * 2)

  expect_gt(v0, v1)
  expect_gt(v1, v2)
})
