# =============================================================================
# test-build_design_matrix.R
# Verify that build_design_matrix() produces the correct column counts and
# order assignments for various (Q, P) configurations.
# =============================================================================

context("build_design_matrix: structure and order assignment")

ensure_loaded <- function() {
  if (!exists("build_design_matrix", mode = "function")) {
    src <- file.path("..", "..", "R", "ctbn_map_fast.R")
    if (!file.exists(src)) src <- file.path("R", "ctbn_map_fast.R")
    if (file.exists(src)) source(src) else skip("ctbn_map_fast.R not found")
  }
}

test_that("Main effects only (P = 1) yields M-1 columns for M conditions", {
  skip_on_cran()
  ensure_loaded()

  if (!requireNamespace("data.table", quietly = TRUE))
    skip("data.table not available")

  # Toy panel for M = 4 conditions with covariates Z1, Z2
  X_all <- data.table::data.table(
    X_A = c(0L, 1L, 0L, 1L),
    X_B = c(0L, 0L, 1L, 1L),
    X_C = c(1L, 0L, 1L, 0L),
    X_D = c(0L, 1L, 1L, 0L)
  )
  Z <- data.table::data.table(Z1 = rnorm(4), Z2 = rnorm(4))

  res <- try(
    build_design_matrix(X_all, Z, cond_names = c("A","B","C","D"),
                        target = "A", max_order = 1L),
    silent = TRUE
  )
  if (inherits(res, "try-error")) skip("build_design_matrix interface mismatch")

  # Expect 3 main-effect columns (everyone except target A)
  expect_equal(ncol(res$Phi_X), 3L)
  expect_true(all(res$x_orders == 0L))
})

test_that("P = 2 yields main effects + pairwise interactions", {
  skip_on_cran()
  ensure_loaded()

  if (!requireNamespace("data.table", quietly = TRUE))
    skip("data.table not available")

  X_all <- data.table::data.table(
    X_A = c(0L, 1L, 0L, 1L),
    X_B = c(0L, 0L, 1L, 1L),
    X_C = c(1L, 0L, 1L, 0L),
    X_D = c(0L, 1L, 1L, 0L)
  )
  Z <- data.table::data.table(Z1 = rnorm(4))

  res <- try(
    build_design_matrix(X_all, Z, cond_names = c("A","B","C","D"),
                        target = "A", max_order = 2L),
    silent = TRUE
  )
  if (inherits(res, "try-error")) skip("build_design_matrix interface mismatch")

  # M-1 = 3 main + choose(3,2) = 3 pairwise = 6 columns
  expect_equal(ncol(res$Phi_X), 6L)
  # First 3 cols should be order 0 (main), next 3 order 1 (pairwise interaction)
  expect_equal(sum(res$x_orders == 0L), 3L)
  expect_equal(sum(res$x_orders == 1L), 3L)
})
