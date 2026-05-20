# =============================================================================
# tests/testthat.R — testthat entry point
# Run with:  devtools::test()  or  Rscript tests/testthat.R
# =============================================================================

if (requireNamespace("testthat", quietly = TRUE) &&
    requireNamespace("CTBNMultimorbidity", quietly = TRUE)) {

  library(testthat)
  library(CTBNMultimorbidity)
  test_check("CTBNMultimorbidity")

} else {
  # Stand-alone mode: source the R files directly so the tests can be run
  # before the package is installed.
  pkg_root <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), ".."))
  if (!dir.exists(file.path(pkg_root, "R"))) pkg_root <- normalizePath(".")

  source(file.path(pkg_root, "R", "ctbn_map_fast.R"))
  library(testthat)
  test_dir(file.path(pkg_root, "tests", "testthat"))
}
