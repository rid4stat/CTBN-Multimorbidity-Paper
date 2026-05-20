#!/usr/bin/env Rscript
# =============================================================================
# run_all.R — One-shot driver to reproduce the manuscript end-to-end
# =============================================================================
#
# Usage:
#   Rscript run_all.R [stage]
#
# where [stage] is one of:
#   sim_map    : Fast MAP simulation study (recommended quick check)
#   sim_stan   : Full Stan simulation study (paper; hours)
#   ukb_stan   : UK Biobank empirical analysis with Stan (paper)
#   ukb_map    : UK Biobank empirical analysis with MAP (scalability check)
#   all        : Run everything in dependency order (default)
#
# The script assumes:
#   * It is invoked from the repository root.
#   * For the `ukb_*` stages, the derived InflAim cohort `DT_wide` is on
#     disk at `data/DT_wide.rds` (NOT shipped — see README §Data availability).
# =============================================================================

suppressPackageStartupMessages({
  library(crayon)
})

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
log_msg <- function(...) {
  cat(crayon::cyan(sprintf("[run_all] %s\n", paste0(..., collapse = ""))))
}

run_script <- function(path) {
  if (!file.exists(path)) stop("Script not found: ", path, call. = FALSE)
  log_msg("→ Sourcing ", path)
  t0 <- Sys.time()
  source(path, chdir = TRUE, echo = FALSE)
  log_msg("  done in ", format(Sys.time() - t0, digits = 3))
}

ensure_dir <- function(d) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
args  <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args) > 0) args[[1L]] else "all"
valid <- c("sim_map", "sim_stan", "ukb_stan", "ukb_map", "all")
if (!stage %in% valid)
  stop("Unknown stage: ", stage, "\nExpected one of: ", paste(valid, collapse = ", "))

log_msg("Selected stage: ", crayon::bold(stage))
ensure_dir("results")

# ---------------------------------------------------------------------------
# Stage runners
# ---------------------------------------------------------------------------
run_sim_map <- function() {
  log_msg(crayon::bold("Stage: sim_map — MAP simulation study"))
  run_script("scripts/simulation/01_map_simulation.R")
  run_script("scripts/simulation/02_map_simulation_results.R")
}

run_sim_stan <- function() {
  log_msg(crayon::bold("Stage: sim_stan — Stan simulation study (paper)"))
  run_script("scripts/simulation/03_stan_simulation.R")
  run_script("scripts/simulation/04_stan_simulation_results.R")
}

run_ukb_stan <- function() {
  log_msg(crayon::bold("Stage: ukb_stan — UK Biobank Stan empirical analysis"))
  if (!file.exists("data/DT_wide.rds"))
    log_msg(crayon::yellow(
      "  ⚠ data/DT_wide.rds not found. The empirical scripts assume the\n",
      "    33,558-participant InflAim cohort has been prepared from\n",
      "    UK Biobank primary care data per Section 4 of the paper."))
  run_script("scripts/ukb_analysis/04_stan_results.R")
}

run_ukb_map <- function() {
  log_msg(crayon::bold("Stage: ukb_map — UK Biobank MAP empirical analysis"))
  run_script("scripts/ukb_analysis/03_map_results.R")
}

# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------
switch(stage,
  "sim_map"  = run_sim_map(),
  "sim_stan" = run_sim_stan(),
  "ukb_stan" = run_ukb_stan(),
  "ukb_map"  = run_ukb_map(),
  "all"      = {
    run_sim_map()
    run_sim_stan()
    run_ukb_map()
    run_ukb_stan()
  })

log_msg(crayon::green("All requested stages completed successfully."))
