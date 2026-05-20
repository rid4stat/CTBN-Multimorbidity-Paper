# =============================================================================
# data-raw/make_toy_DT_wide.R
#
# Generates inst/extdata/toy_DT_wide.rds — a small synthetic CTBN panel for
# four conditions (HTN, IHD, DM, RS) and two continuous covariates, designed
# to be runnable by every example in the package documentation.
#
# Re-run with:    Rscript data-raw/make_toy_DT_wide.R
#
# This script must NEVER be wired up to the UK Biobank data extract.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(2026)

# --- Toy network parameters -------------------------------------------------
CONDS  <- c("HTN", "IHD", "DM", "RS")
TRUE_BETA0 <- setNames(log(c(HTN = 0.10, IHD = 0.05, DM = 0.06, RS = 0.07)), CONDS)
# Simple cascade: DM -> HTN -> IHD, RS independent
TRUE_BETA_MAIN <- list(
  HTN = c(DM = 0.7),
  IHD = c(HTN = 0.6, DM = 0.4),
  DM  = c(HTN = 0.5),
  RS  = numeric(0)
)
N_PAT      <- 1500L
T_HORIZON  <- 10
EVAL_TIMES <- c(1, 2, 5, 7, 10)

# --- Patient covariates -----------------------------------------------------
covariates <- data.table(
  id      = seq_len(N_PAT),
  age40   = round(rnorm(N_PAT, 0, 1), 3),   # standardised age at 40
  bmi_std = round(rnorm(N_PAT, 0, 1), 3)
)

# --- Simulate trajectories --------------------------------------------------
simulate_patient <- function(id, cov_row) {
  state <- setNames(rep(0L, length(CONDS)), CONDS)
  t     <- 0
  events <- list()
  step   <- 0.5  # 6-month panel grid

  while (t < T_HORIZON && any(state == 0L)) {
    t_next <- min(t + step, T_HORIZON)
    dt     <- t_next - t

    for (m in CONDS) {
      if (state[[m]] == 1L) next
      # Build log-intensity
      eta <- TRUE_BETA0[[m]] + 0.10 * cov_row$age40 + 0.08 * cov_row$bmi_std
      for (pa in names(TRUE_BETA_MAIN[[m]]))
        eta <- eta + TRUE_BETA_MAIN[[m]][[pa]] * state[[pa]]

      lambda <- exp(eta)
      N_is   <- rbinom(1, 1, 1 - exp(-lambda * dt))

      events[[length(events) + 1L]] <- data.table(
        id      = id,
        target  = m,
        T_is    = dt,
        N_obs   = N_is,
        X_HTN   = state[["HTN"]],
        X_IHD   = state[["IHD"]],
        X_DM    = state[["DM"]],
        X_RS    = state[["RS"]],
        age40   = cov_row$age40,
        bmi_std = cov_row$bmi_std
      )

      if (N_is == 1L) state[[m]] <- 1L
    }

    t <- t_next
  }

  rbindlist(events)
}

cat(sprintf("Simulating %d patients...\n", N_PAT))
toy_DT_wide <- rbindlist(
  lapply(seq_len(N_PAT), function(i) simulate_patient(i, covariates[id == i]))
)

cat("Toy dataset:\n")
cat("  rows:    ", nrow(toy_DT_wide), "\n")
cat("  targets: ", paste(unique(toy_DT_wide$target), collapse = ", "), "\n")
cat("  events:  ", toy_DT_wide[, sum(N_obs)], "\n")

# Persist as both .rds (for the package) and .csv.gz (for inspection)
dir.create("inst/extdata", showWarnings = FALSE, recursive = TRUE)
saveRDS(toy_DT_wide, "inst/extdata/toy_DT_wide.rds", compress = "xz")
data.table::fwrite(toy_DT_wide, "inst/extdata/toy_DT_wide.csv.gz")

cat("Written:\n")
cat("  inst/extdata/toy_DT_wide.rds\n")
cat("  inst/extdata/toy_DT_wide.csv.gz\n")
