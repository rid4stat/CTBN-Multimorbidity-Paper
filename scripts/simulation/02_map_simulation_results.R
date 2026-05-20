# =============================================================================
# ctbn_map_simulation_results.R
# CTBN MAP Simulation Study — Results Aggregation, Tables and Figures
# =============================================================================
#
# THREE FOCAL SCENARIOS
# ---------------------
# A1 – Underspecified  : P_fit=1 (main effects only), P_true=2 (2-way DGP)
# A2 – Correctly spec. : P_fit=2, P_true=2  [baseline]
# A3 – Overspecified   : P_fit=3, P_true=2  (3-way terms fitted, none in DGP)
#
# FOUR MAP PRIOR FAMILIES
# -----------------------
# spike_slab  |  structured  |  lasso  |  horseshoe
#
# OUTPUTS PRODUCED
# ----------------
# Tables (LaTeX):
#   tab_s0_scenario_key.tex          — scenario definitions
#   tab_s1_recovery_summary.tex      — RMSE/bias/coverage by prior × scenario
#   tab_s2_covariate_recovery.tex    — covariate coefficient recovery
#   tab_s3_selection_metrics.tex     — TPR/FPR/selection-AUC by prior × scenario
#   tab_s4_pred_metrics.tex          — Poisson LL / Brier by prior × scenario
#   tab_s5_master.tex                — master comparison (all metrics)
#
# Figures (PDF + PNG):
#   fig1_rmse_by_prior_scenario      — RMSE by prior × scenario × param_type
#   fig2_bias_coverage               — Bias and 95% CI coverage
#   fig3_pip_covariate               — Covariate PIP by prior × scenario
#   fig4_selection_tpr_fpr           — TPR and FPR for variable selection
#   fig5_pred_metrics                — Poisson LL and Brier score
#   fig6_prior_ranking_heatmap       — Rank heat-map across all metrics
#
# USAGE
# -----
#   source("ctbn_map_fast.R")
#   source("ctbn_map_simulation.R")
#   source("ctbn_map_simulation_results.R")
#
#   # Run everything from scratch:
#   run_and_analyse(n_cores = 20L)
#
#   # Analyse already-completed replicates:
#   build_all_results()
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(xtable)
  library(patchwork)
  library(crayon)
})

# source("ctbn_map_fast.R")       # uncomment if not already sourced
# source("ctbn_map_simulation.R") # uncomment if not already sourced

# =============================================================================
# 0.  Infrastructure
# =============================================================================

RESULTS_LOG <- "sim_results/logs/results_run.log"
MANIFEST    <- "sim_results/output_manifest.txt"

for (d in c("sim_results","sim_results/figures","sim_results/tables",
            "sim_results/logs"))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

results_log <- function(msg, level = "INFO") {
  line <- sprintf("[%s] [%s] %s\n",
                  format(Sys.time(), "%H:%M:%S"), level, msg)
  cat(line)
  tryCatch(cat(line, file = RESULTS_LOG, append = TRUE),
           error = function(e) invisible(NULL))
}

section <- function(name, expr) {
  cat(sprintf("\n%s %s\n%s\n",
              format(Sys.time(), "%H:%M:%S"),
              crayon::bold(name),
              strrep("\u2500", nchar(name) + 10)))
  t0  <- proc.time()["elapsed"]
  out <- force(expr)
  elapsed <- proc.time()["elapsed"] - t0
  results_log(sprintf("%s \u2014 done in %s", name, format_duration(elapsed)))
  cat(sprintf("  \u2713 Done in %s\n", format_duration(elapsed)))
  tryCatch(cat(sprintf("%s | %s\n", format(Sys.time(), "%H:%M:%S"), name),
               file = MANIFEST, append = TRUE),
           error = function(e) invisible(NULL))
  invisible(out)
}

print_progress <- function(i, n, label = "") {
  cat(sprintf("\r  %s  %s", make_progress_bar(i, n, width = 30), label))
  if (i == n) cat("\n")
  utils::flush.console()
}

save_sim_fig <- function(p, name, w = 12, h = 7) {
  for (ext in c("pdf","png")) {
    path <- file.path("sim_results/figures", paste0(name,".",ext))
    if (ext == "pdf")
      ggsave(path, p, width=w, height=h, units="in", device=cairo_pdf)
    else
      ggsave(path, p, width=w, height=h, units="in", dpi=300)
  }
  tryCatch(cat(sprintf("  FIGURE: sim_results/figures/%s.{pdf,png}\n", name),
               file = MANIFEST, append = TRUE), error = function(e) NULL)
  message("  \u2192 figures/", name)
}

save_sim_tex <- function(xt, name, ...) {
  path <- file.path("sim_results/tables", paste0(name,".tex"))
  print(xt, file=path, floating=FALSE, include.rownames=FALSE, ...)
  tryCatch(cat(sprintf("  TABLE:  sim_results/tables/%s.tex\n", name),
               file = MANIFEST, append = TRUE), error = function(e) NULL)
  message("  \u2192 tables/", name)
}

writeLines(sprintf("Results run started: %s\n%s",
                   format(Sys.time()), strrep("=",60)), RESULTS_LOG)
writeLines(sprintf("Output manifest \u2014 %s\n%s",
                   format(Sys.time()), strrep("=",60)), MANIFEST)

# =============================================================================
# 1.  Metadata
# =============================================================================

SCENARIO_META <- data.table(
  scenario    = c("A1",   "A2",   "A3"),
  p_fit       = c(1L,     2L,     3L),
  p_true      = c(2L,     2L,     2L),
  theta       = c(1.0,    1.0,    1.0),
  spec_label  = c("Underspecified\n(P_fit=1, P_true=2)",
                  "Correctly specified\n(P_fit=2, P_true=2)",
                  "Overspecified\n(P_fit=3, P_true=2)"),
  short_label = c("A1: Under (P_fit=1)",
                  "A2: Correct (P_fit=2)",
                  "A3: Over (P_fit=3)")
)

PRIOR_META <- data.table(
  prior       = c("spike_slab","structured","lasso","horseshoe"),
  prior_label = c("Spike-and-Slab","Structured Normal",
                  "Bayesian LASSO","Reg. Horseshoe"),
  prior_short = c("SS","SN","LASSO","HS")
)

SCENARIO_LEVELS <- SCENARIO_META$scenario
PRIOR_LEVELS    <- PRIOR_META$prior
PRIOR_LABELS    <- setNames(PRIOR_META$prior_label, PRIOR_META$prior)

# PRIOR_PALETTE must be keyed by the DISPLAY LABELS that as_prior_factor()
# assigns as factor levels — NOT the raw prior strings.
# scale_colour_manual() matches palette names against the factor's level labels,
# so if the factor levels are "Spike-and-Slab", "Structured Normal", etc.,
# the palette names must match those strings exactly.
PRIOR_PALETTE_RAW <- c(spike_slab="#E41A1C", structured="#377EB8",
                        lasso="#4DAF4A",      horseshoe="#984EA3")
PRIOR_PALETTE     <- setNames(PRIOR_PALETTE_RAW, PRIOR_LABELS[names(PRIOR_PALETTE_RAW)])

PTYPE_ORDER <- c("main_effect","interaction_2way","interaction_3way","covariate")
PTYPE_LABELS <- c(main_effect      = "Main effects",
                  interaction_2way = "2-way interactions",
                  interaction_3way = "3-way interactions",
                  covariate        = "Covariates")

# =============================================================================
# 2.  Orchestration: run all scenarios × priors
# =============================================================================

#' Run all three simulation scenarios across all four MAP prior families.
#'
#' Each (scenario, prior) cell calls run_simulation() from ctbn_map_simulation.R.
#' The critical requirement is that ctbn_map_simulation.R uses
#'   rep_<scenario>_<prior>_<NNN>.rds
#' as the replicate filename so each (scenario, prior) cell writes to a
#' distinct path.  The version delivered alongside this file enforces this.
#' If you have legacy rep_<scenario>_<NNN>.rds files from single-prior runs,
#' call rename_legacy_replicates() first to tag them with the correct prior.
#'
#' @param scenarios   character vector; default all three focal scenarios
#' @param priors      character vector; default all four prior families
#' @param n_rep       integer, replicates per (scenario, prior) cell
#' @param n_sim       integer, subjects per replicate
#' @param k_fold      integer, CV folds
#' @param n_cores     integer, inner doParallel workers per cell
#' @param theta       numeric, order-penalty parameter
#' @param eval_times  numeric vector
#' @param base_seed   integer
run_all_scenarios_all_priors <- function(
    scenarios  = SCENARIO_META$scenario,
    priors     = PRIOR_META$prior,
    n_rep      = N_REP,
    n_sim      = N_SIM,
    k_fold     = K_FOLD,
    n_cores    = 20L,
    theta      = THETA_FIT,
    eval_times = EVAL_TIMES,
    base_seed  = SIM_SEED) {

  grid <- as.data.table(expand.grid(scenario = scenarios, prior = priors,
                                     stringsAsFactors = FALSE))
  grid <- merge(grid, SCENARIO_META[, .(scenario, p_fit, p_true)],
                by = "scenario")

  n_cells <- nrow(grid)
  cat(sprintf("\n%s\nrun_all_scenarios_all_priors:\n", strrep("=", 70)))
  cat(sprintf("  %d scenarios \u00d7 %d priors = %d cells\n",
              length(scenarios), length(priors), n_cells))
  cat(sprintf("  %d replicates \u00d7 %d subjects per cell\n",
              n_rep, n_sim))
  cat(sprintf("  Files: rep_<scenario>_<prior>_<NNN>.rds\n%s\n\n",
              strrep("=", 70)))

  for (i in seq_len(n_cells)) {
    sc <- grid$scenario[i]; pr <- grid$prior[i]
    pf <- grid$p_fit[i];    pt <- grid$p_true[i]
    cat(sprintf("\n[%s] Cell %d/%d: scenario=%s | prior=%s | P_fit=%d | P_true=%d\n",
                format(Sys.time(), "%H:%M:%S"), i, n_cells, sc, pr, pf, pt))
    tryCatch(
      run_simulation(
        scenario   = sc, n_rep = n_rep, n_sim = n_sim, k_fold = k_fold,
        p_fit      = pf, p_true = pt,   prior  = pr,   theta  = theta,
        n_cores    = n_cores, eval_times = eval_times, base_seed = base_seed
      ),
      error = function(e) message("  [ERROR] ", e$message)
    )
  }
  invisible(NULL)
}

# =============================================================================
# 2b.  Legacy file migration helper
# =============================================================================

#' Rename legacy replicate files to include the prior in their filename.
#'
#' If you ran ctbn_map_simulation.R before this path-fix was applied, your
#' replicates are named rep_<scenario>_<NNN>.rds (no prior tag).  Call this
#' function once to rename them to rep_<scenario>_<prior>_<NNN>.rds so that
#' run_all_scenarios_all_priors() can distinguish cells.
#'
#' Only files matching rep_<scenario>_<NNN>.rds that do NOT already have a
#' prior tag are renamed.  Existing prior-tagged files are left untouched.
#'
#' @param prior     character, the prior that was used for the legacy runs
#'                  (e.g. "spike_slab")
#' @param scenarios character vector, scenarios to rename (default: all three)
#' @param n_rep     integer, expected number of replicates
#' @param dry_run   logical; if TRUE, print what would be renamed without doing it
rename_legacy_replicates <- function(prior,
                                      scenarios = SCENARIO_META$scenario,
                                      n_rep     = N_REP,
                                      dry_run   = TRUE) {
  rep_dir <- "sim_results/replicates"
  renamed <- 0L

  for (sc in scenarios) {
    for (r in seq_len(n_rep)) {
      old_path <- file.path(rep_dir, sprintf("rep_%s_%03d.rds", sc, r))
      new_path <- file.path(rep_dir, sprintf("rep_%s_%s_%03d.rds", sc, prior, r))

      if (!file.exists(old_path)) next
      if (file.exists(new_path)) {
        cat(sprintf("  SKIP  %s (target already exists)\n",
                    basename(old_path)))
        next
      }

      if (dry_run) {
        cat(sprintf("  [DRY] %s -> %s\n",
                    basename(old_path), basename(new_path)))
      } else {
        file.rename(old_path, new_path)
        cat(sprintf("  RENAMED %s -> %s\n",
                    basename(old_path), basename(new_path)))
        renamed <- renamed + 1L
      }
    }
  }

  if (dry_run) {
    cat(sprintf("\n  dry_run=TRUE: no files were renamed.\n"))
    cat(sprintf("  Re-run with dry_run=FALSE to apply.\n"))
  } else {
    cat(sprintf("\n  Renamed %d files.\n", renamed))
  }
  invisible(renamed)
}

#' Load replicate RDS files for one (scenario, prior) cell.
#'
#' Naming convention (new):   rep_<scenario>_<prior>_<NNN>.rds
#' Naming convention (legacy): rep_<scenario>_<NNN>.rds
#' Falls back to legacy if new-style files are not found.
load_scenario_prior <- function(scenario, n_rep=N_REP, prior=NULL) {
  rep_dir <- "sim_results/replicates"

  if (!is.null(prior)) {
    paths_new <- file.path(rep_dir,
      sprintf("rep_%s_%s_%03d.rds", scenario, prior, seq_len(n_rep)))
    paths_old <- file.path(rep_dir,
      sprintf("rep_%s_%03d.rds", scenario, seq_len(n_rep)))
    paths <- ifelse(file.exists(paths_new), paths_new, paths_old)
  } else {
    all_files <- list.files(rep_dir,
      pattern = sprintf("^rep_%s_.*\\.rds$", scenario), full.names=TRUE)
    if (!length(all_files)) {
      message(sprintf("  No replicates found for scenario %s", scenario))
      return(list())
    }
    paths <- all_files
  }

  exist <- which(file.exists(paths))
  if (!length(exist)) return(list())
  cat(sprintf("  Loading %d/%d replicates for scenario=%s prior=%s\n",
              length(exist), length(paths), scenario,
              if (is.null(prior)) "ALL" else prior))
  lapply(paths[exist], readRDS)
}

#' Load and combine all (scenario × prior) replicates into one flat list.
#'
#' Attaches metadata slots (scenario, prior, p_fit, p_true, theta) to any
#' replicate objects that lack them (legacy files that predate prior-tagging).
load_all_scenarios <- function(scenarios = SCENARIO_META$scenario,
                                priors    = PRIOR_META$prior,
                                n_rep     = N_REP) {
  all_reps <- list()
  total    <- length(scenarios) * length(priors)
  i        <- 0L

  for (sc in scenarios) {
    meta <- SCENARIO_META[scenario == sc]
    for (pr in priors) {
      i <- i + 1L
      print_progress(i, total, sprintf("%s / %s", sc, pr))
      reps <- load_scenario_prior(sc, n_rep, pr)
      if (!length(reps)) next

      # Stamp metadata onto each replicate where absent
      reps <- lapply(reps, function(x) {
        if (is.null(x$prior))    x$prior    <- pr
        if (is.null(x$scenario)) x$scenario <- sc
        if (is.null(x$p_fit))    x$p_fit    <- meta$p_fit
        if (is.null(x$p_true))   x$p_true   <- meta$p_true
        if (is.null(x$theta))    x$theta    <- meta$theta
        x
      })
      all_reps <- c(all_reps, reps)
    }
  }
  cat("\n")
  results_log(sprintf("Loaded %d replicates across %d scenario x prior cells",
                       length(all_reps), total))
  all_reps
}

# =============================================================================
# 4.  Metric stacking — correct prior propagation
# =============================================================================

#' Stack one named metric slot from all replicates into a single data.table.
#'
#' Run-level metadata (scenario, prior, p_fit, p_true, theta, rep) is
#' ALWAYS stamped onto every row, overwriting any values already in the
#' sub-table.  This avoids the "prior repeats as first value" bug that
#' occurs when the sub-table already has a scenario/theta column from the
#' fold loop but lacks prior, causing set() to only add prior once and
#' recycle.  Explicit overwrite guarantees every row is correctly labelled.
stack_metric <- function(raw_list, field) {
  rbindlist(lapply(seq_along(raw_list), function(i) {
    x <- raw_list[[i]]
    d <- x[[field]]
    if (is.null(d) || nrow(d) == 0L) return(NULL)
    d <- copy(d)
    # Always overwrite metadata — do NOT use the "if absent" guard,
    # because sub-tables may carry stale scenario/theta from the fold loop
    # while lacking prior entirely, causing recycling of the first value.
    set(d, j = "rep",      value = x[["rep"]]      %||% i)
    set(d, j = "scenario", value = x[["scenario"]] %||% NA_character_)
    set(d, j = "prior",    value = x[["prior"]]    %||% NA_character_)
    set(d, j = "p_fit",    value = x[["p_fit"]]    %||% NA_integer_)
    set(d, j = "p_true",   value = x[["p_true"]]   %||% NA_integer_)
    set(d, j = "theta",    value = x[["theta"]]    %||% NA_real_)
    d
  }), fill = TRUE)
}

# =============================================================================
# 5.  Shared theme and factor helpers
# =============================================================================

theme_ctbn <- function(base_size = 11) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text       = element_text(face = "bold", size = rel(0.9)),
      legend.position  = "bottom",
      legend.key.size  = unit(0.45, "cm"),
      plot.title       = element_text(face = "bold", size = rel(1.05)),
      plot.subtitle    = element_text(size = rel(0.88), colour = "grey35"),
      axis.title       = element_text(size = rel(0.92))
    )
}

as_prior_factor <- function(x)
  factor(x, levels = PRIOR_LEVELS, labels = PRIOR_LABELS[PRIOR_LEVELS])

as_scenario_factor <- function(x)
  factor(x, levels = SCENARIO_META$scenario,
         labels = SCENARIO_META$short_label)

as_ptype_factor <- function(x)
  factor(x, levels = PTYPE_ORDER, labels = PTYPE_LABELS[PTYPE_ORDER])

# =============================================================================
# 6.  Main build function
# =============================================================================

#' Aggregate all replicates and produce every table and figure.
#'
#' @param raw_list   list from load_all_scenarios(), or NULL to reload
#' @param scenarios  character vector
#' @param priors     character vector
#' @param n_rep      integer
build_all_results <- function(raw_list  = NULL,
                               scenarios = SCENARIO_META$scenario,
                               priors    = PRIOR_META$prior,
                               n_rep     = N_REP) {

  if (is.null(raw_list)) {
    raw_list <- section("Loading replicates",
                        load_all_scenarios(scenarios, priors, n_rep))
  }
  if (!length(raw_list))
    stop("No replicates loaded. Run run_all_scenarios_all_priors() first.")

  # ── Stack ─────────────────────────────────────────────────────────────────
  recovery_all  <- section("Stacking recovery_dt",  stack_metric(raw_list, "recovery_dt"))
  selection_all <- section("Stacking selection_dt", stack_metric(raw_list, "selection_dt"))
  pred_all      <- section("Stacking pred_dt",      stack_metric(raw_list, "pred_dt"))
  oracle_all    <- section("Stacking oracle_dt",    stack_metric(raw_list, "oracle_dt"))

  # ── Diagnostic: verify prior distribution ─────────────────────────────────
  if (!is.null(recovery_all) && nrow(recovery_all)) {
    prior_counts <- recovery_all[, .N, by = .(scenario, prior)]
    results_log("Prior distribution in recovery_all:")
    print(prior_counts)
  }

  # ── Add ordered factor labels (on the actual objects, not loop copies) ───
  add_labels <- function(dt) {
    if (is.null(dt) || !nrow(dt)) return(dt)
    if ("scenario" %in% names(dt)) dt[, scenario_f := as_scenario_factor(scenario)]
    if ("prior"    %in% names(dt)) dt[, prior_f    := as_prior_factor(prior)]
    if ("param_type" %in% names(dt))
      dt[, ptype_f := as_ptype_factor(param_type)]
    dt
  }
  recovery_all  <- add_labels(recovery_all)
  selection_all <- add_labels(selection_all)
  pred_all      <- add_labels(pred_all)
  oracle_all    <- add_labels(oracle_all)

  # ── Detect AUC column ─────────────────────────────────────────────────────
  # compute_pred_metrics() stores AUC as a column named by the eval time
  # (e.g. "1" for EVAL_TIMES=1L) via as.list(av).  It may also exist as
  # "tdAUC" in older output.  Find whichever is present.
  auc_col <- NULL
  if (!is.null(pred_all) && nrow(pred_all)) {
    non_meta <- setdiff(names(pred_all),
                        c("condition","pll","brier","rep","fold","scenario",
                          "prior","p_fit","p_true","theta",
                          "scenario_f","prior_f"))
    if ("tdAUC" %in% non_meta)       auc_col <- "tdAUC"
    else if (length(non_meta))        auc_col <- non_meta[1L]
    if (!is.null(auc_col))
      results_log(sprintf("AUC column detected: '%s'", auc_col))
    else
      results_log("No AUC column found in pred_all", "WARN")
  }

  # ── Run each section ───────────────────────────────────────────────────────
  section("S0. Scenario key table",
          make_scenario_key_table(scenarios))

  section("S1. RMSE, bias and coverage (primary recovery figure)",
          make_recovery_figures_tables(recovery_all))

  section("S2. Covariate coefficient recovery",
          make_covariate_recovery(recovery_all))

  section("S3. Variable selection: TPR / FPR",
          make_selection_figures_tables(selection_all))

  section("S4. Predictive metrics: Poisson LL / Brier",
          make_pred_figures_tables(pred_all, oracle_all, auc_col))

  section("S5. Master comparison table",
          make_master_table(recovery_all, selection_all, pred_all, auc_col))

  section("S6. Prior ranking heat-map",
          make_prior_ranking_heatmap(recovery_all, selection_all, pred_all, auc_col))

  section("Final manifest",
    cat("\n  All outputs in sim_results/{figures,tables}/\n"))

  invisible(list(recovery_all  = recovery_all,
                 selection_all = selection_all,
                 pred_all      = pred_all,
                 oracle_all    = oracle_all))
}

# =============================================================================
# 7.  S0 — Scenario key table
# =============================================================================

make_scenario_key_table <- function(scenarios = SCENARIO_META$scenario) {
  tbl <- SCENARIO_META[scenario %in% scenarios, .(
    Scenario                     = scenario,
    `$P_{\\text{fit}}$`          = p_fit,
    `$P_{\\text{true}}$`         = p_true,
    `$\\theta$`                  = theta,
    Description                  = gsub("\n"," ",spec_label)
  )]
  save_sim_tex(
    xtable(tbl,
      caption = sprintf(
        paste0("Simulation scenario definitions ($N=%s$ subjects, ",
               "$R=%d$ replicates, $K=%d$ CV folds, $\\theta=1.0$)."),
        formatC(N_SIM, big.mark=","), N_REP, K_FOLD),
      label = "tab:scenarios"),
    "tab_s0_scenario_key",
    sanitize.text.function = identity
  )
}

# =============================================================================
# 8.  S1 — Parameter recovery: RMSE, bias, and 95% CI coverage
# =============================================================================
#
# Three parameter types shown as facets in each figure:
#   1. Main effects     — β for condition influencers
#   2. 2-way interactions — β for pairwise condition products
#   3. Covariates       — γ (age, sex, smoking)
#
# A1 IMPUTATION FOR 2-WAY INTERACTIONS
# ─────────────────────────────────────
# Scenario A1 (P_fit=1) fits no interaction terms, but the true DGP has
# P_true=2.  The model implicitly sets every interaction to exactly 0.
# We formalise this:
#   est_mean = 0,  bias = 0 − true_val,  sq_err = true_val^2,  covered = 0
# Imputed rows are drawn from A2's interaction parameter list (same DGP).
# They are shown with dashed lines / open points to distinguish them from
# estimated rows.

make_recovery_figures_tables <- function(recovery_all) {

  if (is.null(recovery_all) || !nrow(recovery_all)) {
    warning("recovery_all is empty — skipping recovery section.")
    return(NULL)
  }

  FACET3 <- c("main_effect", "interaction_2way", "covariate")

  # ── 2-way: real rows + imputed A1 rows ────────────────────────────────────
  rec_2way_real <- recovery_all[param_type == "interaction_2way" & !is.na(sq_err)]

  a1_imputed <- data.table()
  if (nrow(rec_2way_real) > 0) {
    template <- recovery_all[param_type == "interaction_2way" &
                               scenario == "A2" & !is.na(true_val),
                             .(rep, prior, condition, parameter, true_val)]
    if (nrow(template) > 0) {
      a1_imputed <- copy(template)
      a1_imputed[, `:=`(
        scenario   = "A1",  p_fit = 1L, p_true = 2L, theta = 1.0,
        param_type = "interaction_2way",
        est_mean   = 0.0,   est_var = 0.0,
        ci_lo      = 0.0,   ci_hi   = 0.0,
        bias       = 0.0 - true_val,
        sq_err     = true_val^2,
        covered    = 0.0,
        pip        = NA_real_
      )]
    }
  }

  rec_all3 <- rbindlist(
    list(recovery_all[param_type %in% FACET3 & !is.na(sq_err)],
         a1_imputed),
    fill = TRUE
  )
  # Re-label factors after rbindlist (drops them)
  rec_all3[, ptype_f    := as_ptype_factor(param_type)]
  rec_all3[, scenario_f := as_scenario_factor(scenario)]
  rec_all3[, prior_f    := as_prior_factor(prior)]
  rec_all3 <- rec_all3[param_type %in% FACET3]

  # Flag imputed rows for visual distinction
  rec_all3[, imputed := (scenario == "A1" & param_type == "interaction_2way")]

  # ── Aggregate ─────────────────────────────────────────────────────────────
  rec_agg <- rec_all3[, .(
    N        = .N,
    rmse     = sqrt(mean(sq_err,  na.rm = TRUE)),
    se_rmse  = sd(sqrt(sq_err),   na.rm = TRUE) / sqrt(.N),
    bias     = mean(bias,         na.rm = TRUE),
    se_bias  = sd(bias,           na.rm = TRUE) / sqrt(.N),
    coverage = mean(covered,      na.rm = TRUE),
    se_cov   = sqrt(mean(covered, na.rm = TRUE) *
                    (1 - mean(covered, na.rm = TRUE)) / .N),
    imputed  = any(imputed)
  ), by = .(scenario, prior, param_type, scenario_f, prior_f, ptype_f)]

  # Enforce 3-facet ordering
  facet_levels <- PTYPE_LABELS[FACET3]
  rec_agg[, ptype_f := factor(as.character(ptype_f), levels = facet_levels)]

  # ── Table S1 ──────────────────────────────────────────────────────────────
  tbl_s1 <- rec_agg[order(scenario, prior_f, ptype_f), .(
    Scenario           = scenario,
    Prior              = as.character(prior_f),
    `Param. type`      = as.character(ptype_f),
    `$N$`              = N,
    RMSE               = sprintf("%.4f",         rmse),
    `Bias (SE)`        = sprintf("%.4f (%.4f)", bias,     se_bias),
    `Coverage (95\\%)` = sprintf("%.3f (%.3f)", coverage, se_cov)
  )]
  save_sim_tex(
    xtable(tbl_s1,
      caption = sprintf(
        paste0("MAP parameter recovery ($\\theta=1.0$, $R=%d$, $N=%s$). ",
               "For scenario A1 (P$_{\\text{fit}}$=1), 2-way interaction ",
               "estimates are imputed as 0 (correct by design). ",
               "Nominal 95\\%% CI coverage = 0.95."),
        N_REP, formatC(N_SIM, big.mark = ",")),
      label = "tab:recovery"),
    "tab_s1_recovery_summary",
    sanitize.text.function = identity
  )

  # ── Figure helpers ────────────────────────────────────────────────────────
  # ggplot2 cannot draw a geom_line() with varying linetype AND varying colour
  # along the same group.  The solution is two separate geom_* layers:
  #   Layer 1: estimated rows (imputed == FALSE) — solid lines, filled points
  #   Layer 2: imputed A1 interaction rows       — dashed lines, open points
  # Both layers share the same colour scale so the legend is unified.
  # `show.legend` on the second layer is FALSE to avoid duplicate legend keys.

  rec_est <- rec_agg[imputed == FALSE]
  rec_imp <- rec_agg[imputed == TRUE]

  # Reusable function to build each of the three figures
  make_rec_fig <- function(y_var, ymin_var, ymax_var,
                            y_lab, title_str, subtitle_str,
                            extra_layers = list()) {

    p <- ggplot(mapping = aes(x = scenario_f, colour = prior_f,
                               group = prior_f)) +
      # ── estimated rows: solid lines, filled circles ──
      geom_line(data = rec_est,
                aes(y = .data[[y_var]]),
                linewidth = 0.9, linetype = "solid") +
      geom_point(data = rec_est,
                 aes(y = .data[[y_var]]),
                 shape = 16, size = 2.5) +
      geom_errorbar(data = rec_est,
                    aes(ymin = .data[[ymin_var]],
                        ymax = .data[[ymax_var]]),
                    width = 0.15, linewidth = 0.6) +
      # ── imputed A1 rows: dashed lines, open circles ──
      geom_line(data = rec_imp,
                aes(y = .data[[y_var]]),
                linewidth = 0.75, linetype = "dashed",
                show.legend = FALSE) +
      geom_point(data = rec_imp,
                 aes(y = .data[[y_var]]),
                 shape = 1, size = 2.5,
                 show.legend = FALSE) +
      geom_errorbar(data = rec_imp,
                    aes(ymin = .data[[ymin_var]],
                        ymax = .data[[ymax_var]]),
                    width = 0.15, linewidth = 0.5, linetype = "dashed",
                    show.legend = FALSE) +
      facet_wrap(~ ptype_f, scales = "free_y", nrow = 1) +
      scale_colour_manual(values = PRIOR_PALETTE, name = "Prior") +
      labs(title    = title_str,
           subtitle = subtitle_str,
           x = "Scenario", y = y_lab) +
      theme_ctbn() +
      theme(axis.text.x = element_text(angle = 20, hjust = 1))

    # Add any caller-supplied extra layers (e.g. geom_hline, scale_y_*)
    for (lyr in extra_layers) p <- p + lyr
    p
  }

  # ── Pre-compute ymin/ymax columns for all three figures ───────────────────
  rec_agg[, rmse_lo := pmax(rmse - 1.96 * se_rmse, 0)]
  rec_agg[, rmse_hi := rmse + 1.96 * se_rmse]
  rec_agg[, bias_lo := bias - 1.96 * se_bias]
  rec_agg[, bias_hi := bias + 1.96 * se_bias]
  rec_agg[, cov_lo  := coverage - 1.96 * se_cov]
  rec_agg[, cov_hi  := coverage + 1.96 * se_cov]

  # Refresh split tables after adding columns
  rec_est <- rec_agg[imputed == FALSE]
  rec_imp <- rec_agg[imputed == TRUE]

  # ── Figure 1a: RMSE (3 facets) ────────────────────────────────────────────
  fig1a <- make_rec_fig(
    y_var    = "rmse",
    ymin_var = "rmse_lo",
    ymax_var = "rmse_hi",
    y_lab    = "Root mean squared error",
    title_str    = "Figure a \u2014 RMSE by parameter type, scenario, and prior",
    subtitle_str = paste0("Lower is better. ",
                          "Dashed/open = A1 interaction imputed as 0. ",
                          "Bars = \u00b11.96 SE.")
  )
  save_sim_fig(fig1a, "fig1a_rmse", w = 16, h = 5)

  # ── Figure 1b: Bias (3 facets) ────────────────────────────────────────────
  fig1b <- make_rec_fig(
    y_var    = "bias",
    ymin_var = "bias_lo",
    ymax_var = "bias_hi",
    y_lab    = "Monte Carlo bias",
    title_str    = "Figure b \u2014 Bias by parameter type, scenario, and prior",
    subtitle_str = "Horizontal dashed = zero bias. Bars = \u00b11.96 SE.",
    extra_layers = list(
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50")
    )
  )
  save_sim_fig(fig1b, "fig1b_bias", w = 16, h = 5)

  # ── Figure 2: 95% CI coverage (3 facets) ─────────────────────────────────
  fig2 <- make_rec_fig(
    y_var    = "coverage",
    ymin_var = "cov_lo",
    ymax_var = "cov_hi",
    y_lab    = "Empirical coverage",
    title_str    = "Figure c \u2014 95% CI coverage by parameter type and scenario",
    subtitle_str = "Dashed = nominal 0.95. A1 interaction coverage = 0 (imputed; no CI).",
    extra_layers = list(
      geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey40"),
      scale_y_continuous(labels = percent_format(accuracy = 1),
                         limits = c(0, 1.02))
    )
  )
  save_sim_fig(fig2, "fig2_coverage", w = 16, h = 5)

  # ── Combined paper figure ─────────────────────────────────────────────────
  save_sim_fig(fig1a / fig1b / fig2, "fig1_recovery_combined", w = 16, h = 14)

  invisible(rec_agg)
}

# =============================================================================
# 9.  S2 — Covariate coefficient recovery
# =============================================================================

make_covariate_recovery <- function(recovery_all) {

  if (is.null(recovery_all) || !nrow(recovery_all)) return(NULL)

  # recovery_all uses column `parameter` (e.g. "g_age") not "influencer"
  # Extract covariate name from parameter string: "g_age" -> "age"
  cov_rec <- recovery_all[param_type == "covariate" & !is.na(sq_err)]
  if (!nrow(cov_rec)) {
    warning("No covariate rows in recovery_all."); return(NULL)
  }

  cov_rec[, covariate := sub("^g_","", parameter)]

  cov_agg <- cov_rec[, .(
    N        = .N,
    rmse     = sqrt(mean(sq_err,  na.rm=TRUE)),
    bias     = mean(bias,         na.rm=TRUE),
    se_bias  = sd(bias,           na.rm=TRUE) / sqrt(.N),
    coverage = mean(covered,      na.rm=TRUE),
    se_cov   = sqrt(mean(covered,na.rm=TRUE)*
                    (1-mean(covered,na.rm=TRUE))/.N)
  ), by = .(scenario, prior, covariate, scenario_f, prior_f)]

  # ── Table S2 ──────────────────────────────────────────────────────────────
  tbl_s2 <- cov_agg[order(scenario, prior, covariate), .(
    Scenario           = scenario,
    Prior              = as.character(prior_f),
    Covariate          = covariate,
    `$N$`              = N,
    RMSE               = sprintf("%.4f",         rmse),
    `Bias (SE)`        = sprintf("%.4f (%.4f)", bias,     se_bias),
    `Coverage (95\\%)` = sprintf("%.3f (%.3f)", coverage, se_cov)
  )]
  save_sim_tex(
    xtable(tbl_s2,
      caption = sprintf(
        paste0("Covariate coefficient recovery ($\\theta=1.0$, $R=%d$). ",
               "True values taken from TRUE\\_GAMMA."),
        N_REP),
      label = "tab:cov_recovery"),
    "tab_s2_covariate_recovery",
    sanitize.text.function = identity
  )

  # ── Figure 3: Covariate PIP / pseudo-PIP ──────────────────────────────────
  # Use the `pip` column from recovery_all (stored for covariate rows)
  cov_pip <- recovery_all[param_type == "covariate" & !is.na(pip)]
  if (nrow(cov_pip)) {
    cov_pip[, covariate := sub("^g_","", parameter)]
    cov_pip_agg <- cov_pip[, .(
      mean_pip     = mean(pip,      na.rm=TRUE),
      se_pip       = sd(pip,        na.rm=TRUE) / sqrt(.N),
      pct_selected = mean(pip >= 0.5, na.rm=TRUE)
    ), by = .(scenario, prior, covariate, scenario_f, prior_f)]

    fig3 <- ggplot(cov_pip_agg,
      aes(x=covariate, y=mean_pip,
          ymin=pmax(mean_pip - 1.96*se_pip,0),
          ymax=pmin(mean_pip + 1.96*se_pip,1),
          colour=prior_f, group=prior_f)) +
      geom_hline(yintercept=0.5, linetype="dashed", colour="grey40") +
      geom_line(linewidth=0.8, alpha=0.8) +
      geom_point(size=2.2) +
      geom_errorbar(width=0.2, linewidth=0.5) +
      facet_wrap(~ scenario_f, nrow=1) +
      scale_colour_manual(values=PRIOR_PALETTE, name="Prior") +
      scale_y_continuous(labels=percent_format(accuracy=1), limits=c(0,1)) +
      labs(
        #title    = "Figure 3 \u2014 Covariate PIP / pseudo-PIP by scenario and prior",
        #subtitle = paste0("Dashed = 0.50 selection threshold. ",
        #                  "For horseshoe: pseudo-PIP = 1\u2212\u03ba."),
        x = "Covariate", y = "Mean PIP / pseudo-PIP"
      ) +
      theme_ctbn() +
      theme(axis.text.x = element_text(angle=30, hjust=1))
    save_sim_fig(fig3, "fig3_covariate_pip", w=16, h=6)
  }

  invisible(cov_agg)
}

# =============================================================================
# 10.  S3 — Variable selection: TPR, FPR, selection AUC
# =============================================================================

make_selection_figures_tables <- function(selection_all) {

  if (is.null(selection_all) || !nrow(selection_all)) {
    warning("selection_all is empty — skipping selection section.")
    return(NULL)
  }

  # selection_all columns: level, tpr, fpr, selection_auc,
  #   rep, scenario, prior, p_fit, p_true, theta, scenario_f, prior_f

  sel_agg <- selection_all[, .(
    N        = .N,
    mean_tpr = mean(tpr,           na.rm=TRUE),
    se_tpr   = sd(tpr,             na.rm=TRUE) / sqrt(.N),
    mean_fpr = mean(fpr,           na.rm=TRUE),
    se_fpr   = sd(fpr,             na.rm=TRUE) / sqrt(.N),
    mean_sauc= mean(selection_auc, na.rm=TRUE),
    se_sauc  = sd(selection_auc,   na.rm=TRUE) / sqrt(.N)
  ), by = .(scenario, prior, level, scenario_f, prior_f)]

  # ── Table S3 ──────────────────────────────────────────────────────────────
  tbl_s3 <- sel_agg[order(scenario, level, prior), .(
    Scenario           = scenario,
    Prior              = as.character(prior_f),
    Level              = level,
    `$N$ reps`         = N,
    `TPR (SE)`         = sprintf("%.3f (%.3f)", mean_tpr,  se_tpr),
    `FPR (SE)`         = sprintf("%.3f (%.3f)", mean_fpr,  se_fpr),
    `Sel. AUC (SE)`    = sprintf("%.3f (%.3f)", mean_sauc, se_sauc)
  )]
  save_sim_tex(
    xtable(tbl_s3,
      caption = sprintf(
        paste0("Variable selection metrics ($\\theta=1.0$, ",
               "$R=%d$, PIP threshold $\\geq 0.50$). ",
               "For horseshoe prior, pseudo-PIP $= 1 - \\kappa$ is used. ",
               "Sel.\\ AUC = area under the ROC curve for the PIP score ",
               "vs true non-zero indicator."),
        N_REP),
      label = "tab:selection"),
    "tab_s3_selection_metrics",
    sanitize.text.function = identity
  )

  # ── Figure 4: TPR and FPR side-by-side, scenario on x, prior colours ─────
  sel_main_2way <- sel_agg[level %in% c("main_effect","interaction_2way")]
  sel_main_2way[, level_label := fcase(
    level == "main_effect",      "Main effects",
    level == "interaction_2way", "2-way interactions",
    default = level)]

  fig4a <- ggplot(sel_main_2way,
    aes(x=scenario_f, y=mean_tpr,
        ymin=pmax(mean_tpr-1.96*se_tpr,0),
        ymax=pmin(mean_tpr+1.96*se_tpr,1),
        colour=prior_f, group=prior_f)) +
    geom_line(linewidth=0.9) + geom_point(size=2.5) +
    geom_errorbar(width=0.15, linewidth=0.6) +
    facet_wrap(~ level_label, nrow=1) +
    scale_colour_manual(values=PRIOR_PALETTE, name="Prior") +
    scale_y_continuous(labels=percent_format(accuracy=1), limits=c(0,1)) +
    labs(title="Figure a \u2014 True positive rate (sensitivity) by scenario and prior",
         #subtitle="Higher is better. PIP/pseudo-PIP threshold = 0.50.",
         x="Scenario", y="Mean TPR") +
    theme_ctbn() +
    theme(axis.text.x=element_text(angle=20, hjust=1))

  fig4b <- ggplot(sel_main_2way,
    aes(x=scenario_f, y=mean_fpr,
        ymin=pmax(mean_fpr-1.96*se_fpr,0),
        ymax=pmin(mean_fpr+1.96*se_fpr,1),
        colour=prior_f, group=prior_f)) +
    geom_line(linewidth=0.9) + geom_point(size=2.5) +
    geom_errorbar(width=0.15, linewidth=0.6) +
    facet_wrap(~ level_label, nrow=1) +
    scale_colour_manual(values=PRIOR_PALETTE, name="Prior") +
    scale_y_continuous(labels=percent_format(accuracy=1), limits=c(0,1)) +
    labs(title="Figure b \u2014 False positive rate (1\u2212specificity) by scenario and prior",
         #subtitle="Lower is better. PIP/pseudo-PIP threshold = 0.50.",
         x="Scenario", y="Mean FPR") +
    theme_ctbn() +
    theme(axis.text.x=element_text(angle=20, hjust=1))

  save_sim_fig(fig4a / fig4b, "fig4_selection_tpr_fpr", w=14, h=10)

  # ── Selection AUC figure ──────────────────────────────────────────────────
  fig4c <- ggplot(sel_main_2way,
    aes(x=scenario_f, y=mean_sauc,
        ymin=pmax(mean_sauc-1.96*se_sauc,0),
        ymax=pmin(mean_sauc+1.96*se_sauc,1),
        colour=prior_f, group=prior_f)) +
    geom_hline(yintercept=0.5, linetype="dashed", colour="grey50") +
    geom_line(linewidth=0.9) + geom_point(size=2.5) +
    geom_errorbar(width=0.15, linewidth=0.6) +
    facet_wrap(~ level_label, nrow=1) +
    scale_colour_manual(values=PRIOR_PALETTE, name="Prior") +
    scale_y_continuous(labels=percent_format(accuracy=1), limits=c(0.4,1)) +
    labs(title="Figure c \u2014 Selection AUC (PIP score vs true non-zero)",
         subtitle="Dashed = 0.5 (random). Higher is better.",
         x="Scenario", y="Mean selection AUC") +
    theme_ctbn() +
    theme(axis.text.x=element_text(angle=20, hjust=1))
  save_sim_fig(fig4c, "fig4c_selection_auc", w=14, h=5)

  assign("sel_agg", sel_agg, envir=.GlobalEnv)
  invisible(sel_agg)
}

# =============================================================================
# 11.  S4 — Predictive metrics: Poisson LL and Brier score
# =============================================================================

make_pred_figures_tables <- function(pred_all, oracle_all, auc_col=NULL) {

  if (is.null(pred_all) || !nrow(pred_all)) {
    warning("pred_all is empty — skipping predictive section.")
    return(NULL)
  }

  # pred_all columns: condition, pll, brier, [auc_col], rep, fold,
  #   scenario, prior, p_fit, p_true, theta, scenario_f, prior_f

  # Macro-average pll and brier across conditions per rep
  pred_rep <- pred_all[, .(
    pll   = mean(pll,   na.rm=TRUE),
    brier = mean(brier, na.rm=TRUE)
  ), by = .(rep, scenario, prior, p_fit, p_true, theta, scenario_f, prior_f)]

  pred_agg <- pred_rep[, .(
    mean_pll   = mean(pll,   na.rm=TRUE),
    se_pll     = sd(pll,     na.rm=TRUE) / sqrt(.N),
    mean_brier = mean(brier, na.rm=TRUE),
    se_brier   = sd(brier,   na.rm=TRUE) / sqrt(.N)
  ), by = .(scenario, prior, scenario_f, prior_f)]

  # Per-condition aggregation
  pred_cond_agg <- pred_all[, .(
    mean_pll   = mean(pll,   na.rm=TRUE),
    se_pll     = sd(pll,     na.rm=TRUE) / sqrt(.N),
    mean_brier = mean(brier, na.rm=TRUE),
    se_brier   = sd(brier,   na.rm=TRUE) / sqrt(.N)
  ), by = .(scenario, prior, condition, scenario_f, prior_f)]

  # ── Table S4 ──────────────────────────────────────────────────────────────
  tbl_s4 <- pred_agg[order(scenario, prior), .(
    Scenario    = scenario,
    Prior       = as.character(prior_f),
    `Pois. LL (SE)` = sprintf("%.4f (%.4f)", mean_pll,   se_pll),
    `Brier (SE)`    = sprintf("%.4f (%.4f)", mean_brier, se_brier)
  )]
  save_sim_tex(
    xtable(tbl_s4,
      caption = sprintf(
        paste0("Macro-averaged predictive accuracy ($K=%d$ CV, $\\theta=1.0$, $R=%d$). ",
               "Higher Poisson LL is better; lower Brier is better."),
        K_FOLD, N_REP),
      label = "tab:pred_metrics"),
    "tab_s4_pred_metrics",
    sanitize.text.function = identity
  )

  # ── Figure 5a: Poisson LL by scenario × prior ─────────────────────────────
  fig5a <- ggplot(pred_agg,
    aes(x=scenario_f, y=mean_pll,
        ymin=mean_pll-1.96*se_pll, ymax=mean_pll+1.96*se_pll,
        colour=prior_f, group=prior_f)) +
    geom_line(linewidth=0.9) + geom_point(size=2.5) +
    geom_errorbar(width=0.15, linewidth=0.6) +
    scale_colour_manual(values=PRIOR_PALETTE, name="Prior") +
    labs(#title="Figure 5a \u2014 Macro-averaged Poisson log-likelihood by scenario and prior",
         #subtitle=sprintf("Higher is better. K=%d CV, mean \u00b11.96 SE across %d replicates.",
         #                 K_FOLD, N_REP),
         x="Scenario", y="Mean Poisson log-likelihood") +
    theme_ctbn() +
    theme(axis.text.x=element_text(angle=20, hjust=1))

  # ── Figure 5b: Brier score ────────────────────────────────────────────────
  fig5b <- ggplot(pred_agg,
    aes(x=scenario_f, y=mean_brier,
        ymin=mean_brier-1.96*se_brier, ymax=mean_brier+1.96*se_brier,
        colour=prior_f, group=prior_f)) +
    geom_line(linewidth=0.9) + geom_point(size=2.5) +
    geom_errorbar(width=0.15, linewidth=0.6) +
    scale_colour_manual(values=PRIOR_PALETTE, name="Prior") +
    labs(#title="Figure 5b \u2014 Macro-averaged Brier score by scenario and prior",
         #subtitle="Lower is better.",
         x="Scenario", y="Mean Brier score") +
    theme_ctbn() +
    theme(axis.text.x=element_text(angle=20, hjust=1))

  save_sim_fig(fig5a / fig5b, "fig5_pred_metrics", w=12, h=10)

  # ── Figure 5c: per-condition LL, faceted by scenario ─────────────────────
  fig5c <- ggplot(pred_cond_agg,
    aes(x=condition, y=mean_pll,
        ymin=mean_pll-1.96*se_pll, ymax=mean_pll+1.96*se_pll,
        colour=prior_f, group=prior_f)) +
    geom_line(linewidth=0.7, alpha=0.8) + geom_point(size=1.8) +
    geom_errorbar(width=0.25, linewidth=0.5, alpha=0.6) +
    facet_wrap(~ scenario_f, nrow=1) +
    scale_colour_manual(values=PRIOR_PALETTE, name="Prior") +
    labs(#title="Figure 5c \u2014 Per-condition Poisson LL by scenario and prior",
         x="Condition", y="Mean Poisson LL") +
    theme_ctbn(base_size=9) +
    theme(axis.text.x=element_text(angle=40, hjust=1))
  save_sim_fig(fig5c, "fig5c_pll_per_condition", w=16, h=5)

  # ── Oracle efficiency ──────────────────────────────────────────────────────
  if (!is.null(oracle_all) && nrow(oracle_all)) {
    oracle_rep <- oracle_all[, .(
      oracle_pll   = mean(oracle_pll,   na.rm=TRUE),
      oracle_brier = mean(oracle_brier, na.rm=TRUE)
    ), by = .(rep, condition, scenario, p_fit, p_true)]

    pred_merge <- pred_all[, .(
      pll   = mean(pll,   na.rm=TRUE),
      brier = mean(brier, na.rm=TRUE)
    ), by = .(rep, condition, scenario, prior, p_fit, p_true)]

    eff_dt <- merge(pred_merge, oracle_rep,
                    by=c("rep","condition","scenario","p_fit","p_true"),
                    all.x=TRUE)
    eff_dt[, brier_ratio := brier / oracle_brier]
    eff_dt[, pll_diff    := pll   - oracle_pll]
    eff_dt[, prior_f    := as_prior_factor(prior)]
    eff_dt[, scenario_f := as_scenario_factor(scenario)]

    eff_agg <- eff_dt[, .(
      mean_ratio = mean(brier_ratio, na.rm=TRUE),
      se_ratio   = sd(brier_ratio,   na.rm=TRUE) / sqrt(.N),
      mean_diff  = mean(pll_diff,    na.rm=TRUE),
      se_diff    = sd(pll_diff,      na.rm=TRUE) / sqrt(.N)
    ), by = .(scenario, prior, condition, scenario_f, prior_f)]

    fig5d <- ggplot(eff_dt,
      aes(x=brier_ratio, colour=prior_f, fill=prior_f)) +
      geom_density(alpha=0.2, linewidth=0.7) +
      geom_vline(xintercept=1, linetype="dashed", colour="grey40") +
      facet_grid(scenario_f ~ condition, scales="free_y") +
      scale_colour_manual(values=PRIOR_PALETTE, name="Prior") +
      scale_fill_manual(  values=PRIOR_PALETTE, name="Prior") +
      labs(#title="Figure 5d \u2014 Oracle Brier ratio by scenario, condition, and prior",
           #subtitle="Dashed = 1 (oracle). Values >1 indicate efficiency loss.",
           x="Brier ratio (fitted / oracle)", y="Density") +
      theme_ctbn(base_size=8)
    save_sim_fig(fig5d, "fig5d_oracle_efficiency", w=20, h=9)

    assign("eff_agg", eff_agg, envir=.GlobalEnv)
  }

  assign("pred_agg", pred_agg, envir=.GlobalEnv)
  invisible(pred_agg)
}

# =============================================================================
# 12.  S5 — Master comparison table
# =============================================================================

make_master_table <- function(recovery_all, selection_all, pred_all, auc_col=NULL) {

  rows <- list()

  # ── Recovery: macro RMSE and coverage (main effects only) ─────────────────
  if (!is.null(recovery_all) && nrow(recovery_all)) {
    rec_main <- recovery_all[param_type == "main_effect" & !is.na(sq_err), .(
      rmse     = sqrt(mean(sq_err,  na.rm=TRUE)),
      coverage = mean(covered,      na.rm=TRUE)
    ), by = .(scenario, prior)]
    rows[["recovery"]] <- rec_main
  }

  # ── Selection: main-effect TPR / FPR ──────────────────────────────────────
  if (!is.null(selection_all) && nrow(selection_all)) {
    sel_main <- selection_all[level == "main_effect", .(
      tpr  = mean(tpr,           na.rm=TRUE),
      fpr  = mean(fpr,           na.rm=TRUE),
      sauc = mean(selection_auc, na.rm=TRUE)
    ), by = .(scenario, prior)]
    rows[["selection"]] <- sel_main
  }

  # ── Predictive: macro Poisson LL and Brier ─────────────────────────────────
  if (!is.null(pred_all) && nrow(pred_all)) {
    pred_macro <- pred_all[, .(
      pll   = mean(pll,   na.rm=TRUE),
      brier = mean(brier, na.rm=TRUE)
    ), by = .(scenario, prior)]
    rows[["pred"]] <- pred_macro
  }

  if (!length(rows)) { warning("No data for master table."); return(NULL) }

  master <- Reduce(function(a, b) merge(a, b, by=c("scenario","prior"),
                                         all=TRUE), rows)
  master[, prior_f    := as_prior_factor(prior)]
  master[, scenario_f := as_scenario_factor(scenario)]

  tbl_m <- master[order(scenario, prior_f), .(
    Scenario            = scenario,
    Prior               = as.character(prior_f),
    `RMSE (main)`       = if (!is.null(master$rmse))
                            sprintf("%.4f", rmse)     else NA_character_,
    `Coverage (main)`   = if (!is.null(master$coverage))
                            sprintf("%.3f", coverage)  else NA_character_,
    `TPR (main)`        = if (!is.null(master$tpr))
                            sprintf("%.3f", tpr)       else NA_character_,
    `FPR (main)`        = if (!is.null(master$fpr))
                            sprintf("%.3f", fpr)       else NA_character_,
    `Sel. AUC (main)`   = if (!is.null(master$sauc))
                            sprintf("%.3f", sauc)      else NA_character_,
    `Pois. LL`          = if (!is.null(master$pll))
                            sprintf("%.4f", pll)       else NA_character_,
    `Brier`             = if (!is.null(master$brier))
                            sprintf("%.4f", brier)     else NA_character_
  )]
  # Drop all-NA columns
  all_na_cols <- names(tbl_m)[sapply(tbl_m, function(x) all(is.na(x)))]
  if (length(all_na_cols)) tbl_m[, (all_na_cols) := NULL]

  save_sim_tex(
    xtable(tbl_m,
      caption = sprintf(
        paste0("Master prior comparison ($\\theta=1.0$, $R=%d$, $N=%s$). ",
               "Recovery and selection from full-data fits; ",
               "predictive metrics from $K=%d$ CV. ",
               "RMSE and coverage for main-effect $\\beta$ coefficients only."),
        N_REP, formatC(N_SIM, big.mark=","), K_FOLD),
      label = "tab:master"),
    "tab_s5_master",
    sanitize.text.function = identity
  )

  invisible(master)
}

# =============================================================================
# 13.  S6 — Prior ranking heat-map (aggregates all metrics into ranks)
# =============================================================================

make_prior_ranking_heatmap <- function(recovery_all, selection_all,
                                        pred_all, auc_col=NULL) {

  metrics_list <- list()

  if (!is.null(recovery_all) && nrow(recovery_all)) {
    # RMSE: lower is better
    metrics_list[["RMSE\n(main, lower)"]] <- list(
      dt = recovery_all[param_type=="main_effect" & !is.na(sq_err),
                        .(val=sqrt(mean(sq_err, na.rm=TRUE))),
                        by=.(scenario,prior)],
      higher_better = FALSE)
    # Coverage: closer to 0.95 is better — use abs distance from 0.95
    metrics_list[["Cov. distance\nfrom 0.95 (lower)"]] <- list(
      dt = recovery_all[param_type=="main_effect" & !is.na(covered),
                        .(val=abs(mean(covered,na.rm=TRUE) - 0.95)),
                        by=.(scenario,prior)],
      higher_better = FALSE)
  }

  if (!is.null(selection_all) && nrow(selection_all)) {
    metrics_list[["TPR\n(main, higher)"]] <- list(
      dt = selection_all[level=="main_effect",
                         .(val=mean(tpr, na.rm=TRUE)),
                         by=.(scenario,prior)],
      higher_better = TRUE)
    metrics_list[["FPR\n(main, lower)"]] <- list(
      dt = selection_all[level=="main_effect",
                         .(val=mean(fpr, na.rm=TRUE)),
                         by=.(scenario,prior)],
      higher_better = FALSE)
    metrics_list[["Sel. AUC\n(main, higher)"]] <- list(
      dt = selection_all[level=="main_effect",
                         .(val=mean(selection_auc, na.rm=TRUE)),
                         by=.(scenario,prior)],
      higher_better = TRUE)
  }

  if (!is.null(pred_all) && nrow(pred_all)) {
    metrics_list[["Pois. LL\n(higher)"]] <- list(
      dt = pred_all[, .(val=mean(pll, na.rm=TRUE)), by=.(scenario,prior)],
      higher_better = TRUE)
    metrics_list[["Brier\n(lower)"]] <- list(
      dt = pred_all[, .(val=mean(brier, na.rm=TRUE)), by=.(scenario,prior)],
      higher_better = FALSE)
  }

  if (!length(metrics_list)) { warning("No data for ranking."); return(NULL) }

  rank_list <- lapply(names(metrics_list), function(mn) {
    info <- metrics_list[[mn]]
    d    <- copy(info$dt)
    hb   <- info$higher_better
    d[, rank := frank(if (hb) -val else val, ties.method="average"),
      by = scenario]
    d[, metric := mn]
    d
  })
  rank_all <- rbindlist(rank_list)
  rank_all[, prior_f    := as_prior_factor(prior)]
  rank_all[, scenario_f := as_scenario_factor(scenario)]
  rank_all[, metric_f   := factor(metric, levels=names(metrics_list))]

  n_priors <- length(PRIOR_LEVELS)

  fig6 <- ggplot(rank_all,
    aes(x=metric_f, y=prior_f, fill=rank)) +
    geom_tile(colour="white", linewidth=0.6) +
    geom_text(aes(label=sprintf("%.1f", rank)), size=3.5, fontface="bold") +
    facet_wrap(~ scenario_f, nrow=1) +
    scale_fill_gradient2(
      low      = "#2166AC",
      mid      = "#F7F7F7",
      high     = "#D6604D",
      midpoint = (n_priors + 1) / 2,
      name     = "Rank\n(1=best)",
      limits   = c(1, n_priors)
    ) +
    labs(
      #title    = "Figure 6 \u2014 Prior ranking heat-map by scenario and metric",
      #subtitle = paste0("Rank 1 (blue) = best prior on that metric in that scenario. ",
      #                  "Rank ", n_priors, " (red) = worst."),
      x = "Metric", y = "Prior"
    ) +
    theme_ctbn() +
    theme(legend.position = "right",
          axis.text.x     = element_text(size=8, angle=15, hjust=1))
  save_sim_fig(fig6, "fig6_prior_ranking_heatmap", w=16, h=5)

  invisible(rank_all)
}

# =============================================================================
# 14.  Entry points
# =============================================================================

#' Full pipeline: run simulations then build all results.
run_and_analyse <- function(n_cores=20L, n_rep=N_REP, n_sim=N_SIM, ...) {
  cat(sprintf("\n%s\n  CTBN MAP — Full Simulation Pipeline\n%s\n\n",
              strrep("=",70), strrep("=",70)))
  section("Step 1: Running simulations",
          run_all_scenarios_all_priors(n_cores=n_cores, n_rep=n_rep,
                                       n_sim=n_sim, ...))
  section("Step 2: Building results", build_all_results(n_rep=n_rep))
}

# =============================================================================
# 15.  Example calls
# =============================================================================

# ── Full pipeline ─────────────────────────────────────────────────────────
# source("ctbn_map_fast.R"); source("ctbn_map_simulation (3).R")
# source("ctbn_map_simulation_results (6).R")
# run_and_analyse(n_cores=20L, n_rep = 100, n_sim = 5000, k_fold = 5)

# ── Analyse completed replicates only ────────────────────────────────────
 source("ctbn_map_fast.R"); source("ctbn_map_simulation (3).R")
 source("ctbn_map_simulation_results (6).R")
 build_all_results()

# ── Single scenario × prior ──────────────────────────────────────────────
# run_simulation("A2", p_fit=2L, p_true=2L, prior="spike_slab", n_cores=20L)
# raw <- load_all_scenarios(scenarios="A2", priors="spike_slab")
# build_all_results(raw_list=raw, scenarios="A2", priors="spike_slab")
