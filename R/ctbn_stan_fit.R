# =============================================================================
# ctbn_results.R  (v4 — full revision with interaction effects)
# CTBN UK Biobank Multimorbidity Analysis — Empirical Results
#
# SECTIONS
# --------
# 5.2.1  Effect of Interaction Order on Predictive Performance
# 5.2.2  Sensitivity to the Order Penalty Parameter theta
# 5.2.3  Recommended Model
# 5.2.4  Inferred Multimorbidity Network
# 5.2.5  Conditional Risk: Main Effects (risk of 2nd LTC given 1 existing)
# 5.2.6  Interaction Effects: Synergistic Risk (risk of 3rd LTC given 2 existing)
# 5.2.7  Population-Averaged Risk Trajectories Stratified by Disease Burden
# 5.2.8  Survival and Cumulative Incidence Functions
# 5.2.9  Covariate Associations
#
# REVISION NOTES (v4)
# -------------------
# 1. Section 5.2.5 renamed to "Conditional Risk: Main Effects" and clarified
#    as answering Q1: risk of a second LTC given exactly one existing.
# 2. New Section 5.2.6: Interaction Effects — synergistic risk of a third LTC
#    given two existing conditions.  Implements:
#      a) Full interaction decomposition table: RR_j, RR_k, synergistic
#         multiplier exp(beta^m_{jk}), joint RR, absolute synergistic excess
#         Delta F_m(5), PIP_{jk}.
#      b) Forest plot of synergistic multipliers for all significant triplets.
#      c) Cumulative incidence curves with synergistic excess ribbon for
#         the five most striking triplets (joint F vs additive counterfactual).
# 3. All existing sections (5.2.7–5.2.9) renumbered accordingly.
# 4. compute_interaction_effects(): new function extracting pairwise
#    interaction betas and PIPs from the fitted spike-and-slab object.
# 5. compute_synergistic_excess(): computes Delta F_m(tau) for each
#    significant triplet (j, k -> m).
# =============================================================================

# --------------------------------------------------------------------------- #
# 0.  Packages
# --------------------------------------------------------------------------- #
library(data.table)
library(ggplot2)
library(scales)
library(patchwork)
library(xtable)
library(viridisLite)
library(reshape2)
library(igraph)
library(ggraph)
library(grid)
library(stringr)

dir.create("figures", showWarnings = FALSE)
dir.create("tables",  showWarnings = FALSE)

# ── Shared theme ──────────────────────────────────────────────────────────── #
theme_ctbn <- function(base_size = 11) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.grid.minor  = element_blank(),
      strip.background  = element_rect(fill = "grey92", colour = NA),
      strip.text        = element_text(face = "bold"),
      plot.title        = element_text(face = "bold", hjust = 0),
      plot.subtitle     = element_text(colour = "grey40", hjust = 0),
      legend.background = element_rect(colour = "grey80"),
      legend.key.size   = unit(0.45, "cm")
    )
}

# --------------------------------------------------------------------------- #
# 1.  Condition metadata and label system
# --------------------------------------------------------------------------- #

INFL_CONDS <- c(
  "Rhinitis and Sinusitis",
  "Dermatitis and eczema",
  "Chronic lower respiratory diseases"
)

SHORT_BASE <- c(
  "Hypertensive diseases"                                                           = "HTN",
  "Ischaemic heart diseases"                                                        = "IHD",
  "Rhinitis and Sinusitis"                                                          = "RS",
  "Dermatitis and eczema"                                                           = "DE",
  "Malignant neoplasms"                                                             = "MN",
  "Osteoarthritis"                                                                  = "OA",
  "Chronic lower respiratory diseases"                                              = "CLRD",
  "Diabetes"                                                                        = "DM",
  "Conductive and sensorineural hearing loss"                                       = "HL",
  "Diseases of veins, lymphatic vessels and lymph nodes, note eslewhere classified" = "DVL"
)

COND_SHORT <- c(
  "Hypertensive diseases"                                                           = "Hypertensive Dis.",
  "Ischaemic heart diseases"                                                        = "Ischaemic Heart Dis.",
  "Rhinitis and Sinusitis"                                                          = "Rhinitis & Sinusitis",
  "Dermatitis and eczema"                                                           = "Dermatitis & Eczema",
  "Malignant neoplasms"                                                             = "Malignant Neoplasms",
  "Osteoarthritis"                                                                  = "Osteoarthritis",
  "Chronic lower respiratory diseases"                                              = "Chr. Lower Resp. Dis.",
  "Diabetes"                                                                        = "Diabetes",
  "Conductive and sensorineural hearing loss"                                       = "Hearing Loss",
  "Diseases of veins, lymphatic vessels and lymph nodes, note eslewhere classified" = "Dis. Veins/Lymphatics"
)

make_full_label <- function(fullname) {
  abbr   <- SHORT_BASE[fullname]
  short  <- COND_SHORT[fullname]
  prefix <- ifelse(fullname %in% INFL_CONDS, "[Infl] ", "")
  paste0(prefix, short, " (", abbr, ")")
}

FULL_LABEL <- setNames(
  vapply(names(SHORT_BASE), make_full_label, character(1)),
  names(SHORT_BASE)
)

SHORT_LABELLED                <- SHORT_BASE
SHORT_LABELLED[INFL_CONDS]   <- paste0("[Infl] ", SHORT_BASE[INFL_CONDS])

NON_INFL_LABELS <- as.character(FULL_LABEL[!names(FULL_LABEL) %in% INFL_CONDS])
INFL_LABELS     <- as.character(FULL_LABEL[INFL_CONDS])
COND_ORDER      <- c(NON_INFL_LABELS, INFL_LABELS)

full_label  <- function(x) ifelse(x %in% names(FULL_LABEL),   FULL_LABEL[x],   x)
short_label <- function(x) ifelse(x %in% names(SHORT_LABELLED), SHORT_LABELLED[x], x)
cond_factor <- function(x) factor(x, levels = COND_ORDER)

INFL_COL <- "#e08214"

# ── RR fill scale ─────────────────────────────────────────────────────────── #
rr_fill_scale <- function(name = "RR") {
  scale_fill_gradient2(
    low = "#d73027", mid = "#f7f7f7", high = "#2166ac",
    midpoint = 1, trans = "log",
    breaks = c(0.25, 0.5, 1, 2, 4),
    labels = c("0.25","0.50","1.00","2.00","4.00"),
    na.value = "grey95", name = name
  )
}

# ── Inflammatory strip annotators ─────────────────────────────────────────── #
add_infl_strips <- function(gg, infl_levels = INFL_LABELS, all_levels = COND_ORDER,
                            bar_width = 0.28, gap = 0.62,
                            add_y = TRUE, add_x = TRUE) {
  pos <- which(all_levels %in% infl_levels)
  for (p in pos) {
    if (add_y)
      gg <- gg + annotate("rect",
                          xmin = -(gap), xmax = -(gap - bar_width),
                          ymin = p - 0.44, ymax = p + 0.44,
                          fill = INFL_COL, colour = NA, alpha = 0.90)
    if (add_x)
      gg <- gg + annotate("rect",
                          xmin = p - 0.44, xmax = p + 0.44,
                          ymin = -(gap), ymax = -(gap - bar_width),
                          fill = INFL_COL, colour = NA, alpha = 0.90)
  }
  gg + coord_cartesian(clip = "off")
}

add_infl_strip_x_only <- function(gg, infl_levels = INFL_LABELS,
                                  all_levels = COND_ORDER,
                                  n_rows, gap = 0.62, bar_w = 0.28) {
  pos <- which(all_levels %in% infl_levels)
  for (p in pos)
    gg <- gg + annotate("rect",
                        xmin = p - 0.44, xmax = p + 0.44,
                        ymin = n_rows + gap - bar_w, ymax = n_rows + gap,
                        fill = INFL_COL, colour = NA, alpha = 0.90)
  gg + coord_cartesian(clip = "off")
}

colour_infl_strips <- function(gg_object) {
  gt        <- ggplotGrob(gg_object)
  strip_idx <- grep("^strip", gt$layout$name)
  for (i in strip_idx) {
    tryCatch({
      lbl <- gt$grobs[[i]]$grobs[[1]]$children[[2]]$children[[1]]$label
      if (!is.null(lbl) && grepl("^\\[Infl\\]", lbl))
        gt$grobs[[i]]$grobs[[1]]$children[[1]]$gp$fill <- INFL_COL
    }, error = function(e) invisible(NULL))
  }
  gt
}

# ── Save helpers ──────────────────────────────────────────────────────────── #
save_fig <- function(p, name, w = 9, h = 6) {
  ggsave(file.path("figures", paste0(name, ".pdf")), p,
         width = w, height = h, units = "in", device = cairo_pdf)
  ggsave(file.path("figures", paste0(name, ".png")), p,
         width = w, height = h, units = "in", dpi = 300)
  message("Saved: figures/", name, ".{pdf,png}")
}

save_fig_gt <- function(gt, name, w = 14, h = 9) {
  cairo_pdf(file.path("figures", paste0(name, ".pdf")), width = w, height = h)
  grid::grid.draw(gt); dev.off()
  png(file.path("figures", paste0(name, ".png")),
      width = w, height = h, units = "in", res = 300)
  grid::grid.draw(gt); dev.off()
  message("Saved: figures/", name, ".{pdf,png}")
}

save_tex <- function(xt, name, ...) {
  path <- file.path("tables", paste0(name, ".tex"))
  print(xt, file = path, floating = FALSE, include.rownames = FALSE, ...)
  message("Saved: tables/", name, ".tex")
}

# =============================================================================
# CORE HELPERS: compute_F_m  and  build_interaction_cols
# =============================================================================

build_interaction_cols <- function(data, vars, k) {
  if (length(vars) < 2 || k < 2)
    return(list(data = data, cols = character(0), orders = integer(0)))
  k_eff <- min(k, length(vars))
  cols  <- character(0); orders <- integer(0)
  for (ord in 2:k_eff) {
    for (grp in combn(vars, ord, simplify = FALSE)) {
      cn         <- paste(grp, collapse = "_x_")
      data[[cn]] <- Reduce(`*`, lapply(grp, function(v) as.numeric(data[[v]])))
      cols       <- c(cols, cn); orders <- c(orders, ord - 1L)
    }
  }
  list(data = data, cols = cols, orders = orders)
}

compute_F_m <- function(fit, target, patient_dt, tau, zbar = NULL,
                        debug = FALSE, unit = 1) {
  ca            <- fit$call_args
  all_covs      <- c(ca$fixed_covs, ca$time_varying_covs)
  reserved_base <- c("eid","time_to_event","dt", all_covs, paste0(target,"_event"))
  orig_cols     <- names(patient_dt)
  all_conds     <- setdiff(orig_cols, reserved_base)
  influencers   <- setdiff(all_conds, target)
  
  stan_fit <- fit$stan_fits[[target]]
  if (is.null(stan_fit)) { warning("No Stan fit for: ", target); return(NA_real_) }
  
  post      <- as.data.frame(stan_fit)
  beta_mean <- colMeans(as.matrix(post[, grep("^beta\\[", names(post)), drop = FALSE]))
  gamma_mean <- colMeans(as.matrix(post[, grep("^gamma\\[", names(post)), drop = FALSE]))
  
  if (length(gamma_mean) != length(all_covs)) {
    warning(sprintf("[compute_F_m] gamma mismatch for %s", target)); return(NA_real_)
  }
  names(gamma_mean) <- all_covs
  
  nd     <- as.data.frame(copy(patient_dt))
  x_cols <- influencers
  if (!is.null(ca$max_order) && ca$max_order >= 2) {
    ir     <- build_interaction_cols(nd, influencers, ca$max_order)
    nd     <- ir$data; x_cols <- c(x_cols, ir$cols)
  }
  n_int <- nrow(nd)
  if (n_int == 0) return(0)
  
  if (length(beta_mean) < length(x_cols)) {
    warning(sprintf("[compute_F_m] beta mismatch for %s", target)); return(NA_real_)
  }
  beta_use <- beta_mean[seq_along(x_cols)]
  
  to_num <- function(v) {
    if (is.null(v)) return(rep(0, n_int))
    if (is.logical(v)) return(as.numeric(v))
    if (is.factor(v) || is.character(v)) return(as.numeric(as.factor(v)) - 1)
    as.numeric(v)
  }
  
  X_mat <- matrix(0, n_int, length(x_cols))
  for (j in seq_along(x_cols)) X_mat[, j] <- to_num(nd[[x_cols[j]]])
  
  Z_mat <- matrix(0, n_int, length(all_covs))
  for (k in seq_along(all_covs)) {
    cv <- all_covs[k]
    if (cv == "(Intercept)") { Z_mat[, k] <- 1
    } else if (!is.null(zbar) && cv %in% names(zbar)) { Z_mat[, k] <- as.numeric(zbar[cv])
    } else if (cv %in% names(nd))                     { Z_mat[, k] <- to_num(nd[[cv]])
    }
  }
  
  lp     <- as.numeric(X_mat %*% beta_use + Z_mat %*% gamma_mean)
  q_hat  <- exp(lp)
  t_start <- as.numeric(nd$time_to_event)
  t_end   <- t_start + as.numeric(nd$dt)
  in_win  <- t_start < tau
  t_star  <- pmax(0, pmin(t_end, tau) - t_start)
  pmin(pmax(1 - exp(-sum(q_hat[in_win] * t_star[in_win], na.rm = TRUE) * unit), 0), 1)
}



# =============================================================================
# 5.2.6  Interaction Effects: Synergistic Risk
#         Q2: What is the risk of a THIRD LTC given TWO existing conditions?
#
#  For each significant triplet (j, k -> m) where PIP_{jk}^m > 0.50:
#    RR_j            = exp(beta^m_{j})          [main effect of j alone]
#    RR_k            = exp(beta^m_{k})          [main effect of k alone]
#    syn_mult        = exp(beta^m_{jk})         [synergistic multiplier]
#    joint_RR        = RR_j * RR_k * syn_mult   [total RR given both j,k]
#    Delta F_m(tau)  = F(j,k) - F(j) - F(k) + F(empty)  [absolute synergy]
#    PIP_{jk}        = posterior inclusion prob for interaction term
# =============================================================================
cat("\n--- 5.2.6 ---\n")

# ── Helper: extract pairwise interaction betas and PIPs ───────────────────── #
# In the spike-and-slab Stan model, beta[1..J_main] = main effects,
# beta[J_main+1..J_total] = pairwise interaction terms in the order
# produced by build_interaction_cols(). The pip_beta vector follows the
# same ordering. We extract interaction betas and PIPs by matching
# the interaction column names produced by build_interaction_cols().

compute_interaction_effects <- function(fit, DT_wide, tau = 5,
                                        pip_thresh = 0.50) {
  ca        <- fit$call_args
  all_covs  <- c(ca$fixed_covs, ca$time_varying_covs)
  reserved  <- c("eid","time_to_event","dt", all_covs)
  all_conds <- setdiff(names(DT_wide), reserved)
  
  if (is.null(ca$max_order) || ca$max_order < 2) {
    message("Model fitted with max_order < 2; no pairwise interactions available.")
    return(data.table())
  }
  
  zbar <- sapply(all_covs, function(cv) {
    v <- DT_wide[[cv]]
    if (is.numeric(v)) mean(v, na.rm = TRUE)
    else as.numeric(names(which.max(table(v))))
  })
  names(zbar) <- all_covs
  
  results <- list()
  
  for (target in all_conds) {
    sfit <- fit$stan_fits[[target]]; if (is.null(sfit)) next
    influencers <- setdiff(all_conds, target)
    
    # Build column names in the exact order used during model fitting
    tmp_data <- as.data.frame(matrix(0, 1, length(influencers)))
    colnames(tmp_data) <- influencers
    ir       <- build_interaction_cols(tmp_data, influencers, ca$max_order)
    int_cols <- ir$cols   # e.g. "DM_x_HTN", "DM_x_IHD", ...
    
    if (length(int_cols) == 0) next
    
    post      <- as.data.frame(sfit)
    beta_cols <- grep("^beta\\[", names(post), value = TRUE)
    pip_cols  <- grep("^pip_beta\\[", names(post), value = TRUE)
    
    n_main <- length(influencers)
    n_int  <- length(int_cols)
    
    if (length(beta_cols) < n_main + n_int) next
    
    # Main effect posterior means (positions 1..n_main)
    beta_main_mean <- colMeans(post[, beta_cols[seq_len(n_main)], drop = FALSE])
    names(beta_main_mean) <- influencers
    
    # Interaction posterior means (positions n_main+1 .. n_main+n_int)
    beta_int_mean <- colMeans(post[, beta_cols[n_main + seq_len(n_int)], drop = FALSE])
    names(beta_int_mean) <- int_cols
    
    # Interaction PIPs (same positions in pip_beta)
    if (length(pip_cols) >= n_main + n_int) {
      pip_int_mean <- colMeans(post[, pip_cols[n_main + seq_len(n_int)], drop = FALSE])
      names(pip_int_mean) <- int_cols
    } else {
      pip_int_mean <- setNames(rep(NA_real_, n_int), int_cols)
    }
    
    # Iterate over all pairwise combinations of influencers
    for (pair in combn(influencers, 2, simplify = FALSE)) {
      j <- pair[1]; k <- pair[2]
      col_name <- paste(j, k, sep = "_x_")
      
      # Try both orderings (column name might be k_x_j)
      col_name_rev <- paste(k, j, sep = "_x_")
      actual_col <- if (col_name %in% int_cols) col_name else
        if (col_name_rev %in% int_cols) col_name_rev else NA
      
      if (is.na(actual_col)) next
      
      beta_jk  <- beta_int_mean[actual_col]
      pip_jk   <- pip_int_mean[actual_col]
      beta_j   <- beta_main_mean[j]
      beta_k   <- beta_main_mean[k]
      
      # Skip if interaction not significant
      if (!is.na(pip_thresh) && (is.na(pip_jk) || pip_jk < pip_thresh)) next
      
      # Compute cumulative incidences for the four profiles at covariates = zbar:
      #   F(empty):  all conditions = 0
      #   F(j):      X_j = 1, all others = 0
      #   F(k):      X_k = 1, all others = 0
      #   F(j,k):    X_j = 1, X_k = 1, all others = 0
      make_profile <- function(xj, xk) {
        ref <- as.data.table(as.list(rep(0, length(all_conds))))
        setnames(ref, all_conds)
        ref[, (j)      := as.numeric(xj)]
        ref[, (k)      := as.numeric(xk)]
        ref[, (target) := 0]
        ref[, eid           := 1L]
        ref[, time_to_event := 0]
        ref[, dt            := tau]
        for (cv in all_covs) ref[, (cv) := zbar[[cv]]]
        ref
      }
      
      F_empty <- tryCatch(compute_F_m(fit, target, make_profile(0, 0), tau),
                          error = function(e) NA_real_)
      F_j     <- tryCatch(compute_F_m(fit, target, make_profile(1, 0), tau),
                          error = function(e) NA_real_)
      F_k     <- tryCatch(compute_F_m(fit, target, make_profile(0, 1), tau),
                          error = function(e) NA_real_)
      F_jk    <- tryCatch(compute_F_m(fit, target, make_profile(1, 1), tau),
                          error = function(e) NA_real_)
      
      # Absolute synergistic excess: second-order finite difference
      # Delta F_m = F(j,k) - F(j) - F(k) + F(empty)
      delta_F <- F_jk - F_j - F_k + F_empty
      
      results[[length(results) + 1]] <- data.table(
        target_condition   = target,
        condition_j        = j,
        condition_k        = k,
        # Parameter decomposition
        beta_j             = round(beta_j,  4),
        beta_k             = round(beta_k,  4),
        beta_jk            = round(beta_jk, 4),
        RR_j               = round(exp(beta_j),  3),
        RR_k               = round(exp(beta_k),  3),
        syn_multiplier     = round(exp(beta_jk), 3),
        joint_RR           = round(exp(beta_j + beta_k + beta_jk), 3),
        # Absolute risk quantities
        F_baseline         = round(F_empty, 4),
        F_j_only           = round(F_j,     4),
        F_k_only           = round(F_k,     4),
        F_joint            = round(F_jk,    4),
        # Synergistic excess
        delta_F            = round(delta_F, 4),
        # PIP for interaction term
        pip_jk             = round(pip_jk,  3),
        tau_years          = tau
      )
    }
  }
  
  if (length(results) == 0) return(data.table())
  rbindlist(results)
}



# ============================================================
# Patient-grouped K-fold CV for CTBN prior comparison
# (Parallel edition — folds × models run concurrently)
# ============================================================
#
# PARALLELISATION STRATEGY
# ------------------------
# The two nested loops (folds × models) are the dominant cost.
# Both are embarrassingly parallel: each (fold, model) cell is
# fully independent of every other cell.
#
# Implementation uses {future} + {future.apply} so the same code
# runs on:
#   - multicore workers  (Linux / macOS): plan(multicore, workers = n)
#   - PSOCK clusters     (Windows / HPC): plan(multisession, workers = n)
#   - Sequential fallback (debugging)   : plan(sequential)
#
# Recommended usage:
#
#   library(future)
#   plan(multisession, workers = parallel::detectCores() - 1L)
#   cv_res <- ctbn_cv_compare(DT_wide, fit_fns, ...)
#   plan(sequential)   # restore after you're done
#
# WHAT IS PARALLELISED
# --------------------
#   outer level : fold_k  in 1..k_folds    (future_lapply)
#   inner level : model   in names(fit_fns) (future_lapply, nested)
#
# Because {future} uses a *fork* (multicore) or *socket* (multisession)
# model, all read-only objects (DT_wide, fit_fns, eval_times, …) are
# automatically available inside workers. No explicit export is needed.
#
# PROGRESS REPORTING
# ------------------
# Pass .progress = TRUE (requires {progressr}):
#
#   library(progressr)
#   handlers(global = TRUE)
#   with_progress(cv_res <- ctbn_cv_compare(..., .progress = TRUE))
#
# ============================================================

library(data.table)
library(timeROC)
library(survival)
library(future)
library(future.apply)


# ── Global helper: interaction column builder ─────────────────────────────────

build_interaction_cols <- function(data, vars, k) {
  if (length(vars) < 2 || k < 2)
    return(list(data = data, cols = character(0), orders = integer(0)))
  k_eff  <- min(k, length(vars))
  cols   <- character(0)
  orders <- integer(0)
  for (ord in 2:k_eff) {
    for (grp in combn(vars, ord, simplify = FALSE)) {
      cn         <- paste(grp, collapse = "_x_")
      data[[cn]] <- Reduce(`*`, lapply(grp, function(v) data[[v]]))
      cols       <- c(cols,   cn)
      orders     <- c(orders, ord - 1L)
    }
  }
  list(data = data, cols = cols, orders = orders)
}


# ── Fold assignment ────────────────────────────────────────────────────────────

make_patient_folds <- function(DT, k = 5, seed = 42) {
  eids    <- unique(DT$eid)
  N       <- length(eids)
  set.seed(seed)
  fold_id <- sample(rep(seq_len(k), length.out = N))
  data.table(eid = eids, fold = fold_id)
}


# ── Linear-predictor dispatcher ───────────────────────────────────────────────
#
# All get_lp methods return a named list:
#   $lp     numeric [n]  — log-rate (for ranking in timeROC)
#   $lambda numeric [n]  — rate = exp(lp) (for calibration metrics)

get_lp <- function(fit, newdata, target, ...) UseMethod("get_lp")

# ---- Unified ctbn_fit (class "ctbn_fit") ------------------------------------
get_lp.ctbn_fit <- function(fit, newdata, target, ...) {
  if (!is.null(fit$stan_fits[[target]]))
    return(get_lp.ctbn_fit_bayes(fit, newdata, target, ...))
  if (!is.null(fit$models[[target]]))
    return(get_lp.ctbn_fit_glm(fit, newdata, target, ...))
  warning(sprintf("get_lp: no model found for target '%s'", target))
  list(lp     = rep(NA_real_, nrow(newdata)),
       lambda = rep(NA_real_, nrow(newdata)))
}

# ---- GLM / MLE version -------------------------------------------------------
get_lp.ctbn_fit_glm <- function(fit, newdata, target, ...) {
  mod <- fit$models[[target]]
  if (is.null(mod)) {
    na <- rep(NA_real_, nrow(newdata))
    return(list(lp = na, lambda = na))
  }
  lp     <- as.numeric(predict(mod, newdata = as.data.frame(newdata), type = "link"))
  lambda <- exp(lp)
  list(lp = lp, lambda = lambda)
}

# ---- Stan / Bayesian version ------------------------------------------------
#
#   lp[i]     = X_i * beta_mean + Z_i * gamma_mean
#   lambda[i] = exp(lp[i])
#
# Posterior mean plug-in avoids Jensen bias and exp() overflow.

get_lp.ctbn_fit_bayes <- function(fit, newdata, target, ...) {
  
  newdata  <- copy(newdata)
  
  stan_fit <- fit$stan_fits[[target]]
  if (is.null(stan_fit)) {
    na <- rep(NA_real_, nrow(newdata))
    return(list(lp = na, lambda = na))
  }
  
  post    <- as.data.frame(stan_fit)
  b_samps <- as.matrix(post[, grep("^beta\\[",  names(post)), drop = FALSE])
  g_samps <- as.matrix(post[, grep("^gamma\\[", names(post)), drop = FALSE])
  
  beta_mean  <- colMeans(b_samps)
  gamma_mean <- colMeans(g_samps)
  
  ca       <- fit$call_args
  all_covs <- c(ca$fixed_covs, ca$time_varying_covs)
  reserved <- c("eid", "time_to_event", "dt", all_covs,
                paste0(target, "_event"))
  
  event_cols <- grep("_event$", names(newdata), value = TRUE)
  newdata[, (event_cols) := NULL]
  
  all_conditions_fit <- if (!is.null(ca$all_conditions))
    ca$all_conditions
  else {
    all_covs_reserved <- c("eid", "time_to_event", "dt", all_covs,
                           paste0(target, "_event"))
    setdiff(names(newdata), all_covs_reserved)
  }
  influencers <- setdiff(all_conditions_fit, target)
  
  x_cols   <- influencers
  x_orders <- rep(0L, length(influencers))
  nd       <- as.data.frame(newdata)
  
  if (!is.null(ca$max_order) && ca$max_order >= 2) {
    int_res  <- build_interaction_cols(nd, influencers, ca$max_order)
    nd       <- int_res$data
    x_cols   <- c(x_cols, int_res$cols)
    x_orders <- c(x_orders, int_res$orders)
  }
  
  to_num_col <- function(v) {
    if (is.null(v)) return(rep(0, nrow(nd)))
    if (is.factor(v) || is.character(v)) return(as.numeric(as.factor(v)) - 1)
    as.numeric(v)
  }
  
  X_new <- matrix(0, nrow = nrow(nd), ncol = length(x_cols))
  for (j in seq_along(x_cols)) X_new[, j] <- to_num_col(nd[[x_cols[j]]])
  
  if (length(all_covs) > 0) {
    Z_new <- matrix(0, nrow = nrow(nd), ncol = length(all_covs))
    for (k in seq_along(all_covs)) Z_new[, k] <- to_num_col(nd[[all_covs[k]]])
  } else {
    Z_new <- matrix(1, nrow = nrow(nd), ncol = 1)
  }
  
  if (length(beta_mean) != ncol(X_new) || length(gamma_mean) != ncol(Z_new)) {
    warning(sprintf(
      "get_lp [%s]: dimension mismatch — beta_mean %d vs X %d; gamma_mean %d vs Z %d",
      target, length(beta_mean), ncol(X_new), length(gamma_mean), ncol(Z_new)))
    na <- rep(NA_real_, nrow(nd))
    return(list(lp = na, lambda = na))
  }
  
  lp     <- as.numeric(X_new %*% beta_mean + Z_new %*% gamma_mean)
  lambda <- exp(lp)
  
  list(lp = lp, lambda = lambda)
}


# ── Per-(fold, model, target) evaluation cell ─────────────────────────────────
#
# Extracted from the inner loops so it can be called from any parallel
# backend (future, parallel::mclapply, etc.) without duplication.

.eval_cell <- function(fold_k, model_name, fit, DT_te,
                       all_conds, eval_times) {
  
  if (is.null(fit)) return(NULL)
  
  cell_results <- list()
  
  for (target in all_conds) {
    
    ec  <- paste0(target, "_event")
    atr <- DT_te[get(target) == 0 & dt > 0]
    if (nrow(atr) == 0) next
    
    pred <- tryCatch(
      get_lp(fit, atr, target),
      error = function(e) {
        warning(sprintf("get_lp failed [fold=%d, %s, %s]: %s",
                        fold_k, model_name, target, e$message))
        list(lp     = rep(NA_real_, nrow(atr)),
             lambda = rep(NA_real_, nrow(atr)))
      }
    )
    
    y      <- atr[[ec]]
    dt_    <- atr[["dt"]]
    lp     <- pred$lp
    lambda <- pred$lambda
    
    mu_hat <- lambda * pmax(dt_, 1e-10)
    
    poisson_ll <- mean(
      dpois(y, lambda = pmax(mu_hat, 1e-15), log = TRUE),
      na.rm = TRUE
    )
    
    p_hat <- 1 - exp(-pmax(mu_hat, 0))
    brier <- mean((y - p_hat)^2, na.rm = TRUE)
    
    pat_dt <- data.table(
      event   = y,
      t_event = atr$time_to_event,
      marker  = lp + log(pmax(dt_, 1e-10))
    )[is.finite(marker)]
    
    if (!is.null(eval_times) &&
        nrow(pat_dt) > 1 &&
        length(unique(pat_dt$event)) == 2) {
      
      td <- tryCatch(
        timeROC::timeROC(
          T      = pat_dt$t_event,
          delta  = pat_dt$event,
          marker = pat_dt$marker,
          cause  = 1,
          times  = eval_times,
          iid    = FALSE
        ),
        error = function(e) {
          warning(sprintf("timeROC failed [fold=%d, %s, %s]: %s",
                          fold_k, model_name, target, e$message))
          NULL
        }
      )
      
      for (t_idx in seq_along(eval_times)) {
        cell_results[[length(cell_results) + 1]] <- data.table(
          fold        = fold_k,
          model       = model_name,
          target      = target,
          poisson_ll  = poisson_ll,
          brier       = brier,
          eval_time   = eval_times[t_idx],
          tdauc       = if (!is.null(td)) td$AUC[t_idx] else NA_real_,
          n_test_rows = nrow(atr),
          n_test_pats = uniqueN(atr$eid)
        )
      }
      
    } else {
      cell_results[[length(cell_results) + 1]] <- data.table(
        fold        = fold_k,
        model       = model_name,
        target      = target,
        poisson_ll  = poisson_ll,
        brier       = brier,
        eval_time   = NA_real_,
        tdauc       = NA_real_,
        n_test_rows = nrow(atr),
        n_test_pats = uniqueN(atr$eid)
      )
    }
  }
  
  rbindlist(cell_results)
}


# ── Per-fold worker ────────────────────────────────────────────────────────────
#
# Called by future_lapply() — one future per fold.
# Fits all models, evaluates all targets, returns a list of data.tables.
#
# .progress_fn: a zero-argument function that signals one unit of work
#   (used by progressr when .progress = TRUE).

.fold_worker <- function(fold_k, DT_wide, fit_fns, all_conds,
                         eval_times, extra_args, .progress_fn = NULL) {
  
  DT_tr <- copy(DT_wide[fold != fold_k])[, fold := NULL]
  DT_te <- copy(DT_wide[fold == fold_k])[, fold := NULL]
  
  # Event indicators on test slice
  setorder(DT_te, eid, time_to_event)
  for (cond in all_conds) {
    ec <- paste0(cond, "_event")
    DT_te[, (ec) := as.numeric(
      get(cond) == 0 & shift(get(cond), type = "lead") == 1
    ), by = eid]
    DT_te[is.na(get(ec)), (ec) := 0L]
  }
  
  fold_results <- list()
  
  for (model_name in names(fit_fns)) {
    
    fit <- tryCatch(
      do.call(fit_fns[[model_name]],
              c(list(DT_tr, target_conditions = all_conds), extra_args)),
      error = function(e) {
        message(sprintf("  [ERROR] Fold %d, model '%s' FAILED: %s",
                        fold_k, model_name, conditionMessage(e)))
        NULL
      }
    )
    
    cell_dt <- .eval_cell(fold_k, model_name, fit, DT_te,
                          all_conds, eval_times)
    if (!is.null(cell_dt) && nrow(cell_dt) > 0)
      fold_results[[model_name]] <- cell_dt
    
    # Signal one completed (fold, model) unit to progressr
    if (!is.null(.progress_fn)) .progress_fn()
  }
  
  rbindlist(fold_results)
}


# ── Main CV loop (parallel) ────────────────────────────────────────────────────

ctbn_cv_compare <- function(DT_wide,
                            fit_fns,
                            fixed_covs        = character(0),
                            time_varying_covs = character(0),
                            k_folds           = 3,
                            targets           = NULL,
                            seed              = 42,
                            eval_times        = NULL,
                            .progress         = FALSE,
                            ...) {
  
  # ── Validation ───────────────────────────────────────────────────────────────
  dots <- list(...)
  if ("target_conditions" %in% names(dots))
    stop(paste0(
      "Do not pass 'target_conditions' via '...' or inside fit_fns lambdas.\n",
      "ctbn_cv_compare() injects target_conditions automatically from 'targets'.\n",
      "Remove target_conditions from your fit_fns definition."
    ))
  
  # ── Fold assignment ───────────────────────────────────────────────────────────
  fold_dt <- make_patient_folds(DT_wide, k = k_folds, seed = seed)
  DT_wide <- merge(copy(DT_wide), fold_dt, by = "eid")
  
  all_covs       <- c(fixed_covs, time_varying_covs)
  reserved       <- c("eid", "time_to_event", "dt", "fold", all_covs)
  all_conds_full <- setdiff(names(DT_wide), reserved)
  
  if (is.null(targets)) {
    all_conds <- all_conds_full
  } else {
    missing_targets <- setdiff(targets, all_conds_full)
    if (length(missing_targets))
      stop(sprintf("targets not found in DT_wide: %s",
                   paste(missing_targets, collapse = ", ")))
    all_conds <- targets
  }
  
  # ── Progress setup (optional — requires {progressr}) ────────────────────────
  # Total units = k_folds × n_models.  Each unit = one (fold, model) fit.
  p_fn <- NULL
  if (.progress) {
    if (!requireNamespace("progressr", quietly = TRUE))
      warning(".progress = TRUE requires the {progressr} package; ignoring.")
    else {
      p <- progressr::progressor(steps = k_folds * length(fit_fns))
      p_fn <- function() p()
    }
  }
  
  # ── Parallel fold dispatch ───────────────────────────────────────────────────
  # future_lapply() honours whatever plan() the caller has set:
  #   plan(sequential)               — single-threaded (default, safe)
  #   plan(multisession, workers=N)  — N background R sessions (Windows-safe)
  #   plan(multicore,    workers=N)  — forked workers (Linux/macOS only)
  #
  # future.seed ensures reproducible RNG across workers.
  
  closure_vars <- local({
    out <- list()
    for (fn in fit_fns) {
      ge  <- globalenv()
      nms <- ls(envir = ge, all.names = FALSE)
      for (nm in nms) {
        if (nm %in% names(out)) next
        val <- tryCatch(get(nm, envir = ge, inherits = FALSE),
                        error = function(x) NULL)
        if (!is.null(val) &&
            !inherits(val, "data.table") &&
            !inherits(val, "data.frame") &&
            !is.environment(val) &&
            (is.numeric(val) || is.character(val) ||
             is.logical(val) || is.list(val)))
          out[[nm]] <- val
      }
    }
    out
  })
  
  worker_globals <- c(
    list(
      #DT_wide        = DT_wide,
      fit_fns        = fit_fns,
      all_conds      = all_conds,
      all_conds_full = all_conds_full,
      eval_times     = eval_times
    ),
    closure_vars,
    list(
      .fold_worker           = .fold_worker,
      .eval_cell             = .eval_cell,
      get_lp                 = get_lp,
      get_lp.ctbn_fit        = get_lp.ctbn_fit,
      get_lp.ctbn_fit_bayes  = get_lp.ctbn_fit_bayes,
      get_lp.ctbn_fit_glm    = get_lp.ctbn_fit_glm,
      build_interaction_cols = build_interaction_cols,
      `%||%`                 = `%||%`,
      ctbn_fit               = ctbn_fit,
      .get_stan_model        = .get_stan_model,
      .ctbn_model_cache      = .ctbn_model_cache,
      .STAN_SPIKE_SLAB       = .STAN_SPIKE_SLAB,
      .STAN_STRUCTURED       = .STAN_STRUCTURED,
      .STAN_LASSO            = .STAN_LASSO,
      .STAN_HORSESHOE        = .STAN_HORSESHOE
    )
  )
  
  fold_list <- future.apply::future_lapply(
    X              = seq_len(k_folds),
    FUN            = function(fold_k) {
      .fold_worker(
        fold_k      = fold_k,
        DT_wide     = DT_wide,
        fit_fns     = fit_fns,
        all_conds   = all_conds,
        eval_times  = eval_times,
        extra_args  = dots,
        .progress_fn = p_fn
      )
    },
    future.globals  = worker_globals,
    future.seed    = seed,       # reproducible RNG in workers
    future.packages = c("rstan", "data.table", "timeROC", "survival")
  )
  
  rbindlist(fold_list)
}


# ── Summary helpers ────────────────────────────────────────────────────────────

summarise_cv <- function(cv_results) {
  
  scalar_dt <- unique(cv_results[, .(fold, model, target, poisson_ll, brier,
                                     n_test_rows, n_test_pats)])
  
  scalar_summary <- scalar_dt[, .(
    mean_poisson_ll = mean(poisson_ll, na.rm = TRUE),
    se_poisson_ll   = sd(poisson_ll,   na.rm = TRUE) / sqrt(.N),
    mean_brier      = mean(brier,      na.rm = TRUE),
    se_brier        = sd(brier,        na.rm = TRUE) / sqrt(.N),
    n_folds         = .N
  ), by = .(model, target)][order(target, -mean_poisson_ll)]
  
  auc_summary <- cv_results[!is.na(eval_time) & !is.na(tdauc), .(
    mean_tdauc = mean(tdauc, na.rm = TRUE),
    se_tdauc   = sd(tdauc,   na.rm = TRUE) / sqrt(.N),
    n_folds    = .N
  ), by = .(model, target, eval_time)][order(target, eval_time, model)]
  
  list(scalar = scalar_summary, tdauc = auc_summary)
}


# ── Plot time-dependent AUC ────────────────────────────────────────────────────

plot_tdauc <- function(cv_results, targets = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 required for plot_tdauc()")
  library(ggplot2)
  
  summ   <- summarise_cv(cv_results)
  auc_dt <- summ$tdauc
  if (!is.null(targets)) auc_dt <- auc_dt[target %in% targets]
  
  ggplot(auc_dt, aes(x = eval_time, y = mean_tdauc,
                     colour = model, fill = model, group = model)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_ribbon(aes(ymin = mean_tdauc - 1.96 * se_tdauc,
                    ymax = mean_tdauc + 1.96 * se_tdauc),
                alpha = 0.15, colour = NA) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
    facet_wrap(~target, scales = "free_y") +
    scale_y_continuous(limits = c(0.4, 1), labels = scales::percent_format(1)) +
    labs(
      title    = "Time-dependent AUC by model and target (K-fold CV)",
      subtitle = "Ribbon = \u00b11.96 SE across folds",
      x        = "Evaluation time",
      y        = "AUC\u207f(I/D)",
      colour   = "Prior / Model",
      fill     = "Prior / Model"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
}


# =============================================================================
# CTBN - UNIFIED BAYESIAN PRIOR ESTIMATION  (rstan)
# =============================================================================
#
# Supports four prior methods via the `prior` argument:
#
#   "spike_slab"  -- Marginalised normal mixture (spike-and-slab). DEFAULT.
#                    Beta[j] drawn from a two-component normal mixture;
#                    inclusion indicator delta_j analytically marginalised.
#                    PIP recovered per draw via Bayes' theorem.
#
#   "structured"  -- Order-penalised normal (continuous shrinkage).
#                    sigma2_p ~ Inv-Gamma(a0, b0) per order p.
#                    var(beta_j) = sigma2_p * exp(-theta * order_j).
#                    No selection — all coefficients retained; use posterior
#                    credible intervals or magnitude thresholds downstream.
#
#   "lasso"       -- Bayesian LASSO (double-exponential / Laplace).
#                    lambda2_p ~ Inv-Gamma(a0, b0) per order p.
#                    beta[j] ~ DE(0, lambda_p * exp(-theta * order_j)).
#                    The DE prior is equivalent to a scale mixture of normals
#                    with an Exp(1/(2*lambda^2)) mixing distribution, giving
#                    the exact Bayesian LASSO of Park & Casella (2008).
#                    No closed-form PIP; use credible intervals.
#
#   "horseshoe"   -- Regularised horseshoe (Piironen & Vehtari 2017).
#                    Global scale tau ~ Half-Cauchy(0, tau0).
#                    Local scales lambda[j] ~ Half-Cauchy(0, 1).
#                    Slab variance c2 ~ Inv-Gamma(slab_df/2, slab_df*slab_scale^2/2).
#                    beta[j] ~ N(0, (tau * lambda_tilde[j])^2) where
#                      lambda_tilde[j]^2 = c2*lambda[j]^2 / (c2 + tau^2*lambda[j]^2)
#                    tau further penalised by order: tau_j = tau * exp(-theta*order_j).
#                    This is the recommended "regularised" form — avoids the
#                    pathological tails of the plain horseshoe under weak data,
#                    while still giving near-zero shrinkage for large signals.
#
# PRIOR COMPARISON GUIDE
# ----------------------
#  Method        Sparsity  Selection  PIP  Order-penalised  Recommended when
#  ------------- --------- ---------  ---  ---------------  -----------------
#  spike_slab    Hard      Yes        Yes  Yes              Interpretability + selection
#  structured    Soft      No         No   Yes              Smooth shrinkage, no selection
#  lasso         Moderate  No*        No   Yes              Many small effects
#  horseshoe     Adaptive  No*        No   Yes              Few large + many near-zero
#
#  * Use 95% credible interval excluding 0 as a selection proxy for lasso/horseshoe.
#
# SHARED HYPERPARAMETERS
# ----------------------
#  pi0       : base inclusion probability  (spike_slab only)
#  theta     : order penalty; higher => stronger penalty on interactions
#  a0, b0    : Inv-Gamma shape/scale for slab / LASSO lambda variance
#  spike_var : spike variance              (spike_slab only)
#  tau0      : global scale prior          (horseshoe only)
#  slab_df, slab_scale : regularised-slab params (horseshoe only)
#
# =============================================================================

library(data.table)
library(rstan)

# ---------------------------------------------------------------------------
# Stan model strings
# ---------------------------------------------------------------------------

# ---- 1. Spike-and-slab (marginalised mixture) ----------------------------
.STAN_SPIKE_SLAB <- "
data {
  int<lower=1> N_rows;
  int<lower=1> Q;
  int<lower=1> P;
  array[N_rows] int<lower=0> N_obs;
  matrix[N_rows, Q] X;
  matrix[N_rows, P] Z;
  vector<lower=0>[N_rows] T;
  array[Q] int<lower=0> beta_order;
  real<lower=0, upper=1> pi0;
  real<lower=0>          theta;
  real<lower=0>          a0;
  real<lower=0>          b0;
  real<lower=0>          spike_var;
  int<lower=1>           max_order;
}
transformed data {
  vector[Q] pi_beta;
  for (j in 1:Q)
    pi_beta[j] = pi0 * exp(-theta * beta_order[j]);
  real pi_gamma = pi0;
}
parameters {
  vector[Q]                       beta;
  vector[P]                       gamma;
  vector<lower=0>[max_order + 1]  sigma2;
  real<lower=0>                   sigma2_gamma;
}
model {
  for (p in 0:max_order)
    sigma2[p + 1] ~ inv_gamma(a0, b0);
  sigma2_gamma ~ inv_gamma(a0, b0);

  for (j in 1:Q)
    target += log_mix(pi_beta[j],
                      normal_lpdf(beta[j]  | 0, sqrt(sigma2[beta_order[j] + 1])),
                      normal_lpdf(beta[j]  | 0, sqrt(spike_var)));

  for (k in 1:P)
    target += log_mix(pi_gamma,
                      normal_lpdf(gamma[k] | 0, sqrt(sigma2_gamma)),
                      normal_lpdf(gamma[k] | 0, sqrt(spike_var)));

  vector[N_rows] log_mu = X * beta + Z * gamma + log(T);
  target += poisson_log_lpmf(N_obs | log_mu);
}
generated quantities {
  vector[Q] pip_beta;
  vector[P] pip_gamma;
  for (j in 1:Q) {
    real ls = log(pi_beta[j])   + normal_lpdf(beta[j]  | 0, sqrt(sigma2[beta_order[j] + 1]));
    real lp = log1m(pi_beta[j]) + normal_lpdf(beta[j]  | 0, sqrt(spike_var));
    pip_beta[j] = exp(ls - log_sum_exp(ls, lp));
  }
  for (k in 1:P) {
    real ls = log(pi_gamma)   + normal_lpdf(gamma[k] | 0, sqrt(sigma2_gamma));
    real lp = log1m(pi_gamma) + normal_lpdf(gamma[k] | 0, sqrt(spike_var));
    pip_gamma[k] = exp(ls - log_sum_exp(ls, lp));
  }
}
"

# ---- 2. Structured (order-penalised normal) ------------------------------
.STAN_STRUCTURED <- "
data {
  int<lower=1> N_rows;
  int<lower=1> Q;
  int<lower=1> P;
  array[N_rows] int<lower=0> N_obs;
  matrix[N_rows, Q] X;
  matrix[N_rows, P] Z;
  vector<lower=0>[N_rows] T;
  array[Q] int<lower=0> beta_order;
  real<lower=0> theta;
  real<lower=0> a0;
  real<lower=0> b0;
  int<lower=1>  max_order;
}
parameters {
  vector[Q]                      beta;
  vector[P]                      gamma;
  vector<lower=0>[max_order + 1] sigma2;
  real<lower=0>                  sigma2_gamma;
}
model {
  // Order-specific slab variances; factorial penalty on the prior itself
  // captures the increasing implausibility of higher-order effects.
  for (p in 0:max_order)
    sigma2[p + 1] ~ inv_gamma(a0, b0);
  sigma2_gamma ~ inv_gamma(a0, b0);

  for (j in 1:Q) {
    // Effective variance shrinks exponentially with interaction order
    real eff_var = sigma2[beta_order[j] + 1] * exp(-theta * beta_order[j]);
    beta[j] ~ normal(0, sqrt(eff_var));
  }

  gamma ~ normal(0, sqrt(sigma2_gamma));

  vector[N_rows] log_mu = X * beta + Z * gamma + log(T);
  target += poisson_log_lpmf(N_obs | log_mu);
}
generated quantities {
  // No PIP for continuous shrinkage; return posterior SD per coefficient
  // for downstream credible-interval selection.
  vector[Q] post_sd_beta;
  for (j in 1:Q)
    post_sd_beta[j] = sqrt(sigma2[beta_order[j] + 1] * exp(-theta * beta_order[j]));
}
"

# ---- 3. Bayesian LASSO (Park & Casella 2008) -----------------------------
.STAN_LASSO <- "
data {
  int<lower=1> N_rows;
  int<lower=1> Q;
  int<lower=1> P;
  array[N_rows] int<lower=0> N_obs;
  matrix[N_rows, Q] X;
  matrix[N_rows, P] Z;
  vector<lower=0>[N_rows] T;
  array[Q] int<lower=0> beta_order;
  real<lower=0> theta;
  real<lower=0> a0;
  real<lower=0> b0;
  int<lower=1>  max_order;
}
parameters {
  vector[Q]                      beta;
  vector[P]                      gamma;
  // lambda2[p] is the per-order squared LASSO penalty parameter.
  // Placing Inv-Gamma on lambda^2 gives the Bayesian LASSO hierarchical model
  // of Park & Casella (2008) where the DE marginal arises from
  //   beta | tau2 ~ N(0, tau2),  tau2 ~ Exp(lambda^2/2).
  vector<lower=0>[max_order + 1] lambda2;
  real<lower=0>                  lambda2_gamma;
}
model {
  for (p in 0:max_order)
    lambda2[p + 1] ~ inv_gamma(a0, b0);
  lambda2_gamma ~ inv_gamma(a0, b0);

  for (j in 1:Q) {
    // Order-penalised Laplace scale: shrinks higher-order interactions more
    real scale_j = sqrt(lambda2[beta_order[j] + 1]) * exp(-theta * beta_order[j] / 2.0);
    beta[j] ~ double_exponential(0, scale_j);
  }

  gamma ~ double_exponential(0, sqrt(lambda2_gamma));

  vector[N_rows] log_mu = X * beta + Z * gamma + log(T);
  target += poisson_log_lpmf(N_obs | log_mu);
}
generated quantities {
  // Effective shrinkage factor kappa_j = 1 / (1 + 2*lambda2_eff_j)
  // (analogous to the ridge factor; near 1 => heavy shrinkage)
  vector[Q] kappa;
  for (j in 1:Q) {
    real lam2_eff = lambda2[beta_order[j] + 1] * exp(-theta * beta_order[j]);
    kappa[j] = 1.0 / (1.0 + 2.0 * lam2_eff);
  }
}
"

# ---- 4. Regularised Horseshoe (Piironen & Vehtari 2017) ------------------
.STAN_HORSESHOE <- "
data {
  int<lower=1> N_rows;
  int<lower=1> Q;
  int<lower=1> P;
  array[N_rows] int<lower=0> N_obs;
  matrix[N_rows, Q] X;
  matrix[N_rows, P] Z;
  vector<lower=0>[N_rows] T;
  array[Q] int<lower=0> beta_order;
  real<lower=0> theta;
  real<lower=0> tau0;        // global scale prior (set to p0/Q * 1/sqrt(N))
  real<lower=0> slab_df;     // regularising slab degrees of freedom (e.g. 4)
  real<lower=0> slab_scale;  // regularising slab scale (e.g. 2)
  int<lower=1>  max_order;
}
parameters {
  // Non-centred parameterisation for beta (better geometry in HMC)
  vector[Q]          z_beta;           // standard normal auxiliary
  vector[P]          gamma;
  vector<lower=0>[Q] lambda;           // local horseshoe scales
  real<lower=0>      tau;              // global horseshoe scale
  real<lower=0>      c2;              // slab variance (regularisation)
}
transformed parameters {
  vector[Q] beta;
  {
    // Regularised local scales: lambda_tilde^2 = c2*lambda^2 / (c2 + tau^2*lambda^2)
    // This caps the effective local scale at sqrt(c2), preventing unbounded tails.
    for (j in 1:Q) {
      real tau_j       = tau * exp(-theta * beta_order[j]);  // order penalty
      real lambda2_j   = square(lambda[j]);
      real lambda2_t   = c2 * lambda2_j / (c2 + square(tau_j) * lambda2_j);
      beta[j]          = z_beta[j] * sqrt(lambda2_t) * tau_j;
    }
  }
}
model {
  // Global scale: Half-Cauchy(0, tau0)
  tau ~ cauchy(0, tau0);

  // Slab variance: Inv-Gamma(slab_df/2, slab_df*slab_scale^2/2)
  // => c2 has a (scaled) Inv-chi-squared distribution, consistent with RHS paper
  c2 ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df * square(slab_scale));

  // Local scales: Half-Cauchy(0, 1)
  lambda ~ cauchy(0, 1);

  // Non-centred standard normals
  z_beta ~ std_normal();

  gamma ~ normal(0, 1);

  vector[N_rows] log_mu = X * beta + Z * gamma + log(T);
  target += poisson_log_lpmf(N_obs | log_mu);
}
generated quantities {
  // Shrinkage factor kappa_j = 1 / (1 + N_rows * tau_j^2 * lambda_tilde_j^2)
  // (Carvalho et al. 2010); kappa close to 1 => near-zero signal
  vector[Q] kappa;
  for (j in 1:Q) {
    real tau_j      = tau * exp(-theta * beta_order[j]);
    real lambda2_j  = square(lambda[j]);
    real lambda2_t  = c2 * lambda2_j / (c2 + square(tau_j) * lambda2_j);
    kappa[j] = 1.0 / (1.0 + N_rows * square(tau_j) * lambda2_t);
  }
}
"

# ---------------------------------------------------------------------------
# Stan model cache (one compiled model per prior type)
# ---------------------------------------------------------------------------

.ctbn_model_cache <- list()

.get_stan_model <- function(prior) {
  if (!is.null(.ctbn_model_cache[[prior]])) return(.ctbn_model_cache[[prior]])
  code <- switch(prior,
                 spike_slab = .STAN_SPIKE_SLAB,
                 structured = .STAN_STRUCTURED,
                 lasso      = .STAN_LASSO,
                 horseshoe  = .STAN_HORSESHOE,
                 stop("Unknown prior: '", prior, "'. Choose one of: spike_slab, structured, lasso, horseshoe.")
  )
  message(sprintf("Compiling Stan model for prior='%s' (once per session)...", prior))
  sm <- stan_model(model_code = code,
                   model_name = paste0("ctbn_", prior),
                   verbose    = FALSE)
  .ctbn_model_cache[[prior]] <<- sm
  sm
}

# =============================================================================
# ctbn_fit()  —  Unified entry point
# =============================================================================
#
# Arguments
# ---------
# DT_wide           : data.table, wide format; columns: eid, time_to_event, dt,
#                     [fixed/time-varying covariates], [condition columns]
# prior             : character, one of "spike_slab" | "structured" | "lasso" | "horseshoe"
# max_order         : maximum interaction order to include (0 = main effects only)
# fixed_covs        : character vector of fixed covariate column names
# time_varying_covs : character vector of time-varying covariate column names
# variable_select   : logical; if TRUE, zero out intensity for pairs below pip_threshold
# pip_threshold     : PIP threshold for selection (spike_slab only; ignored for others)
# alpha             : kept for backward compatibility (not used)
#
# Spike-slab hyperparameters (prior = "spike_slab"):
#   pi0, theta, a0, b0, spike_var
#
# Structured hyperparameters (prior = "structured"):
#   theta, a0, b0
#
# LASSO hyperparameters (prior = "lasso"):
#   theta, a0, b0
#
# Horseshoe hyperparameters (prior = "horseshoe"):
#   theta, tau0, slab_df, slab_scale
#
# Stan sampling:
#   chains, iter, warmup, seed, verbose
#
# =============================================================================

ctbn_fit <- function(DT_wide,
                     prior             = c("spike_slab", "structured", "lasso", "horseshoe"),
                     max_order         = 1,
                     fixed_covs        = character(0),
                     time_varying_covs = character(0),
                     target_conditions = NULL,   # NULL = fit all; else named subset
                     variable_select   = FALSE,
                     pip_threshold     = 0.5,
                     alpha             = 0.05,
                     # Spike-slab hyperparameters
                     pi0               = 0.5,
                     spike_var         = 0.01,
                     # Shared hyperparameters
                     theta             = 1.0,
                     a0                = 2.0,
                     b0                = 1.0,
                     # Horseshoe hyperparameters
                     tau0              = 1.0,
                     slab_df           = 4.0,
                     slab_scale        = 2.0,
                     # Stan
                     chains            = 4,
                     iter              = 2000,
                     warmup            = 1000,
                     seed              = 42L,
                     verbose           = TRUE) {
  
  prior <- match.arg(prior)
  
  # ---- Validate -----------------------------------------------------------
  stopifnot(is.data.table(DT_wide))
  stopifnot("eid"           %in% names(DT_wide))
  stopifnot("time_to_event" %in% names(DT_wide))
  stopifnot("dt"            %in% names(DT_wide))
  
  if (prior == "spike_slab") {
    if (pi0 <= 0 || pi0 >= 1) stop("pi0 must be strictly in (0, 1).")
    if (spike_var <= 0)        stop("spike_var must be positive.")
  }
  if (prior %in% c("spike_slab", "structured", "lasso")) {
    if (a0 <= 0 || b0 <= 0)   stop("a0 and b0 must be positive.")
  }
  
  all_covs       <- c(fixed_covs, time_varying_covs)
  reserved       <- c("eid", "time_to_event", "dt", all_covs)
  
  # all_conditions: every condition column present in DT_wide.
  # Used as the full influencer pool for every target model.
  all_conditions <- setdiff(names(DT_wide), reserved)
  n_cond         <- length(all_conditions)
  if (n_cond < 2) stop("Need at least 2 condition columns.")
  
  # target_loop: which conditions to model as outcomes.
  # target_conditions = NULL  => fit all (default behaviour).
  # target_conditions = c("Diabetes", "AF")  => fit only those two as
  #   outcomes while still using all_conditions as the influencer pool.
  if (is.null(target_conditions)) {
    target_loop <- all_conditions
  } else {
    bad <- setdiff(target_conditions, all_conditions)
    if (length(bad))
      stop(sprintf("target_conditions not found in DT_wide: %s",
                   paste(bad, collapse = ", ")))
    target_loop <- target_conditions
  }
  
  if (verbose) {
    message(sprintf("CTBN (%s prior): fitting %d / %d conditions as targets | max_order=%d",
                    prior, length(target_loop), n_cond, max_order))
    if (!is.null(target_conditions))
      message(sprintf("  Target subset : %s", paste(target_loop, collapse = ", ")))
    switch(prior,
           spike_slab = message(sprintf(
             "  Hyperparams: pi0=%.2f | theta=%.2f | spike_var=%.4f | Inv-Gamma(%.1f,%.1f)",
             pi0, theta, spike_var, a0, b0)),
           structured = message(sprintf(
             "  Hyperparams: theta=%.2f | Inv-Gamma(%.1f,%.1f)", theta, a0, b0)),
           lasso      = message(sprintf(
             "  Hyperparams: theta=%.2f | Inv-Gamma(%.1f,%.1f) [LASSO]", theta, a0, b0)),
           horseshoe  = message(sprintf(
             "  Hyperparams: theta=%.2f | tau0=%.2f | slab_df=%.1f | slab_scale=%.2f",
             theta, tau0, slab_df, slab_scale))
    )
  }
  
  # ---- Output containers --------------------------------------------------
  beta_matrix      <- matrix(0,        n_cond, n_cond,
                             dimnames = list(all_conditions, all_conditions))
  intensity_matrix <- matrix(0,        n_cond, n_cond,
                             dimnames = list(all_conditions, all_conditions))
  pip_matrix       <- matrix(NA_real_, n_cond, n_cond,
                             dimnames = list(all_conditions, all_conditions))
  # For non-spike-slab: store posterior shrinkage (kappa) instead of PIP
  kappa_matrix     <- matrix(NA_real_, n_cond, n_cond,
                             dimnames = list(all_conditions, all_conditions))
  pvalue_matrix    <- matrix(NA_real_, n_cond, n_cond,
                             dimnames = list(all_conditions, all_conditions))
  pip_cov_list     <- vector("list", n_cond)
  names(pip_cov_list)  <- all_conditions
  fitted_models    <- vector("list", n_cond)
  names(fitted_models) <- all_conditions
  stan_fits        <- vector("list", n_cond)
  names(stan_fits)     <- all_conditions
  ref_profiles     <- vector("list", n_cond)
  names(ref_profiles)  <- all_conditions
  
  # ---- Prepare data -------------------------------------------------------
  DT <- copy(DT_wide)
  setorder(DT, eid, time_to_event)
  
  # Compute event indicators for ALL conditions — needed because any condition
  # can appear as an at-risk filter (target == 0) even if it is not in target_loop.
  for (cond in all_conditions) {
    ec <- paste0(cond, "_event")
    DT[, (ec) := as.numeric(get(cond) == 0 &
                              shift(get(cond), type = "lead") == 1), by = eid]
    DT[is.na(get(ec)), (ec) := 0]
  }
  
  # ---- Helpers ------------------------------------------------------------
  build_interaction_cols <- function(data, vars, k) {
    if (length(vars) < 2 || k < 2)
      return(list(data = data, cols = character(0), orders = integer(0)))
    k_eff <- min(k, length(vars))
    cols <- character(0); orders <- integer(0)
    for (ord in 2:k_eff) {
      for (grp in combn(vars, ord, simplify = FALSE)) {
        cn         <- paste(grp, collapse = "_x_")
        data[[cn]] <- Reduce(`*`, lapply(grp, function(v) data[[v]]))
        cols   <- c(cols,   cn)
        orders <- c(orders, ord - 1L)
      }
    }
    list(data = data, cols = cols, orders = orders)
  }
  
  to_num <- function(col, n) {
    if (is.null(col)) return(rep(0, n))
    if (is.factor(col) || is.character(col)) return(as.numeric(as.factor(col)) - 1)
    as.numeric(col)
  }
  
  build_ref_profile <- function(atr) {
    ref <- data.frame(dt = 1)
    for (cond in all_conditions) ref[[cond]] <- 0
    for (cv in all_covs) {
      col <- atr[[cv]]
      if (is.null(col)) next
      if (is.factor(col) || is.character(col)) {
        tbl       <- sort(table(col), decreasing = TRUE)
        ref[[cv]] <- names(tbl)[1]
        if (is.factor(col)) ref[[cv]] <- factor(ref[[cv]], levels = levels(col))
      } else {
        ref[[cv]] <- mean(col, na.rm = TRUE)
      }
    }
    ref
  }
  
  sm <- .get_stan_model(prior)
  
  # ==========================================================================
  # Main loop: one Stan model per target condition
  # Iterates over target_loop only (subset when target_conditions specified).
  # influencers always drawn from the full all_conditions pool.
  # ==========================================================================
  for (target in target_loop) {
    
    ec          <- paste0(target, "_event")
    influencers <- setdiff(all_conditions, target)   # full pool as influencers
    
    atr <- DT[get(target) == 0 & dt > 0]
    if (nrow(atr) == 0) {
      warning(sprintf("No at-risk rows for '%s' — skipping.", target)); next
    }
    
    df <- as.data.frame(atr)
    n  <- nrow(df)
    
    # ---- Build X -----------------------------------------------------------
    int_res  <- build_interaction_cols(df, influencers, max_order)
    df       <- int_res$data
    x_cols   <- c(influencers, int_res$cols)
    x_orders <- c(rep(0L, length(influencers)), int_res$orders)
    
    X_mat <- matrix(0, nrow = n, ncol = length(x_cols))
    for (j in seq_along(x_cols)) X_mat[, j] <- to_num(df[[x_cols[j]]], n)
    colnames(X_mat) <- x_cols
    Q_stan <- ncol(X_mat)
    
    # ---- Build Z -----------------------------------------------------------
    if (length(all_covs) > 0) {
      z_cols <- all_covs
      Z_mat  <- matrix(0, nrow = n, ncol = length(z_cols))
      for (k in seq_along(z_cols)) Z_mat[, k] <- to_num(df[[z_cols[k]]], n)
      colnames(Z_mat) <- z_cols
    } else {
      z_cols <- "__intercept__"
      Z_mat  <- matrix(1, nrow = n, ncol = 1)
      colnames(Z_mat) <- z_cols
    }
    P_stan <- ncol(Z_mat)
    
    # ---- Assemble prior-specific Stan data ---------------------------------
    stan_data <- list(
      N_rows     = n,
      Q          = Q_stan,
      P          = P_stan,
      N_obs      = as.integer(df[[ec]]),
      X          = X_mat,
      Z          = Z_mat,
      T          = pmax(df[["dt"]], 1e-10),
      beta_order = as.array(x_orders),
      max_order  = as.integer(max_order),
      theta      = theta
    )
    
    # Add prior-specific fields
    if (prior == "spike_slab") {
      stan_data$pi0       <- pi0
      stan_data$a0        <- a0
      stan_data$b0        <- b0
      stan_data$spike_var <- spike_var
    } else if (prior %in% c("structured", "lasso")) {
      stan_data$a0 <- a0
      stan_data$b0 <- b0
    } else if (prior == "horseshoe") {
      stan_data$tau0       <- tau0
      stan_data$slab_df    <- slab_df
      stan_data$slab_scale <- slab_scale
    }
    
    # ---- Sample -----------------------------------------------------------
    fit <- tryCatch(
      sampling(sm,
               data    = stan_data,
               chains  = chains,
               iter    = iter,
               warmup  = warmup,
               seed    = seed,
               refresh = if (verbose) max(100L, iter %/% 10L) else 0L),
      error = function(e) {
        warning(sprintf("Stan FAILED for '%s': %s", target, e$message)); NULL
      }
    )
    if (is.null(fit)) next
    
    stan_fits[[target]]     <- fit
    fitted_models[[target]] <- fit
    
    post <- as.data.frame(fit)
    
    # ---- Extract posterior summaries -------------------------------------
    beta_samps  <- post[, grep("^beta\\[",  names(post), value = TRUE), drop = FALSE]
    gamma_samps <- post[, grep("^gamma\\[", names(post), value = TRUE), drop = FALSE]
    beta_means  <- colMeans(beta_samps)
    
    # Prior-specific: PIP (spike_slab) or kappa (lasso/horseshoe/structured)
    if (prior == "spike_slab") {
      pip_b_samps <- post[, grep("^pip_beta\\[",  names(post), value = TRUE), drop = FALSE]
      pip_g_samps <- post[, grep("^pip_gamma\\[", names(post), value = TRUE), drop = FALSE]
      pip_b_means <- colMeans(pip_b_samps)
      pip_g_means <- colMeans(pip_g_samps)
      
      for (idx in seq_along(influencers)) {
        inf <- influencers[idx]
        pip_matrix[inf, target] <- pip_b_means[idx]
      }
      pip_cov_list[[target]] <- setNames(pip_g_means, z_cols)
      
    } else if (prior %in% c("lasso", "horseshoe")) {
      kap_samps  <- post[, grep("^kappa\\[", names(post), value = TRUE), drop = FALSE]
      kap_means  <- colMeans(kap_samps)
      for (idx in seq_along(influencers)) {
        inf <- influencers[idx]
        kappa_matrix[inf, target] <- kap_means[idx]
      }
    }
    # structured: no selection statistic (use credible intervals externally)
    
    # ---- Store main-effect betas -----------------------------------------
    for (idx in seq_along(influencers)) {
      inf <- influencers[idx]
      beta_matrix[inf, target] <- beta_means[idx]
    }
    
    # ---- Reference profile -----------------------------------------------
    ref <- build_ref_profile(atr)
    ref_profiles[[target]] <- ref
    
    # ---- Intensities: E[exp(eta)] at ref with inf = 1 --------------------
    for (idx in seq_along(influencers)) {
      inf <- influencers[idx]
      
      # Selection gate: PIP for spike_slab; kappa<0.5 (signal) for others;
      # no gate for structured (all effects retained)
      skip <- FALSE
      if (variable_select) {
        if (prior == "spike_slab" &&
            !is.na(pip_matrix[inf, target]) &&
            pip_matrix[inf, target] < pip_threshold) skip <- TRUE
        if (prior %in% c("lasso", "horseshoe") &&
            !is.na(kappa_matrix[inf, target]) &&
            kappa_matrix[inf, target] > (1 - pip_threshold))  # kappa > 0.5 => heavy shrinkage
          skip <- TRUE
      }
      if (skip) next
      
      ref_inf        <- ref
      ref_inf[[inf]] <- 1
      df_new         <- as.data.frame(ref_inf)
      
      int_new <- build_interaction_cols(df_new, influencers, max_order)
      df_new  <- int_new$data
      
      x_new <- sapply(x_cols, function(cv) {
        v <- df_new[[cv]]; if (is.null(v)) 0 else as.numeric(v)[1]
      })
      z_new <- if ("__intercept__" %in% z_cols) {
        rep(1, P_stan)
      } else {
        sapply(z_cols, function(cv) {
          v <- df_new[[cv]]; if (is.null(v)) 0 else as.numeric(v)[1]
        })
      }
      
      eta_samps <- as.matrix(beta_samps) %*% x_new +
        as.matrix(gamma_samps) %*% z_new
      intensity_matrix[inf, target] <- mean(exp(eta_samps))
    }
    
    if (verbose) message(sprintf("  [OK] target: %s", target))
  }
  
  # ---- Cleanup event columns ----------------------------------------------
  ec_cols <- grep("_event$", names(DT), value = TRUE)
  if (length(ec_cols)) DT[, (ec_cols) := NULL]
  
  structure(
    list(
      prior            = prior,
      beta_matrix      = beta_matrix,
      intensity_matrix = intensity_matrix,
      pip_matrix       = pip_matrix,       # spike_slab only; else NA
      kappa_matrix     = kappa_matrix,     # lasso / horseshoe; else NA
      pip_cov_list     = pip_cov_list,     # spike_slab covariate PIP
      pvalue_matrix    = pvalue_matrix,    # NA (backward compat)
      models           = fitted_models,
      stan_fits        = stan_fits,
      ref_profiles     = ref_profiles,
      call_args        = list(
        prior             = prior,
        max_order         = max_order,
        fixed_covs        = fixed_covs,
        time_varying_covs = time_varying_covs,
        target_conditions = target_conditions,  # NULL = all
        all_conditions    = all_conditions,      # full influencer pool
        variable_select   = variable_select,
        pip_threshold     = pip_threshold,
        alpha             = alpha,
        pi0               = pi0,
        spike_var         = spike_var,
        theta             = theta,
        a0                = a0,
        b0                = b0,
        tau0              = tau0,
        slab_df           = slab_df,
        slab_scale        = slab_scale,
        chains            = chains,
        iter              = iter,
        warmup            = warmup
      )
    ),
    class = "ctbn_fit"
  )
}

# =============================================================================
# print.ctbn_fit
# =============================================================================

print.ctbn_fit <- function(x, digits = 4, ...) {
  ca <- x$call_args
  cat(sprintf("\n=== CTBN Fit (%s prior) ===\n", toupper(x$prior)))
  cat(sprintf("  Conditions      : %d\n", nrow(x$beta_matrix)))
  cat(sprintf("  Max order       : %d\n", ca$max_order))
  cat(sprintf("  theta (penalty) : %.2f\n", ca$theta))
  cat(sprintf("  Chains / iter   : %d / %d\n", ca$chains, ca$iter))
  
  switch(x$prior,
         spike_slab = {
           cat(sprintf("  pi0             : %.2f\n", ca$pi0))
           cat(sprintf("  spike_var       : %.4f\n", ca$spike_var))
           cat(sprintf("  Inv-Gamma(a0,b0): (%.1f, %.1f)\n", ca$a0, ca$b0))
         },
         structured = cat(sprintf(
           "  Inv-Gamma(a0,b0): (%.1f, %.1f)\n", ca$a0, ca$b0)),
         lasso = cat(sprintf(
           "  Inv-Gamma(a0,b0): (%.1f, %.1f) [LASSO lambda^2]\n", ca$a0, ca$b0)),
         horseshoe = cat(sprintf(
           "  tau0=%.2f | slab_df=%.1f | slab_scale=%.2f\n",
           ca$tau0, ca$slab_df, ca$slab_scale))
  )
  
  cat("\nPosterior mean beta (influencer -> target):\n")
  print(round(x$beta_matrix, digits))
  
  if (x$prior == "spike_slab") {
    cat("\nInfluencer PIP [E(delta_j|data)]:\n")
    print(round(x$pip_matrix, digits))
  } else if (x$prior %in% c("lasso", "horseshoe")) {
    cat("\nShrinkage kappa (0=signal, 1=shrunk to zero):\n")
    print(round(x$kappa_matrix, digits))
  }
  
  cat("\nIntensity matrix at reference (other conditions = 0):\n")
  print(round(x$intensity_matrix, digits))
  
  invisible(x)
}

# =============================================================================
# summary.ctbn_fit
# =============================================================================

summary.ctbn_fit <- function(object,
                             pip_threshold   = NULL,
                             show_covariates = FALSE,
                             ...) {
  
  if (is.null(pip_threshold)) {
    pip_threshold <- object$call_args$pip_threshold
    if (is.null(pip_threshold)) pip_threshold <- 0.5
  }
  
  pm    <- object$pip_matrix
  km    <- object$kappa_matrix
  bm    <- object$beta_matrix
  im    <- object$intensity_matrix
  pr    <- object$prior
  conds <- rownames(bm)
  
  # ---- Selection criterion -------------------------------------------------
  # spike_slab: PIP >= threshold
  # lasso/horseshoe: kappa <= (1 - threshold)  [i.e. signal dominates]
  # structured: all pairs shown (no selection statistic)
  is_active <- function(inf, tgt) {
    if (pr == "spike_slab") {
      !is.na(pm[inf, tgt]) && pm[inf, tgt] >= pip_threshold
    } else if (pr %in% c("lasso", "horseshoe")) {
      !is.na(km[inf, tgt]) && km[inf, tgt] <= (1 - pip_threshold)
    } else {
      TRUE  # structured: all
    }
  }
  
  rows <- list()
  for (inf in conds) {
    for (tgt in conds) {
      if (inf == tgt) next
      if (!is_active(inf, tgt)) next
      row <- data.frame(
        influencer = inf,
        target     = tgt,
        log_RR     = round(bm[inf, tgt], 4),
        RR         = round(exp(bm[inf, tgt]), 4),
        intensity  = round(im[inf, tgt], 6),
        stringsAsFactors = FALSE
      )
      if (pr == "spike_slab") {
        row$PIP   <- round(pm[inf, tgt], 3)
      } else if (pr %in% c("lasso", "horseshoe")) {
        row$kappa <- round(km[inf, tgt], 3)
      }
      rows[[length(rows) + 1]] <- row
    }
  }
  
  sel_label <- switch(pr,
                      spike_slab = sprintf("PIP >= %.2f", pip_threshold),
                      lasso      = sprintf("kappa <= %.2f", 1 - pip_threshold),
                      horseshoe  = sprintf("kappa <= %.2f", 1 - pip_threshold),
                      structured = "all pairs"
  )
  cat(sprintf("\n=== Active Influencer Pairs (%s, %s) ===\n", pr, sel_label))
  
  if (length(rows) == 0) {
    cat("  None\n")
    result_inf <- NULL
  } else {
    result_inf <- do.call(rbind, rows)
    sort_col   <- if (pr == "spike_slab") "PIP" else if (pr %in% c("lasso","horseshoe")) "kappa" else "log_RR"
    result_inf <- result_inf[order(result_inf[[sort_col]],
                                   decreasing = (pr == "spike_slab")), ]
    rownames(result_inf) <- NULL
    print(result_inf, ...)
  }
  
  # ---- Covariate PIPs (spike_slab only) -----------------------------------
  result_cov <- NULL
  if (show_covariates && pr == "spike_slab") {
    cov_rows <- list()
    for (tgt in names(object$pip_cov_list)) {
      pips <- object$pip_cov_list[[tgt]]
      if (is.null(pips)) next
      for (cv in names(pips)) {
        if (!is.na(pips[cv]) && pips[cv] >= pip_threshold)
          cov_rows[[length(cov_rows) + 1]] <- data.frame(
            covariate = cv, target = tgt,
            PIP = round(pips[cv], 3),
            stringsAsFactors = FALSE)
      }
    }
    if (length(cov_rows)) {
      result_cov <- do.call(rbind, cov_rows)
      result_cov <- result_cov[order(-result_cov$PIP), ]
      rownames(result_cov) <- NULL
      cat(sprintf("\n=== Active Covariates (PIP >= %.2f) ===\n", pip_threshold))
      print(result_cov, ...)
    }
  }
  
  invisible(list(influencers = result_inf, covariates = result_cov))
}


# =============================================================================
# Example calls (uncomment to run)
# =============================================================================

# --- Spike-and-slab (default, as before) ---
# result_ss <- ctbn_fit(DT_wide = DT_wide2[, -14:-15],
#                       prior      = "spike_slab",
#                       max_order  = 1,
#                       pi0        = 0.5,
#                       theta      = 1.0,
#                       a0         = 2.0,
#                       b0         = 1.0,
#                       spike_var  = 0.01,
#                       chains     = 1, iter = 1000, warmup = 500)

# --- Structured (order-penalised normal) ---
# result_str <- ctbn_fit(DT_wide = DT_wide2[, -14:-15],
#                        prior     = "structured",
#                        max_order = 1,
#                        theta     = 1.0,
#                        a0        = 2.0,
#                        b0        = 1.0,
#                        chains    = 1, iter = 1000, warmup = 500)

# --- Bayesian LASSO ---
# result_las <- ctbn_fit(DT_wide = DT_wide2[, -14:-15],
#                        prior     = "lasso",
#                        max_order = 1,
#                        theta     = 1.0,
#                        a0        = 2.0,
#                        b0        = 1.0,
#                        chains    = 1, iter = 1000, warmup = 500)

# --- Regularised Horseshoe ---
# result_hs <- ctbn_fit(DT_wide = DT_wide2[, -14:-15],
#                       prior      = "horseshoe",
#                       max_order  = 1,
#                       theta      = 1.0,
#                       tau0       = 1.0,
#                       slab_df    = 4.0,
#                       slab_scale = 2.0,
#                       chains     = 1, iter = 1000, warmup = 500)
