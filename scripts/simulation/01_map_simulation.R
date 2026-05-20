# =============================================================================
# ctbn_map_simulation.R
# CTBN MAP Simulation Study — Data Generation, Model Fitting, and Metrics
# =============================================================================
#
# REPLICATION OF ctbn_simulation_v3.R USING ctbn_map_fast
# --------------------------------------------------------
# This file is a complete self-contained replacement of ctbn_simulation_v3.R
# that uses ctbn_map_fast() (L-BFGS-B MAP + Laplace) instead of ctbn_fit()
# (Stan MCMC). Every function signature, output slot, and metric is preserved
# so that downstream result-analysis scripts work unchanged.
#
# WHAT IS IDENTICAL TO ctbn_simulation_v3.R
# ------------------------------------------
#   - True network parameters (TRUE_BETA0, TRUE_GRAPH, TRUE_BETA_MAIN,
#     TRUE_BETA_INT2, TRUE_BETA_INT3, TRUE_GAMMA)
#   - Data generation: simulate_patient(), generate_dataset(), prepare_wide()
#   - Smoking dynamics: sim_smoking_step(), smk_to_covs()
#   - All metric functions: compute_oracle_metrics(), compute_selection_metrics(),
#     compute_recovery_metrics(), compute_pred_metrics()
#   - Progress monitoring: SimProgress R6 class, make_progress_bar(),
#     format_duration(), rep_log(), inner_progress()
#   - Replicate caching and resumability via .rds files
#   - run_simulation() orchestrator with parallel workers
#   - run_simulation_with_live_progress() live-dashboard wrapper
#   - load_scenario()
#
# WHAT IS DIFFERENT
# -----------------
#   1. run_one_replicate() calls ctbn_map_fast() instead of ctbn_fit().
#      All four priors (spike_slab, structured, lasso, horseshoe) work
#      without Stan compilation.
#   2. compute_recovery_metrics() reads from map_fits instead of stan_fits.
#   3. compute_selection_metrics() uses pip_hat / kappa_hat from map_fits
#      instead of posterior MCMC draws.
#   4. get_lp dispatch uses get_lp.ctbn_map (MAP plug-in) — same interface,
#      no posterior samples needed.
#   5. Scenario defaults include all four priors:
#      PRIOR_FIT <- c("spike_slab", "structured", "lasso", "horseshoe")
#
# GLOBAL SETTINGS (edit as needed)
# ---------------------------------
SIM_SEED    <- 2026L
N_SIM       <- 5000L
N_REP       <- 3L
K_FOLD      <- 3L
EVAL_TIMES  <- 1L
T_HORIZON   <- 10
P_MAX_TRUE  <- 3L
P_FIT_GRID  <- c(1L, 2L, 3L)
THETA_GRID  <- c(0.5, 1.0, 2.0)
THETA_FIT   <- 1.0
PRIOR_FIT   <- c("spike_slab", "structured", "lasso", "horseshoe")
PIP_THRESH  <- 0.50

# =============================================================================

# --------------------------------------------------------------------------- #
# 0.  Packages
# --------------------------------------------------------------------------- #
suppressPackageStartupMessages({
  library(data.table)
  library(foreach)
  library(doParallel)
  library(pROC)
  library(timeROC)
  library(R6)
  library(jsonlite)
  library(survival)
  library(crayon)
})

# Source the MAP estimator (adjust path as needed)
# source("ctbn_map_fast.R")   # uncomment if not already sourced

# Directories
for (d in c("sim_results","sim_results/replicates","sim_results/figures",
            "sim_results/tables","sim_results/logs","sim_results/progress"))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# --------------------------------------------------------------------------- #
# 1.  True network parameters  (identical to ctbn_simulation_v3.R)
# --------------------------------------------------------------------------- #

CONDS <- c("HTN","IHD","MN","OA","DM","HL","DVL","RS","DE","CLRD")
M     <- length(CONDS)

TRUE_BETA0 <- setNames(log(c(
  HTN=0.150,IHD=0.080,MN=0.005,OA=0.060,DM=0.050,
  HL=0.040,DVL=0.040,RS=0.080,DE=0.070,CLRD=0.050)), CONDS)

TRUE_GRAPH <- list(
  HTN=c("DM","IHD","RS"), IHD=c("HTN","DM"), MN=c("HTN"),
  OA=c("DM","HTN"),       DM=c("HTN","IHD"),
  HL=c("OA","HTN"),       DVL=c("OA","DM"),
  RS=c("DE","CLRD"),      DE=c("RS","OA"),
  CLRD=c("RS","DM","HTN"))

TRUE_BETA_MAIN <- c(
  "HTN|DM"=0.75,"HTN|IHD"=0.35,"HTN|RS"=0.30,
  "IHD|HTN"=0.35,"IHD|DM"=0.55,
  "MN|HTN"=0.28,
  "OA|DM"=0.38,"OA|HTN"=0.32,
  "DM|HTN"=0.30,"DM|IHD"=0.35,
  "HL|OA"=0.32,"HL|HTN"=0.28,
  "DVL|OA"=0.25,"DVL|DM"=0.30,
  "RS|DE"=0.40,"RS|CLRD"=0.45,
  "DE|RS"=0.42,"DE|OA"=0.20,
  "CLRD|RS"=0.50,"CLRD|DM"=0.35,"CLRD|HTN"=0.28)

make_true_interactions_2way <- function(seed = SIM_SEED) {
  set.seed(seed)
  result <- list()
  for (m in names(TRUE_GRAPH)) {
    pa <- TRUE_GRAPH[[m]]
    if (length(pa) < 2) next
    for (pair in combn(pa, 2, simplify = FALSE)) {
      key <- paste(m, pair[1], pair[2], sep = "|")
      result[[key]] <- if (runif(1) < 0.40)
        (if (runif(1) < 0.70) rnorm(1,0,0.3) else rnorm(1,1.2,0.1)) else 0
    }
  }
  result
}

make_true_interactions_3way <- function(seed = SIM_SEED + 99L) {
  set.seed(seed)
  result <- list()
  for (m in names(TRUE_GRAPH)) {
    pa <- TRUE_GRAPH[[m]]
    if (length(pa) < 3) next
    for (triple in combn(pa, 3, simplify = FALSE)) {
      key <- paste(c(m, triple), collapse = "|")
      result[[key]] <- if (runif(1) < 0.30)
        (if (runif(1) < 0.70) rnorm(1,0,0.2) else rnorm(1,0.9,0.1)) else 0
    }
  }
  result
}

TRUE_BETA_INT2 <- make_true_interactions_2way()
TRUE_BETA_INT3 <- make_true_interactions_3way()

TRUE_GAMMA <- list(
  HTN =c(age=0.35,sex_male=0.20,smk_current=0.25,smk_former=0.20),
  IHD =c(age=0.30,sex_male=0.25,smk_current=0.40,smk_former=0.25),
  MN  =c(age=0.25,sex_male=0.25,smk_current=0.35,smk_former=0.20),
  OA  =c(age=0.45,sex_male=0.20,smk_current=0.25,smk_former=0.25),
  DM  =c(age=0.25,sex_male=0.25,smk_current=0.30,smk_former=0.22),
  HL  =c(age=0.40,sex_male=0.30,smk_current=0.20,smk_former=0.20),
  DVL =c(age=0.20,sex_male=0.20,smk_current=0.20,smk_former=0.28),
  RS  =c(age=0.20,sex_male=0.25,smk_current=0.25,smk_former=0.25),
  DE  =c(age=0.25,sex_male=-0.20,smk_current=0.20,smk_former=0.25),
  CLRD=c(age=0.30,sex_male=0.20,smk_current=0.60,smk_former=0.20))

# --------------------------------------------------------------------------- #
# 2.  Progress monitoring system  (identical to ctbn_simulation_v3.R)
# --------------------------------------------------------------------------- #

make_progress_bar <- function(current, total, width = 40) {
  frac   <- min(current / max(total, 1), 1)
  filled <- round(frac * width)
  bar    <- paste0(strrep("\u2588", filled), strrep("\u2591", width - filled))
  sprintf("[%s] %3.0f%%  %d/%d", bar, frac * 100, current, total)
}

format_duration <- function(secs) {
  secs <- as.numeric(secs)
  if (secs < 60)   return(sprintf("%.0fs", secs))
  if (secs < 3600) return(sprintf("%.0fm %.0fs", secs %/% 60, secs %% 60))
  sprintf("%.0fh %.0fm", secs %/% 3600, (secs %% 3600) %/% 60)
}

SimProgress <- R6::R6Class("SimProgress",
  private = list(
    state_file = NULL, lock_file = NULL, scenario = NULL,
    n_rep = NULL, k_fold = NULL, n_cond = NULL, start_time = NULL,

    acquire_lock = function(timeout = 30) {
      t0 <- proc.time()["elapsed"]
      while (file.exists(private$lock_file)) {
        if (proc.time()["elapsed"] - t0 > timeout) {
          warning("Progress lock timeout — proceeding anyway"); break
        }
        Sys.sleep(0.05)
      }
      writeLines(as.character(Sys.getpid()), private$lock_file)
    },

    release_lock = function() {
      if (file.exists(private$lock_file))
        try(file.remove(private$lock_file), silent = TRUE)
    },

    read_state = function() {
      if (!file.exists(private$state_file))
        return(list(completed=0L, failed=0L, rep_times=numeric(0), errors=list()))
      tryCatch(jsonlite::fromJSON(private$state_file),
               error = function(e)
                 list(completed=0L, failed=0L, rep_times=numeric(0), errors=list()))
    },

    write_state = function(state) {
      tryCatch(jsonlite::write_json(state, private$state_file, auto_unbox=TRUE),
               error = function(e) invisible(NULL))
    }
  ),
  public = list(
    initialize = function(scenario, n_rep, k_fold, n_cond, prior = "unknown") {
      private$scenario   <- scenario; private$n_rep  <- n_rep
      private$k_fold     <- k_fold;   private$n_cond <- n_cond
      private$start_time <- proc.time()["elapsed"]
      # Include prior in filenames so parallel cells don't collide
      key <- paste(scenario, prior, sep = "_")
      private$state_file <- file.path("sim_results/progress",
                                       sprintf("state_%s.json", key))
      private$lock_file  <- file.path("sim_results/progress",
                                       sprintf(".lock_%s", key))
      invisible(self)
    },

    increment = function(rep_idx, wall_time) {
      private$acquire_lock(); on.exit(private$release_lock())
      state <- private$read_state()
      state$completed <- state$completed + 1L
      state$rep_times <- c(tail(state$rep_times, 9), wall_time)
      private$write_state(state); invisible(self)
    },

    record_failure = function(rep_idx, msg) {
      private$acquire_lock(); on.exit(private$release_lock())
      state <- private$read_state()
      state$failed <- state$failed + 1L
      state$errors <- c(state$errors, list(list(
        rep=rep_idx, msg=as.character(msg), time=format(Sys.time()))))
      private$write_state(state); invisible(self)
    },

    print_dashboard = function() {
      state   <- private$read_state()
      done    <- state$completed; failed <- state$failed
      n       <- private$n_rep
      elapsed <- proc.time()["elapsed"] - private$start_time
      bar     <- make_progress_bar(done, n)
      eta_str <- ""
      if (done > 0 && length(state$rep_times) > 0) {
        mean_t  <- mean(state$rep_times)
        eta_str <- sprintf("  ETA: %s", format_duration(mean_t * (n - done)))
      }
      cat("\r\033[K")
      cat(sprintf("%s [Scen %s]  Elapsed: %s%s  Failed: %d",
                  bar, private$scenario, format_duration(elapsed),
                  eta_str, failed))
      utils::flush.console()
      invisible(state)
    },

    final_report = function() {
      state   <- private$read_state()
      elapsed <- proc.time()["elapsed"] - private$start_time
      cat(sprintf("\n\n%s Scenario %s complete\n",
                  format(Sys.time(), "%Y-%m-%d %H:%M:%S"), private$scenario))
      cat(sprintf("  Replicates completed : %d / %d\n",
                  state$completed, private$n_rep))
      cat(sprintf("  Replicates failed    : %d\n", state$failed))
      cat(sprintf("  Total elapsed        : %s\n", format_duration(elapsed)))
      if (length(state$rep_times) > 0)
        cat(sprintf("  Mean rep time        : %s\n",
                    format_duration(mean(state$rep_times))))
      if (state$failed > 0) {
        cat("  Errors:\n")
        for (e in state$errors)
          cat(sprintf("    Rep %03d [%s]: %s\n", e$rep, e$time,
                      substr(e$msg, 1, 120)))
      }
      invisible(state)
    },

    rep_done = function(rep_idx, scenario, prior = "unknown") {
      file.exists(file.path("sim_results/replicates",
                             sprintf("rep_%s_%s_%03d.rds", scenario, prior, rep_idx)))
    },

    n_remaining = function(scenario, prior = "unknown") {
      paths <- file.path("sim_results/replicates",
        sprintf("rep_%s_%s_%03d.rds", scenario, prior, seq_len(private$n_rep)))
      sum(!file.exists(paths))
    }
  )
)

rep_log <- function(scenario, rep_idx, msg, level = "INFO") {
  log_path <- file.path("sim_results/logs",
                         sprintf("scenario_%s.log", scenario))
  line <- sprintf("[%s] [%s] [Rep %03d] %s\n",
                  format(Sys.time(), "%H:%M:%S"), level, rep_idx, msg)
  tryCatch(cat(line, file = log_path, append = TRUE),
           error = function(e) invisible(NULL))
}

inner_progress <- function(scenario, rep_idx, phase, item, total_items) {
  rep_log(scenario, rep_idx,
          sprintf("%s | %s (%d/%d)", phase, item,
                  total_items[[item]], total_items[["total"]]),
          level = "DETAIL")
}

# --------------------------------------------------------------------------- #
# 3.  Data generation  (identical to ctbn_simulation_v3.R)
# --------------------------------------------------------------------------- #

compute_q_true <- function(m, X_state, covs, include_3way = TRUE) {
  lp <- TRUE_BETA0[m]
  pa <- TRUE_GRAPH[[m]]
  for (j in pa) {
    key <- paste(m, j, sep = "|")
    if (key %in% names(TRUE_BETA_MAIN))
      lp <- lp + TRUE_BETA_MAIN[key] * X_state[j]
  }
  if (length(pa) >= 2)
    for (pair in combn(pa, 2, simplify = FALSE)) {
      key <- paste(m, pair[1], pair[2], sep = "|")
      if (key %in% names(TRUE_BETA_INT2) && TRUE_BETA_INT2[[key]] != 0)
        lp <- lp + TRUE_BETA_INT2[[key]] * X_state[pair[1]] * X_state[pair[2]]
    }
  if (include_3way && length(pa) >= 3)
    for (triple in combn(pa, 3, simplify = FALSE)) {
      key <- paste(c(m, triple), collapse = "|")
      if (key %in% names(TRUE_BETA_INT3) && TRUE_BETA_INT3[[key]] != 0)
        lp <- lp + TRUE_BETA_INT3[[key]] *
          X_state[triple[1]] * X_state[triple[2]] * X_state[triple[3]]
    }
  gamma_m <- TRUE_GAMMA[[m]]
  lp <- lp + sum(gamma_m * covs[names(gamma_m)])
  exp(lp)
}

sim_smoking_step <- function(smk, dt, q01 = 0.04, q12 = 0.06) {
  if (smk == 2L) return(2L)
  if (smk == 0L) return(if (runif(1) < 1-exp(-q01*dt)) 1L else 0L)
  if (runif(1) < 1-exp(-q12*dt)) 2L else 1L
}

smk_to_covs <- function(smk)
  c(smk_current = as.numeric(smk == 1L), smk_former = as.numeric(smk == 2L))

simulate_patient <- function(eid, t_horizon = T_HORIZON, dt_step = 5,
                              include_3way = TRUE) {
  init_prev <- c(HTN=0.16,IHD=0.08,MN=0.06,OA=0.06,DM=0.05,
                 HL=0.04,DVL=0.04,RS=0.08,DE=0.07,CLRD=0.05)
  X_state   <- setNames(rbinom(M, 1, init_prev), CONDS)
  age0      <- log(rnorm(1, 57, 9.39))
  sex_male  <- rbinom(1, 1, 0.56)
  smk       <- sample(0:2, 1, prob=c(0.47, 0.11, 0.42))
  records   <- list(); t_now <- 0

  while (t_now < t_horizon) {
    t_next   <- min(t_now + dt_step, t_horizon)
    dt       <- t_next - t_now
    smk_covs <- smk_to_covs(smk)
    covs     <- c(age=age0+t_now, sex_male=sex_male, smk_covs)
    at_risk  <- CONDS[X_state == 0]

    if (length(at_risk) == 0) {
      smk <- sim_smoking_step(smk, dt); t_now <- t_next; next
    }

    q_vec  <- setNames(
      sapply(at_risk, compute_q_true, X_state=X_state, covs=covs,
             include_3way=include_3way), at_risk)
    Lambda <- sum(q_vec)
    t_cand <- rexp(1, max(Lambda, 1e-10))

    if (t_cand < dt) {
      m_ev <- sample(at_risk, 1, prob = q_vec / Lambda)
      t_ev <- t_now + t_cand
      records[[length(records)+1]] <- data.table(
        eid=eid, time_to_event=t_now, dt=t_ev-t_now,
        matrix(t(X_state), nrow=1, dimnames=list(NULL, CONDS)),
        age=age0+t_now, sex_male=sex_male,
        smk_current=smk_covs["smk_current"],
        smk_former =smk_covs["smk_former"],
        event_cond=m_ev)
      X_state[m_ev] <- 1L; t_now <- t_ev
      if (all(X_state==1)) break
    } else {
      records[[length(records)+1]] <- data.table(
        eid=eid, time_to_event=t_now, dt=dt,
        matrix(t(X_state), nrow=1, dimnames=list(NULL, CONDS)),
        age=age0+t_now, sex_male=sex_male,
        smk_current=smk_covs["smk_current"],
        smk_former =smk_covs["smk_former"],
        event_cond=NA_character_)
      smk <- sim_smoking_step(smk, dt); t_now <- t_next
    }
  }
  if (length(records)==0) return(NULL)
  rbindlist(records)
}

generate_dataset <- function(n = N_SIM, seed = NULL, include_3way = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  rbindlist(lapply(seq_len(n), simulate_patient, include_3way=include_3way))
}

prepare_wide <- function(DT) {
  DT[, event_cond := NULL]
  DT[, (CONDS) := lapply(.SD, as.integer), .SDcols=CONDS]
  DT
}

# --------------------------------------------------------------------------- #
# 4.  Oracle metrics  (identical to ctbn_simulation_v3.R)
# --------------------------------------------------------------------------- #

compute_oracle_metrics <- function(DT_test, include_3way = TRUE) {
  results <- list()
  for (m in CONDS) {
    ec  <- paste0(m,"_event")
    atr <- DT_test[get(m)==0 & dt>0]
    if (nrow(atr)==0 || !ec %in% names(atr)) next
    y <- atr[[ec]]; dt_ <- atr$dt
    q_true <- vapply(seq_len(nrow(atr)), function(ri) {
      row  <- atr[ri]
      X_st <- unlist(row[, ..CONDS])
      covs <- c(age=row$age, sex_male=row$sex_male,
                smk_current=row$smk_current, smk_former=row$smk_former)
      compute_q_true(m, X_st, covs, include_3way)
    }, numeric(1))
    mu  <- q_true * pmax(dt_, 1e-10)
    p   <- 1 - exp(-mu)
    results[[m]] <- data.table(
      condition    = m,
      oracle_pll   = mean(dpois(y, pmax(mu, 1e-15), log=TRUE), na.rm=TRUE),
      oracle_brier = mean((y-p)^2, na.rm=TRUE))
  }
  rbindlist(results)
}

# --------------------------------------------------------------------------- #
# 5.  Selection metrics  (adapted for MAP — reads pip_hat / kappa_hat)
# --------------------------------------------------------------------------- #

permn_3 <- function(x)
  matrix(c(x[1],x[2],x[3], x[1],x[3],x[2],
           x[2],x[1],x[3], x[2],x[3],x[1],
           x[3],x[1],x[2], x[3],x[2],x[1]),
         ncol=3, byrow=TRUE)

#' Compute selection metrics from a ctbn_map fit
#'
#' Reads pip_hat (spike_slab/structured/lasso) or kappa_hat (horseshoe)
#' directly from the per-node MAP fits — no MCMC posterior samples needed.
compute_selection_metrics <- function(fit, pip_thresh = PIP_THRESH, p_true = 2L) {
  true_main_nz <- names(TRUE_BETA_MAIN)[TRUE_BETA_MAIN != 0]
  true_2way_nz <- names(TRUE_BETA_INT2)[unlist(TRUE_BETA_INT2) != 0]
  true_3way_nz <- if (p_true >= 3L)
    names(TRUE_BETA_INT3)[unlist(TRUE_BETA_INT3) != 0] else character(0)

  is_map <- inherits(fit, "ctbn_map")
  prior  <- fit$prior
  p_fit  <- fit$call_args$max_order

  # Score function: returns the selection score for position col_pos in node m
  mk_score <- function(mfit, col_pos) {
    if (!is_map || is.null(mfit)) return(NA_real_)
    if (prior == "horseshoe") {
      kh <- mfit$kappa_hat
      if (!is.null(kh) && col_pos <= length(kh)) return(1.0 - kh[col_pos])
    } else {
      ph <- mfit$pip_hat
      if (!is.null(ph) && col_pos <= length(ph)) return(ph[col_pos])
    }
    NA_real_
  }

  is_selected <- function(score) {
    if (is.na(score)) return(FALSE)
    if (prior == "horseshoe") score >= pip_thresh   # pseudo-PIP = 1 - kappa
    else                      score >= pip_thresh
  }

  r_main <- r_2way <- r_3way <- list()

  for (m in CONDS) {
    mfit <- fit$map_fits[[m]]; if (is.null(mfit)) next
    infl  <- setdiff(CONDS, m)
    n_m   <- length(infl)
    x_cols   <- mfit$x_cols
    x_orders <- mfit$x_orders

    # Main effects
    for (idx in seq_along(infl)) {
      j   <- infl[idx]; key <- paste(m, j, sep="|")
      sc  <- mk_score(mfit, idx)
      r_main[[length(r_main)+1]] <- data.table(
        target=m, influencer=j,
        true_nz  = key %in% true_main_nz,
        selected = is_selected(sc),
        score    = sc)
    }

    # 2-way interactions
    if (p_fit >= 2L) {
      pairs    <- combn(infl, 2, simplify=FALSE)
      int_pos  <- which(x_orders == 1L)
      for (pi in seq_along(pairs)) {
        j <- pairs[[pi]][1]; k <- pairs[[pi]][2]
        cn1 <- paste(j,k,sep="_x_"); cn2 <- paste(k,j,sep="_x_")
        ac  <- if (cn1 %in% x_cols) cn1 else if (cn2 %in% x_cols) cn2 else NA
        if (is.na(ac)) next
        pos <- which(x_cols == ac)
        sc  <- mk_score(mfit, pos)
        k1  <- paste(m,j,k,sep="|"); k2 <- paste(m,k,j,sep="|")
        r_2way[[length(r_2way)+1]] <- data.table(
          target=m, cond_j=j, cond_k=k,
          true_nz  = k1 %in% true_2way_nz || k2 %in% true_2way_nz,
          selected = is_selected(sc),
          score    = sc)
      }
    }

    # 3-way interactions
    if (p_fit >= 3L && p_true >= 3L) {
      triples <- combn(infl, 3, simplify=FALSE)
      for (ti in seq_along(triples)) {
        j<-triples[[ti]][1]; k<-triples[[ti]][2]; l<-triples[[ti]][3]
        perms <- apply(permn_3(c(j,k,l)), 1,
                       function(x) paste(x, collapse="_x_"))
        ac    <- perms[perms %in% x_cols]
        if (!length(ac)) next
        pos   <- which(x_cols == ac[1])
        sc    <- mk_score(mfit, pos)
        keys  <- apply(permn_3(c(j,k,l)), 1,
                       function(x) paste(c(m,x), collapse="|"))
        r_3way[[length(r_3way)+1]] <- data.table(
          target=m, cond_j=j, cond_k=k, cond_l=l,
          true_nz  = any(keys %in% true_3way_nz),
          selected = is_selected(sc),
          score    = sc)
      }
    }
  }

  calc_sel <- function(dt) {
    if (is.null(dt) || nrow(dt)==0)
      return(list(tpr=NA_real_, fpr=NA_real_, auc=NA_real_))
    nz  <- dt[true_nz==TRUE]; z <- dt[true_nz==FALSE]
    tpr <- if (nrow(nz)>0) mean(nz$selected, na.rm=TRUE) else NA_real_
    fpr <- if (nrow(z)>0)  mean(z$selected,  na.rm=TRUE) else NA_real_
    av  <- tryCatch(
      if (!all(is.na(dt$score)) && length(unique(dt$true_nz))==2)
        as.numeric(auc(roc(as.numeric(dt$true_nz), dt$score,
                           quiet=TRUE, direction="<")))
      else NA_real_,
      error=function(e) NA_real_)
    list(tpr=tpr, fpr=fpr, auc=av)
  }

  ms <- calc_sel(rbindlist(r_main, fill=TRUE))
  ts <- calc_sel(rbindlist(r_2way, fill=TRUE))
  hs <- calc_sel(rbindlist(r_3way, fill=TRUE))

  data.table(
    level=c("main_effect","interaction_2way","interaction_3way"),
    tpr=c(ms$tpr,ts$tpr,hs$tpr),
    fpr=c(ms$fpr,ts$fpr,hs$fpr),
    selection_auc=c(ms$auc,ts$auc,hs$auc))
}

# --------------------------------------------------------------------------- #
# 6.  Recovery metrics  (adapted for MAP — reads beta_hat from map_fits)
# --------------------------------------------------------------------------- #

#' Compute parameter recovery metrics from a ctbn_map fit
#'
#' Reads MAP point estimates and Laplace SEs directly from map_fits.
#' Returns the same data.table structure as the Stan version in
#' ctbn_simulation_v3.R so downstream analysis scripts work unchanged.
compute_recovery_metrics <- function(fit, p_true = 2L) {
  p_fit   <- fit$call_args$max_order
  results <- list()

  mk <- function(pt, pn, tv, e, v, lo, hi)
    data.table(condition=NA_character_, parameter=pn, param_type=pt,
               true_val=tv, est_mean=e, est_var=v,
               ci_lo=lo, ci_hi=hi,
               bias=e-tv, sq_err=(e-tv)^2,
               covered=as.numeric(tv>=lo & tv<=hi))

  for (m in CONDS) {
    mfit <- fit$map_fits[[m]]; if (is.null(mfit)) next
    bh   <- mfit$beta_hat; seh <- mfit$se_hat
    gh   <- mfit$gamma_hat
    x_cols   <- mfit$x_cols
    x_orders <- mfit$x_orders
    n_beta   <- mfit$n_beta; n_gamma <- mfit$n_gamma

    # Use the dedicated se_gamma slot from ctbn_map_fast (v2.0+).
    # Fall back to extracting from the full se_hat vector for backward compat.
    gse <- if (!is.null(mfit$se_gamma))
      mfit$se_gamma
    else if (!is.null(mfit$se_hat) && length(mfit$se_hat) > n_beta)
      mfit$se_hat[n_beta + seq_len(n_gamma)]
    else
      rep(NA_real_, n_gamma)

    # Covariate PIPs (pip_gamma_hat slot from ctbn_map_fast v2.0+)
    pg_vec <- mfit$pip_gamma_hat %||% rep(NA_real_, n_gamma)
    kg_vec <- mfit$kappa_gamma_hat %||% rep(NA_real_, n_gamma)

    # ── Baseline (intercept-like first covariate coefficient)
    r_base <- mk("baseline","beta0", TRUE_BETA0[m],
                 gh[1], pmax(gse[1],0,na.rm=TRUE)^2,
                 gh[1] - 1.96*pmax(gse[1],0,na.rm=TRUE),
                 gh[1] + 1.96*pmax(gse[1],0,na.rm=TRUE))
    r_base$condition <- m
    r_base$pip <- pg_vec[1]
    results[[length(results)+1]] <- r_base

    infl <- setdiff(CONDS, m)

    # ── Main effects
    for (idx in seq_along(infl)) {
      j  <- infl[idx]; key <- paste(m, j, sep="|")
      tv <- if (key %in% names(TRUE_BETA_MAIN)) TRUE_BETA_MAIN[key] else 0
      s  <- pmax(seh[idx], 0, na.rm=TRUE)
      r  <- mk("main_effect", paste0("b_",j), tv,
               bh[idx], s^2, bh[idx]-1.96*s, bh[idx]+1.96*s)
      r$condition <- m; 
      r$pip <- if (idx <= length(infl)) mfit$pip_hat[idx] else NA_real_
      results[[length(results)+1]] <- r
    }

    # ── 2-way interactions
    if (p_fit >= 2L) {
      int_pos <- which(x_orders == 1L)
      for (ic in seq_along(int_pos)) {
        pos <- int_pos[ic]; cn <- x_cols[pos]
        pts <- strsplit(cn, "_x_")[[1]]; j <- pts[1]; k <- pts[2]
        k1  <- paste(m,j,k,sep="|"); k2 <- paste(m,k,j,sep="|")
        tv  <- if (k1 %in% names(TRUE_BETA_INT2)) TRUE_BETA_INT2[[k1]]
               else if (k2 %in% names(TRUE_BETA_INT2)) TRUE_BETA_INT2[[k2]] else 0
        s   <- pmax(seh[pos], 0, na.rm=TRUE)
        r   <- mk("interaction_2way", paste0("b2_",j,"_",k), tv,
                  bh[pos], s^2, bh[pos]-1.96*s, bh[pos]+1.96*s)
        r$condition <- m; 
        r$pip <- if (ic <= length(int_pos)) mfit$pip_hat[ic] else NA_real_
        results[[length(results)+1]] <- r
      }
    }

    # ── 3-way interactions
    if (p_fit >= 3L && p_true >= 3L) {
      three_pos <- which(x_orders == 2L)
      for (ic in seq_along(three_pos)) {
        pos <- three_pos[ic]; cn <- x_cols[pos]
        pts <- strsplit(cn, "_x_")[[1]]; j<-pts[1]; k<-pts[2]; l<-pts[3]
        keys <- apply(permn_3(c(j,k,l)), 1,
                      function(x) paste(c(m,x), collapse="|"))
        mk_k <- keys[keys %in% names(TRUE_BETA_INT3)]
        tv   <- if (length(mk_k)>0) TRUE_BETA_INT3[[mk_k[1]]] else 0
        s    <- pmax(seh[pos], 0, na.rm=TRUE)
        r    <- mk("interaction_3way", paste0("b3_",j,"_",k,"_",l), tv,
                   bh[pos], s^2, bh[pos]-1.96*s, bh[pos]+1.96*s)
        r$condition <- m; 
        r$pip <- if (ic <= length(three_pos)) mfit$pip_hat[ic] else NA_real_
        results[[length(results)+1]] <- r
      }
    }

    # ── Covariate coefficients  (gamma_1..gamma_P)
    # Use z_cols from mfit for accurate naming; fall back to positional names.
    z_names   <- if (!is.null(mfit$z_cols)) mfit$z_cols
                 else c("intercept","sex_male","smk_current","smk_former","age")
    cov_names <- c("intercept","sex_male","smk_current","smk_former","age")
    true_g    <- c(TRUE_BETA0[m], TRUE_GAMMA[[m]][-1], TRUE_GAMMA[[m]][1]); names(true_g) <- cov_names
    
    
    for (gi in 1:n_gamma) {
      # Map z_cols name → true_g name (z_cols may not include "intercept")
      zn  <- if (gi <= length(z_names)) z_names[gi] else paste0("gamma_",gi)
      tvg <- if (zn %in% names(true_g)) true_g[zn] else NA_real_
      s_g <- pmax(gse[gi], 0, na.rm = TRUE)
      r   <- mk("covariate", paste0("g_", zn), tvg,
                gh[gi], s_g^2, gh[gi]-1.96*s_g, gh[gi]+1.96*s_g)
      r$condition <- m
      r$pip <- if (gi <= length(pg_vec)) pg_vec[gi] else NA_real_
      results[[length(results)+1]] <- r
    }
  }

  rbindlist(results, fill = TRUE)
}

# --------------------------------------------------------------------------- #
# 7.  Predictive metrics  (identical to ctbn_simulation_v3.R; uses get_lp S3)
# --------------------------------------------------------------------------- #

compute_pred_metrics <- function(fit, DT_test, eval_times = EVAL_TIMES) {
  results <- list()
  for (m in CONDS) {
    ec  <- paste0(m,"_event")
    atr <- DT_test[get(m)==0 & dt>0]
    if (nrow(atr)==0 || !ec %in% names(atr)) next
    pred  <- tryCatch(get_lp(fit, atr, m),
                      error=function(e)
                        list(lp=rep(NA,nrow(atr)), lambda=rep(NA,nrow(atr))))
    y     <- atr[[ec]]; dt_ <- atr$dt
    lam   <- pred$lambda; lp <- pred$lp
    mu    <- lam * pmax(dt_, 1e-10)
    p_hat <- 1 - exp(-pmax(mu, 0))
    pll   <- mean(dpois(y, pmax(mu,1e-15), log=TRUE), na.rm=TRUE)
    brier <- mean((y-p_hat)^2, na.rm=TRUE)
    pd    <- data.table(event=y, t_ev=atr$time_to_event,
                         marker=lp+log(pmax(dt_,1e-10)))[is.finite(marker)]
    av    <- setNames(rep(NA_real_,length(eval_times)), as.character(eval_times))
    if (nrow(pd)>1 && length(unique(pd$event))==2) {
      td <- tryCatch(
        timeROC::timeROC(T=pd$t_ev, delta=pd$event, marker=pd$marker,
                         cause=1, times=eval_times, iid=FALSE),
        error=function(e) NULL)
      if (!is.null(td)) av <- setNames(td$AUC, as.character(eval_times))
    }
    # av already handles NA fallback and uses the correct eval_time values
    results[[m]] <- data.table(condition=m, pll=pll, brier=brier, as.list(av))
  }
  rbindlist(results, fill=TRUE)
}

# --------------------------------------------------------------------------- #
# 8.  One-replicate worker  (ctbn_map_fast replaces ctbn_fit)
# --------------------------------------------------------------------------- #

run_one_replicate <- function(r, scenario, p_fit, p_true, theta, prior,
                               k_fold, n_sim, eval_times, base_seed) {
  seed_r       <- base_seed + r * 1000L
  include_3way <- p_true >= 3L
  rep_path     <- file.path("sim_results/replicates",
                             sprintf("rep_%s_%s_%03d.rds", scenario, prior, r))
  t_start      <- proc.time()["elapsed"]

  rep_log(scenario, r, sprintf(
    "START | prior=%s | p_fit=%d p_true=%d theta=%.1f seed=%d",
    prior, p_fit, p_true, theta, seed_r))

  # ── Generate data ──────────────────────────────────────────────────────── #
  set.seed(seed_r)
  rep_log(scenario, r, "Generating dataset...")
  DT_raw  <- generate_dataset(n=n_sim, seed=seed_r, include_3way=include_3way)
  DT_wide <- prepare_wide(DT_raw)
  DT_wide$intercept <- rep(1, nrow(DT_wide))

  setorder(DT_wide, eid, time_to_event)
  for (cond in CONDS) {
    ec <- paste0(cond, "_event")
    DT_wide[, (ec) := as.numeric(
      get(cond)==0 & shift(get(cond), type="lead")==1), by=eid]
    DT_wide[is.na(get(ec)), (ec) := 0L]
  }

  # ── Full-data MAP fit: recovery + selection ────────────────────────────── #
  rep_log(scenario, r, "Fitting full-data MAP model...")
  fit_full <- tryCatch(
    ctbn_map_fast(
      DT_wide           = DT_wide,
      prior             = prior,
      max_order         = p_fit,
      fixed_covs        = c("intercept","sex_male","smk_current","smk_former"),
      time_varying_covs = "age",
      variable_select   = FALSE,
      theta             = theta,
      lbfgs_maxit       = 500L,
      compute_se        = TRUE,
      verbose           = FALSE
    ),
    error = function(e) {
      rep_log(scenario, r, conditionMessage(e), "WARN"); NULL
    }
  )

  recovery_dt <- if (!is.null(fit_full)) {
    rep_log(scenario, r, "Computing recovery metrics...")
    tryCatch(compute_recovery_metrics(fit_full, p_true=p_true),
             error=function(e) { rep_log(scenario,r,e$message,"WARN"); NULL })
  } else NULL

  selection_dt <- if (!is.null(fit_full)) {
    rep_log(scenario, r, "Computing selection metrics...")
    tryCatch(compute_selection_metrics(fit_full, p_true=p_true),
             error=function(e) { rep_log(scenario,r,e$message,"WARN"); NULL })
  } else NULL

  # ── K-fold CV: predictive + oracle ───────────────────────────────────── #
  eids     <- unique(DT_wide$eid)
  set.seed(seed_r + 1L)
  fold_ids <- sample(rep(seq_len(k_fold), length.out=length(eids)))
  fold_map <- data.table(eid=eids, fold=fold_ids)
  DT_wide  <- merge(DT_wide, fold_map, by="eid")

  pred_list <- oracle_list <- list()

  for (k in seq_len(k_fold)) {
    rep_log(scenario, r, sprintf("Fold %d/%d — fitting MAP...", k, k_fold))
    DT_tr <- DT_wide[fold != k][, fold := NULL]
    DT_te <- DT_wide[fold == k][, fold := NULL]

    fit_k <- tryCatch(
      ctbn_map_fast(
        DT_wide           = DT_tr,
        prior             = prior,
        max_order         = p_fit,
        fixed_covs        = c("intercept","sex_male","smk_current","smk_former"),
        time_varying_covs = "age",
        variable_select   = FALSE,
        theta             = theta,
        lbfgs_maxit       = 500L,
        compute_se        = FALSE,   # skip SE in CV folds — not needed
        verbose           = FALSE
      ),
      error=function(e) { rep_log(scenario,r,conditionMessage(e),"WARN"); NULL }
    )

    if (!is.null(fit_k)) {
      rep_log(scenario, r, sprintf("Fold %d/%d — computing metrics...", k, k_fold))
      pm <- compute_pred_metrics(fit_k, DT_te, eval_times)
      pm[, `:=`(fold=k, scenario=scenario, p_fit=p_fit,
                p_true=p_true, theta=theta)]
      pred_list[[k]] <- pm

      om <- compute_oracle_metrics(DT_te, include_3way=include_3way)
      om[, `:=`(fold=k, scenario=scenario, p_fit=p_fit,
                p_true=p_true, theta=theta)]
      oracle_list[[k]] <- om
    }
  }
  DT_wide[, fold := NULL]

  wall_time <- proc.time()["elapsed"] - t_start
  rep_log(scenario, r, sprintf("DONE in %s", format_duration(wall_time)))

  result <- list(
    rep          = r,
    scenario     = scenario,
    prior        = prior,
    p_fit        = p_fit,
    p_true       = p_true,
    theta        = theta,
    wall_time    = wall_time,
    recovery_dt  = recovery_dt,
    selection_dt = selection_dt,
    pred_dt      = if (length(pred_list)>0)
      rbindlist(pred_list, fill=TRUE) else NULL,
    oracle_dt    = if (length(oracle_list)>0)
      rbindlist(oracle_list, fill=TRUE) else NULL
  )
  saveRDS(result, rep_path)
  result
}

# --------------------------------------------------------------------------- #
# 9.  Main orchestrator with progress monitoring
# --------------------------------------------------------------------------- #

run_simulation <- function(scenario    = "A2",
                            n_rep      = N_REP,
                            n_sim      = N_SIM,
                            k_fold     = K_FOLD,
                            p_fit      = 2L,
                            p_true     = 2L,
                            prior      = PRIOR_FIT[1],
                            theta      = THETA_FIT,
                            n_cores    = 20L,
                            eval_times = EVAL_TIMES,
                            base_seed  = SIM_SEED,
                            poll_secs  = 5) {

  rep_paths <- file.path("sim_results/replicates",
    sprintf("rep_%s_%s_%03d.rds", scenario, prior, seq_len(n_rep)))
  todo     <- which(!file.exists(rep_paths))
  done_pre <- n_rep - length(todo)

  cat(sprintf(
    "\n%s === Scenario %s | Prior=%s | P_fit=%d | P_true=%d | theta=%.1f ===\n",
    format(Sys.time(),"%H:%M:%S"), scenario, prior, p_fit, p_true, theta))
  cat(sprintf("  Replicates: %d total | %d already done | %d to run\n\n",
              n_rep, done_pre, length(todo)))

  if (length(todo) == 0) {
    cat("  Nothing to do — all replicates already completed.\n")
    return(invisible(load_scenario(scenario, n_rep, prior)))
  }

  prog     <- SimProgress$new(scenario, n_rep, k_fold, M, prior)
  n_actual <- min(n_cores, length(todo))
  cl       <- makeCluster(n_actual,
                outfile = file.path("sim_results/logs",
                                     sprintf("cluster_%s_%s.log", scenario, prior)))
  registerDoParallel(cl)

  # All functions and constants needed by workers
  export_vars <- c(
    "CONDS","M","TRUE_BETA0","TRUE_GRAPH","TRUE_BETA_MAIN",
    "TRUE_BETA_INT2","TRUE_BETA_INT3","TRUE_GAMMA","T_HORIZON",
    "generate_dataset","prepare_wide","simulate_patient",
    "compute_q_true","sim_smoking_step","smk_to_covs",
    "compute_recovery_metrics","compute_selection_metrics",
    "compute_pred_metrics","compute_oracle_metrics",
    "run_one_replicate","permn_3","rep_log","format_duration",
    "get_lp","get_lp.ctbn_map",
    "PIP_THRESH","EVAL_TIMES","N_SIM",
    # MAP estimator and all its helpers
    "ctbn_map_fast","build_design_matrix","make_map_objective",
    "poisson_loglik","poisson_grad","poisson_neg_hess",
    "laplace_se","compute_selection_stats","compute_gamma_selection_stats",
    "init_hyperparams_eb",
    "prior_logdens_structured","grad_prior_structured","hess_diag_prior_structured",
    "prior_logdens_lasso","grad_prior_lasso","hess_diag_prior_lasso",
    "prior_logdens_spikeslab","grad_prior_spikeslab","hess_diag_prior_spikeslab",
    "prior_logdens_horseshoe","grad_prior_horseshoe","hess_diag_prior_horseshoe",
    ".hs_eff_var","hs_c2_mean","%||%"
  )

  cat(sprintf("  Starting %d parallel workers...\n\n", n_actual))

  future_results <- foreach(
    r          = todo,
    .packages  = c("data.table","timeROC","pROC","survival"),
    .export    = export_vars,
    .errorhandling = "pass"
  ) %dopar% {
    tryCatch({
      result <- run_one_replicate(
        r=r, scenario=scenario, p_fit=p_fit, p_true=p_true,
        theta=theta, prior=prior, k_fold=k_fold, n_sim=n_sim,
        eval_times=eval_times, base_seed=base_seed)
      prog$increment(r, result$wall_time)
      result
    }, error = function(e) {
      prog$record_failure(r, conditionMessage(e))
      list(rep=r, error=conditionMessage(e))
    })
  }

  stopCluster(cl); registerDoSEQ()
  cat("\n"); prog$final_report()

  all_results <- load_scenario(scenario, n_rep, prior)
  cat(sprintf("\n  Loaded %d results.\n", length(all_results)))

  summary_path <- file.path("sim_results/logs",
    sprintf("summary_%s_%s_p%d_t%02d.log", scenario, prior, p_fit, round(theta*10)))
  state <- prog$final_report()
  write(sprintf(
    "Scenario %s | Prior=%s | P_fit=%d | P_true=%d | theta=%.1f\nCompleted: %d/%d | Failed: %d",
    scenario, prior, p_fit, p_true, theta,
    state$completed + done_pre, n_rep, state$failed),
    summary_path)

  invisible(all_results)
}

# --------------------------------------------------------------------------- #
# 10. Live-progress wrapper  (identical interface to ctbn_simulation_v3.R)
# --------------------------------------------------------------------------- #

run_simulation_with_live_progress <- function(script_path,
                                               scenario    = "A2",
                                               p_fit       = 2L,
                                               p_true      = 2L,
                                               theta       = 1.0,
                                               prior       = PRIOR_FIT[1],
                                               n_rep       = N_REP,
                                               n_sim       = N_SIM,
                                               k_fold      = K_FOLD,
                                               n_cores     = 20L,
                                               eval_times  = EVAL_TIMES,
                                               base_seed   = SIM_SEED,
                                               poll_secs   = 5,
                                               map_script  = "ctbn_map_fast.R") {

  if (!requireNamespace("callr", quietly = TRUE))
    stop("Install callr first:  install.packages('callr')")

  script_path <- normalizePath(script_path, mustWork = TRUE)
  map_script  <- normalizePath(map_script,  mustWork = TRUE)

  for (d in c("sim_results","sim_results/replicates","sim_results/figures",
              "sim_results/tables","sim_results/logs","sim_results/progress"))
    dir.create(d, showWarnings = FALSE, recursive = TRUE)

  bg_args <- list(
    script_path = script_path,
    map_script  = map_script,
    scenario    = scenario,
    p_fit       = p_fit,
    p_true      = p_true,
    theta       = theta,
    prior       = prior,
    n_rep       = n_rep,
    n_sim       = n_sim,
    k_fold      = k_fold,
    n_cores     = n_cores,
    eval_times  = eval_times,
    base_seed   = base_seed
  )

  bg_func <- function(args) {
    source(args$map_script,  local = FALSE)   # loads ctbn_map_fast + helpers
    source(args$script_path, local = FALSE)   # loads this simulation script
    run_simulation(
      scenario   = args$scenario,
      n_rep      = args$n_rep,
      n_sim      = args$n_sim,
      k_fold     = args$k_fold,
      p_fit      = args$p_fit,
      p_true     = args$p_true,
      prior      = args$prior,
      theta      = args$theta,
      n_cores    = args$n_cores,
      eval_times = args$eval_times,
      base_seed  = args$base_seed
    )
  }

  bg <- callr::r_bg(
    func      = bg_func,
    args      = list(args = bg_args),
    supervise = TRUE,
    stdout    = file.path("sim_results/logs",
                           sprintf("bg_stdout_%s_%s.log", scenario, prior)),
    stderr    = file.path("sim_results/logs",
                           sprintf("bg_stderr_%s_%s.log", scenario, prior))
  )

  cat(sprintf("\n%s  Scenario %s / Prior %s launched (PID %d)\n",
              format(Sys.time(),"%H:%M:%S"), scenario, prior, bg$get_pid()))
  cat(sprintf("  MAP script  : %s\n", map_script))
  cat(sprintf("  stdout → sim_results/logs/bg_stdout_%s_%s.log\n\n",
              scenario, prior))

  prog <- SimProgress$new(scenario=scenario, n_rep=n_rep, k_fold=k_fold,
                           n_cond=M, prior=prior)
  Sys.sleep(2)

  while (bg$is_alive()) {
    prog$print_dashboard()
    Sys.sleep(poll_secs)
  }
  prog$print_dashboard(); cat("\n")

  exit_ok <- tryCatch({
    bg$get_result(); TRUE
  }, error = function(e) {
    cat(sprintf("\n  [ERROR] Background process failed: %s\n",
                conditionMessage(e)))
    cat(sprintf("  Check: sim_results/logs/bg_stderr_%s_%s.log\n",
                scenario, prior))
    FALSE
  })

  prog$final_report()
  if (exit_ok) cat("\n  Background process completed successfully.\n")
  invisible(load_scenario(scenario, n_rep, prior))
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# --------------------------------------------------------------------------- #
# 11. Load completed replicates
# --------------------------------------------------------------------------- #

load_scenario <- function(scenario, n_rep = N_REP, prior = NULL) {
  rep_dir <- "sim_results/replicates"
  if (!is.null(prior)) {
    # New naming: rep_<scenario>_<prior>_<NNN>.rds
    paths_new <- file.path(rep_dir,
      sprintf("rep_%s_%s_%03d.rds", scenario, prior, seq_len(n_rep)))
    # Legacy fallback: rep_<scenario>_<NNN>.rds (single-prior runs)
    paths_old <- file.path(rep_dir,
      sprintf("rep_%s_%03d.rds", scenario, seq_len(n_rep)))
    paths <- ifelse(file.exists(paths_new), paths_new, paths_old)
  } else {
    # No prior specified — return all replicates found for this scenario
    paths <- list.files(rep_dir,
      pattern = sprintf("^rep_%s_.*\\.rds$", scenario),
      full.names = TRUE)
    if (!length(paths)) {
      cat(sprintf("  No replicates found for scenario %s\n", scenario))
      return(list())
    }
  }
  exists_idx <- which(file.exists(paths))
  cat(sprintf("  Loading %d/%d replicates for scenario=%s prior=%s\n",
              length(exists_idx), length(paths), scenario,
              if (is.null(prior)) "ALL" else prior))
  lapply(paths[exists_idx], readRDS)
}

# --------------------------------------------------------------------------- #
# 12. Run scenarios  (uncomment as needed)
# --------------------------------------------------------------------------- #
#
# ── Blocking (simple, no live dashboard) ─────────────────────────────────── #
#
# Scenario A1: truth P=2, fit P=1  — all four priors
# for (pr in PRIOR_FIT)
#   run_simulation("B", p_fit=2L, p_true=2L, theta=1.0, prior=pr, n_cores=20L)
#
# Scenario A2: truth P=2, fit P=2  — BASE CASE
# for (pr in PRIOR_FIT)
#   run_simulation("A2", p_fit=2L, p_true=2L, theta=1.0, prior=pr, n_cores=20L)
#
# Scenario A3: truth P=3, fit P=3
# for (pr in PRIOR_FIT)
#   run_simulation("A3", p_fit=3L, p_true=3L, theta=1.0, prior=pr, n_cores=20L)
#
# Scenario C: truth P=3, fit P=2  — MISSPECIFIED
# for (pr in PRIOR_FIT)
#   run_simulation("C", p_fit=2L, p_true=3L, theta=1.0, prior=pr, n_cores=20L)
#
# Theta sensitivity:
# for (th in c(0.5, 2.0)) for (pr in PRIOR_FIT)
#   run_simulation("A2", p_fit=2L, p_true=2L, theta=th, prior=pr, n_cores=20L)
#
# ── Live dashboard (requires callr) ────────────────────────────────────────── #
#
#start_time <- Sys.time()
#res_A2 <- run_simulation_with_live_progress(
#  script_path = "ctbn_map_simulation.R",   # path to THIS file
#  map_script  = "ctbn_map_fast.R",         # path to MAP estimator
#  scenario    = "A2",
#  p_fit       = 2L,
#  p_true      = 2L,
#  theta       = 1.0,
#  prior       = "spike_slab",              # change to run other priors
#  n_cores     = 20L,
#  poll_secs   = 5
#)
#end_time <- Sys.time()
#message(sprintf("Simulation completed in %.1f minutes",
#                as.numeric(difftime(end_time, start_time, units="mins"))))
