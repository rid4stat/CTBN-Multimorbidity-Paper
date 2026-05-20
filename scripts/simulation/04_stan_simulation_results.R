# =============================================================================
# ctbn_simulation_results_v3.R  (progress monitoring edition)
# CTBN Simulation Study — Results Aggregation and LaTeX Preparation
#
# KEY CHANGES FROM v2
# -------------------
# 1. Progress monitoring for the results aggregation pipeline itself:
#    - Section-level progress bar printed to console as each block runs.
#    - Timing reported per section.
#    - Results aggregation log written to sim_results/logs/results_run.log
#
# 2. load_all_scenarios() with progress bar shows RDS loading progress.
#
# 3. stack_metric() reports per-scenario row counts to catch missing data.
#
# 4. All outputs (figures, tables) tracked in a manifest
#    sim_results/output_manifest.txt for easy checking.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(xtable)
  library(patchwork)
})

source("ctbn_all_codes.R")   # theme_ctbn() etc.
source("ctbn_simulation_v3.R") # load_scenario(), format_duration() etc.

# --------------------------------------------------------------------------- #
# 0.  Progress infrastructure for results pipeline
# --------------------------------------------------------------------------- #

RESULTS_LOG <- "sim_results/logs/results_run.log"
MANIFEST    <- "sim_results/output_manifest.txt"

results_log <- function(msg, level = "INFO") {
  line <- sprintf("[%s] [%s] %s\n",
                  format(Sys.time(), "%H:%M:%S"), level, msg)
  cat(line)
  cat(line, file = RESULTS_LOG, append = TRUE)
}

# Section timer: prints name, runs expr, reports elapsed
section <- function(name, expr) {
  cat(sprintf("\n%s %s\n%s\n",
              format(Sys.time(), "%H:%M:%S"),
              crayon::bold(name),
              strrep("\u2500", nchar(name) + 10)))
  t0  <- proc.time()["elapsed"]
  out <- force(expr)
  elapsed <- proc.time()["elapsed"] - t0
  results_log(sprintf("%s — done in %s", name, format_duration(elapsed)))
  cat(sprintf("  \u2713 Done in %s\n", format_duration(elapsed)))
  cat(file = MANIFEST,
      sprintf("%s | %s\n", format(Sys.time(), "%H:%M:%S"), name),
      append = TRUE)
  invisible(out)
}

# Simple inline progress bar for loops
print_progress <- function(i, n, label = "") {
  bar <- make_progress_bar(i, n, width = 30)
  cat(sprintf("\r  %s  %s", bar, label))
  if (i == n) cat("\n")
  utils::flush.console()
}

save_sim_fig <- function(p, name, w=10, h=6) {
  for (ext in c("pdf","png")) {
    path <- file.path("sim_results/figures", paste0(name,".",ext))
    if (ext == "pdf")
      ggsave(path, p, width=w, height=h, units="in", device=cairo_pdf)
    else
      ggsave(path, p, width=w, height=h, units="in", dpi=300)
  }
  cat(file = MANIFEST,
      sprintf("  FIGURE: sim_results/figures/%s.{pdf,png}\n", name),
      append = TRUE)
  message("  \u2192 figures/", name)
}

save_sim_tex <- function(xt, name, ...) {
  path <- file.path("sim_results/tables", paste0(name,".tex"))
  print(xt, file=path, floating=FALSE, include.rownames=FALSE, ...)
  cat(file = MANIFEST,
      sprintf("  TABLE:  sim_results/tables/%s.tex\n", name),
      append = TRUE)
  message("  \u2192 tables/", name)
}

# Initialise log and manifest
writeLines(sprintf("Results run started: %s\n%s",
                   format(Sys.time()), strrep("=",60)), RESULTS_LOG)
writeLines(sprintf("Output manifest — %s\n%s",
                   format(Sys.time()), strrep("=",60)), MANIFEST)

# --------------------------------------------------------------------------- #
# 1.  Load all scenario results
# --------------------------------------------------------------------------- #

scenario_meta <- data.table(
  scenario = c("A1","A2","A3","C"),
  p_true   = c(2L, 2L, 3L, 3L),
  p_fit    = c(1L, 2L, 3L, 2L),
  theta    = c(1.0,1.0,1.0,1.0),
  label    = c("P=1 fit (under-specified)",
               "P=2 fit (base case)",
               "P=3 fit (extended)",
               "P=2 fit, P=3 truth (misspecified)"))

theta_meta <- rbindlist(lapply(c("A2","A3"), function(sc) {
  rbindlist(lapply(c(0.5, 2.0), function(th) {
    data.table(scenario=sc,
               p_true=if(sc=="A2") 2L else 3L,
               p_fit =if(sc=="A2") 2L else 3L,
               theta=th,
               label=sprintf("P=%d, theta=%.1f",
                             if(sc=="A2")2L else 3L, th))
  }))
}))

all_meta <- rbindlist(list(scenario_meta, theta_meta), fill=TRUE)

raw_list <- section("1. Loading replicate results", {
  result <- list()
  total  <- nrow(all_meta)
  for (i in seq_len(total)) {
    sc <- all_meta$scenario[i]; pf <- all_meta$p_fit[i]; th <- all_meta$theta[i]
    key <- paste(sc, pf, round(th*10), sep="_")
    print_progress(i, total,
      sprintf("Loading Scenario %s (P=%d, theta=%.1f)", sc, pf, th))
    raw <- load_scenario(sc, N_REP)
    ok  <- vapply(raw, function(x) is.list(x) && !is.null(x$pred_dt), logical(1))
    results_log(sprintf("Scenario %s p=%d t=%.1f: %d/%d OK",
                        sc, pf, th, sum(ok), length(raw)))
    result[[key]] <- raw[ok]
  }
  result
})

# --------------------------------------------------------------------------- #
# 2.  Stack all metrics
# --------------------------------------------------------------------------- #

stack_metric <- function(field) {
  section(sprintf("  Stacking '%s'", field), {
    rbindlist(lapply(seq_len(nrow(all_meta)), function(i) {
      sc  <- all_meta$scenario[i]; pf <- all_meta$p_fit[i]; th <- all_meta$theta[i]
      key <- paste(sc, pf, round(th*10), sep="_")
      raw <- raw_list[[key]]; if (is.null(raw)) return(NULL)
      R_eff <- length(raw)
      dt <- rbindlist(lapply(seq_len(R_eff), function(r) {
        print_progress(r, R_eff, sprintf("Scen %s rep %d", sc, r))
        d <- raw[[r]][[field]]
        if (!is.null(d)) { d[, rep := r]; d } else NULL
      }), fill=TRUE)
      if (nrow(dt)==0) {
        results_log(sprintf("  [WARN] No data for %s field=%s", key, field), "WARN")
        return(NULL)
      }
      results_log(sprintf("  Stacked %s %s: %d rows from %d reps",
                          key, field, nrow(dt), R_eff))
      dt[, `:=`(scenario=sc, p_fit=pf, p_true=all_meta$p_true[i],
                theta=th, label=all_meta$label[i])]
      dt
    }), fill=TRUE)
  })
}

recovery_all  <- stack_metric("recovery_dt")
selection_all <- stack_metric("selection_dt")
pred_all      <- stack_metric("pred_dt")
oracle_all    <- stack_metric("oracle_dt")

auc_cols <- as.character(c(1, 2, 5, 7, 10))

# --------------------------------------------------------------------------- #
# 3.  Scenario key table
# --------------------------------------------------------------------------- #

section("3. Scenario key table (Table S0)", {
  save_sim_tex(xtable(all_meta[, .(
    Scenario = scenario,
    `$P^\\text{true}$` = p_true,
    `$P^\\text{fit}$`  = p_fit,
    `$\\theta$`        = theta,
    Description        = label)],
    caption = "Simulation scenario key ($R=100$ replicates, $N=5{,}000$ each).",
    label   = "tab:sim_scenarios"),
    "tab_s0_scenario_key", sanitize.text.function=identity)
})

# --------------------------------------------------------------------------- #
# 4.  Predictive AUC: primary cross-P figures
# --------------------------------------------------------------------------- #

section("4. Primary AUC comparison (Figures S0, S0b)", {

  pred_rep <- pred_all[, c(
    list(pll=mean(pll,na.rm=TRUE), brier=mean(brier,na.rm=TRUE)),
    lapply(setNames(auc_cols,auc_cols), function(col) {
      if (col %in% names(pred_all)) mean(.SD[[col]],na.rm=TRUE) else NA_real_
    })
  ), by=.(rep,condition,scenario,p_fit,p_true,theta,label)]

  auc_long <- rbindlist(lapply(auc_cols, function(col) {
    if (!col %in% names(pred_rep)) return(NULL)
    pred_rep[, .(tau=as.numeric(col),
                 mean_auc=mean(get(col),na.rm=TRUE),
                 se_auc  =sd(get(col),  na.rm=TRUE)/sqrt(.N)),
             by=.(scenario,p_fit,p_true,theta,label,condition)]
  }))

  auc_macro <- auc_long[, .(
    mean_auc = mean(mean_auc,na.rm=TRUE),
    se_auc   = sqrt(sum(se_auc^2,na.rm=TRUE))/.N
  ), by=.(scenario,p_fit,p_true,theta,label,tau)]

  # ── Figure S0: cross-P primary comparison ─────────────────────────────── #
  auc_main <- auc_macro[theta==1.0 & scenario %in% c("A1","A2","A3","C")]
  auc_main[, p_label := paste0("P=",p_fit,
    ifelse(scenario=="C"," (misspecified)",""))]

  fig_s0 <- ggplot(auc_main,
         aes(x=tau, y=mean_auc,
             ymin=mean_auc-1.96*se_auc,
             ymax=mean_auc+1.96*se_auc,
             colour=p_label, fill=p_label, group=p_label)) +
    geom_ribbon(alpha=0.12, colour=NA) +
    geom_line(linewidth=0.95) + geom_point(size=2.2) +
    geom_hline(yintercept=0.5, linetype="dashed", colour="grey50") +
    scale_colour_brewer(palette="Dark2", name="Fitted model") +
    scale_fill_brewer(  palette="Dark2", name="Fitted model") +
    scale_x_continuous(breaks=c(1,2,5,7,10)) +
    scale_y_continuous(labels=percent_format(1), limits=c(0.40,1.0)) +
    labs(title    = "Figure S0 \u2014 Macro-averaged tdAUC by interaction order P",
         subtitle = "A1=P1 fit; A2=P2 fit (base); A3=P3 fit; C=P2 fit on P3 truth\nMean \u00b1 1.96 SE | \u03b8=1.0",
         x = "Evaluation horizon (years)", y = "Macro-averaged tdAUC") +
    theme_ctbn()
  save_sim_fig(fig_s0, "fig_s0_auc_by_order", w=11, h=6)

  # ── Figure S0b: theta sensitivity ─────────────────────────────────────── #
  auc_th <- auc_macro[scenario %in% c("A2","A3")]
  auc_th[, theta_label := sprintf("\u03b8=%.1f",theta)]
  auc_th[, p_label     := paste0("P=",p_fit)]

  fig_s0b <- ggplot(auc_th,
         aes(x=tau, y=mean_auc,
             ymin=mean_auc-1.96*se_auc, ymax=mean_auc+1.96*se_auc,
             colour=theta_label, fill=theta_label, group=theta_label)) +
    geom_ribbon(alpha=0.12, colour=NA) +
    geom_line(linewidth=0.9) + geom_point(size=2) +
    geom_hline(yintercept=0.5, linetype="dashed", colour="grey50") +
    facet_wrap(~p_label, labeller=label_both) +
    scale_colour_brewer(palette="Set1", name="\u03b8") +
    scale_fill_brewer(  palette="Set1", name="\u03b8") +
    scale_x_continuous(breaks=c(1,2,5,7,10)) +
    scale_y_continuous(labels=percent_format(1), limits=c(0.40,1.0)) +
    labs(title    = "Figure S0b \u2014 Theta sensitivity of macro-averaged tdAUC",
         subtitle = "Left = A2 (P=2); Right = A3 (P=3)",
         x = "Horizon (years)", y = "Macro-averaged tdAUC") +
    theme_ctbn()
  save_sim_fig(fig_s0b, "fig_s0b_theta_sensitivity", w=12, h=5)

  assign("pred_rep", pred_rep, envir=.GlobalEnv)
  assign("auc_long", auc_long, envir=.GlobalEnv)
  assign("auc_macro", auc_macro, envir=.GlobalEnv)
})

# --------------------------------------------------------------------------- #
# 5.  Parameter recovery
# --------------------------------------------------------------------------- #

section("5. Parameter recovery (Table S1, Figure S1)", {

  PTYPE_ORDER <- c("baseline","main_effect","interaction_2way",
                   "interaction_3way","covariate")

  rec_type <- recovery_all[, .(
    R        = .N,
    bias     = mean(bias,    na.rm=TRUE),
    se_bias  = sd(bias,      na.rm=TRUE)/sqrt(.N),
    rmse     = sqrt(mean(sq_err,na.rm=TRUE)),
    coverage = mean(covered, na.rm=TRUE),
    se_cov   = sqrt(mean(covered,na.rm=TRUE)*(1-mean(covered,na.rm=TRUE))/.N)
  ), by=.(scenario,p_fit,p_true,theta,param_type)]
  rec_type[, param_type := factor(param_type, levels=PTYPE_ORDER)]

  tbl_s1 <- rec_type[theta==1.0 & scenario %in% c("A1","A2","A3"), .(
    Scenario          = scenario,
    `$P^{\\text{fit}}$` = p_fit,
    `Parameter type`  = param_type,
    `$N$`             = R,
    Bias              = sprintf("%.3f (%.3f)", bias,     se_bias),
    RMSE              = sprintf("%.3f",         rmse),
    `Coverage (95\\%)`= sprintf("%.3f (%.3f)", coverage, se_cov)
  )][order(Scenario, param_type)]
  save_sim_tex(xtable(tbl_s1,
    caption="Monte Carlo parameter recovery ($\\theta=1.0$; $R=100$).
             3-way rows only for A3. Nominal coverage = 0.95.",
    label="tab:sim_recovery"),
    "tab_s1_param_recovery", sanitize.text.function=identity)

  fig_s1 <- ggplot(
    rec_type[theta==1.0 & scenario %in% c("A1","A2","A3") & !is.na(param_type)],
    aes(x=param_type, colour=scenario)) +
    geom_pointrange(aes(y=bias,
                        ymin=bias-1.96*se_bias, ymax=bias+1.96*se_bias),
                    position=position_dodge(0.5), size=0.7) +
    geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
    facet_wrap(~paste0("Scenario ",scenario,", P=",p_fit), nrow=1) +
    scale_colour_brewer(palette="Set2", guide="none") +
    labs(title="Figure S1a \u2014 Bias by parameter type and scenario",
         x="Parameter type", y="Monte Carlo bias") +
    theme_ctbn(base_size=9) +
    theme(axis.text.x=element_text(angle=30,hjust=1))

  fig_s1b <- ggplot(
    rec_type[theta==1.0 & scenario %in% c("A1","A2","A3") & !is.na(param_type)],
    aes(x=param_type, colour=scenario)) +
    geom_pointrange(aes(y=coverage,
                        ymin=coverage-1.96*se_cov, ymax=coverage+1.96*se_cov),
                    position=position_dodge(0.5), size=0.7) +
    geom_hline(yintercept=0.95, linetype="dashed", colour="grey40") +
    facet_wrap(~paste0("Scenario ",scenario,", P=",p_fit), nrow=1) +
    scale_colour_brewer(palette="Set2", guide="none") +
    scale_y_continuous(limits=c(0.70,1.0), labels=percent_format(1)) +
    labs(title="Figure S1b \u2014 Coverage by parameter type and scenario",
         x="Parameter type", y="Empirical 95% CI coverage") +
    theme_ctbn(base_size=9) +
    theme(axis.text.x=element_text(angle=30,hjust=1))

  save_sim_fig(fig_s1 / fig_s1b, "fig_s1_param_recovery", w=14, h=9)
})

# --------------------------------------------------------------------------- #
# 6.  Variable selection
# --------------------------------------------------------------------------- #

section("6. Variable selection (Table S2, Figure S2)", {

  sel_agg <- selection_all[theta==1.0, .(
    R         = .N,
    mean_tpr  = mean(tpr,           na.rm=TRUE),
    se_tpr    = sd(tpr,             na.rm=TRUE)/sqrt(.N),
    mean_fpr  = mean(fpr,           na.rm=TRUE),
    se_fpr    = sd(fpr,             na.rm=TRUE)/sqrt(.N),
    mean_auc  = mean(selection_auc, na.rm=TRUE),
    se_auc    = sd(selection_auc,   na.rm=TRUE)/sqrt(.N)
  ), by=.(scenario,p_fit,level)]

  tbl_s2 <- sel_agg[scenario %in% c("A1","A2","A3"), .(
    Scenario            = scenario,
    `$P^{\\text{fit}}$` = p_fit,
    Level               = level,
    `$N$ reps`          = R,
    TPR                 = sprintf("%.3f (%.3f)", mean_tpr, se_tpr),
    FPR                 = sprintf("%.3f (%.3f)", mean_fpr, se_fpr),
    `Sel. AUC`          = sprintf("%.3f (%.3f)", mean_auc, se_auc)
  )][order(Scenario, Level)]
  save_sim_tex(xtable(tbl_s2,
    caption="Variable selection ($\\theta=1.0$; $R=100$; PIP $\\geq 0.50$).
             3-way row for A3 only.",
    label="tab:sim_selection"),
    "tab_s2_variable_selection", sanitize.text.function=identity)

  sel_long <- melt(
    selection_all[theta==1.0 & scenario %in% c("A1","A2","A3")],
    id.vars=c("rep","scenario","p_fit","level"),
    measure.vars=c("tpr","fpr","selection_auc"),
    variable.name="metric", value.name="value")
  sel_long[, p_label := paste0("P=",p_fit," (Scen.",scenario,")")]

  fig_s2 <- ggplot(sel_long,
         aes(x=value, fill=p_label, colour=p_label)) +
    geom_density(alpha=0.25, linewidth=0.7) +
    facet_grid(level ~ metric,
               labeller=labeller(
                 level  = c(main_effect="Main effects",
                            interaction_2way="2-way interactions",
                            interaction_3way="3-way interactions"),
                 metric = c(tpr="TPR",fpr="FPR",selection_auc="Sel. AUC"))) +
    scale_fill_brewer(  palette="Dark2", name="Scenario") +
    scale_colour_brewer(palette="Dark2", name="Scenario") +
    labs(title="Figure S2 \u2014 Variable selection distribution",
         x="Metric value", y="Density") +
    theme_ctbn(base_size=9)
  save_sim_fig(fig_s2, "fig_s2_selection", w=14, h=8)

  assign("sel_agg", sel_agg, envir=.GlobalEnv)
})

# --------------------------------------------------------------------------- #
# 7.  Predictive accuracy tables and per-condition AUC figure
# --------------------------------------------------------------------------- #

section("7. Predictive accuracy (Table S3, Figure S3)", {

  pred_agg <- pred_rep[theta==1.0, .(
    mean_pll   = mean(pll,   na.rm=TRUE),
    se_pll     = sd(pll,     na.rm=TRUE)/sqrt(.N),
    mean_brier = mean(brier, na.rm=TRUE),
    se_brier   = sd(brier,   na.rm=TRUE)/sqrt(.N)
  ), by=.(scenario,p_fit,condition)]

  tbl_s3 <- pred_agg[scenario %in% c("A1","A2","A3"), .(
    Scenario            = scenario,
    `$P^{\\text{fit}}$` = p_fit,
    Condition           = condition,
    `Pois. LL`          = sprintf("%.4f (%.4f)", mean_pll,   se_pll),
    Brier               = sprintf("%.4f (%.4f)", mean_brier, se_brier)
  )][order(Scenario, Condition)]
  save_sim_tex(xtable(tbl_s3,
    caption="Predictive accuracy ($K=5$ CV; $\\theta=1.0$; $R=100$).",
    label="tab:sim_pred"),
    "tab_s3_predictive_accuracy", sanitize.text.function=identity)

  auc_cond <- auc_long[theta==1.0 & scenario %in% c("A1","A2","A3")]
  auc_cond[, p_label := paste0("Scenario ",scenario," (P=",p_fit,")")]

  fig_s3a <- ggplot(auc_cond,
         aes(x=tau, y=mean_auc,
             ymin=mean_auc-1.96*se_auc, ymax=mean_auc+1.96*se_auc,
             colour=condition, fill=condition, group=condition)) +
    geom_ribbon(alpha=0.10, colour=NA) +
    geom_line(linewidth=0.85) +
    geom_hline(yintercept=0.5, linetype="dashed", colour="grey50") +
    facet_wrap(~p_label, nrow=1) +
    scale_colour_brewer(palette="Paired", name="Condition") +
    scale_fill_brewer(  palette="Paired", name="Condition") +
    scale_x_continuous(breaks=c(1,2,5,10)) +
    scale_y_continuous(labels=percent_format(1), limits=c(0.40,1.0)) +
    labs(title="Figure S3a \u2014 tdAUC by condition and P",
         x="Horizon (years)", y="Mean tdAUC") +
    theme_ctbn(base_size=9) +
    guides(colour=guide_legend(ncol=2), fill=guide_legend(ncol=2))
  save_sim_fig(fig_s3a, "fig_s3a_auc_curves", w=16, h=5)
})

# --------------------------------------------------------------------------- #
# 8.  Oracle efficiency
# --------------------------------------------------------------------------- #

section("8. Oracle efficiency (Table S4, Figure S4)", {

  oracle_rep <- oracle_all[theta==1.0, .(
    oracle_pll   = mean(oracle_pll,   na.rm=TRUE),
    oracle_brier = mean(oracle_brier, na.rm=TRUE)
  ), by=.(rep,condition,scenario,p_fit,p_true)]

  eff_dt <- merge(
    pred_rep[theta==1.0, .(rep,condition,scenario,p_fit,p_true,pll,brier)],
    oracle_rep, by=c("rep","condition","scenario","p_fit","p_true"))
  eff_dt[, brier_ratio := brier / oracle_brier]
  eff_dt[, pll_diff    := pll   - oracle_pll]

  eff_agg <- eff_dt[, .(
    mean_ratio = mean(brier_ratio, na.rm=TRUE),
    se_ratio   = sd(brier_ratio,   na.rm=TRUE)/sqrt(.N),
    mean_diff  = mean(pll_diff,    na.rm=TRUE),
    se_diff    = sd(pll_diff,      na.rm=TRUE)/sqrt(.N)
  ), by=.(scenario,p_fit,p_true,condition)]

  tbl_s4 <- eff_agg[scenario %in% c("A1","A2","A3","C"), .(
    Scenario              = scenario,
    `$P^{\\text{fit}}$`   = p_fit,
    `$P^{\\text{true}}$`  = p_true,
    Condition             = condition,
    `Brier ratio`         = sprintf("%.3f (%.3f)", mean_ratio, se_ratio),
    `PLL difference`      = sprintf("%.4f (%.4f)", mean_diff,  se_diff)
  )][order(Scenario, Condition)]
  save_sim_tex(xtable(tbl_s4,
    caption="Oracle efficiency ($\\theta=1.0$; $R=100$).
             Brier ratio $>1$ = efficiency loss. Scenario C = misspecified.",
    label="tab:sim_oracle"),
    "tab_s4_oracle_efficiency", sanitize.text.function=identity)

  eff_dt[, p_label := paste0("Scen.",scenario," (fit P=",p_fit,
                              ", true P=",p_true,")")]
  fig_s4 <- ggplot(eff_dt[scenario %in% c("A1","A2","A3","C")],
         aes(x=brier_ratio, fill=p_label, colour=p_label)) +
    geom_density(alpha=0.25, linewidth=0.7) +
    geom_vline(xintercept=1, linetype="dashed", colour="grey40") +
    facet_wrap(~condition, scales="free_y", ncol=5) +
    scale_fill_brewer(  palette="Dark2", name="Scenario") +
    scale_colour_brewer(palette="Dark2", name="Scenario") +
    labs(title="Figure S4 \u2014 Oracle Brier ratio by scenario",
         x="Brier ratio (estimated/oracle)", y="Density") +
    theme_ctbn(base_size=9)
  save_sim_fig(fig_s4, "fig_s4_oracle_efficiency", w=16, h=7)

  assign("eff_dt",  eff_dt,  envir=.GlobalEnv)
  assign("eff_agg", eff_agg, envir=.GlobalEnv)
})

# --------------------------------------------------------------------------- #
# 9.  Master summary table
# --------------------------------------------------------------------------- #

section("9. Master summary table (Table S5)", {

  macro_pred <- pred_rep[theta==1.0, .(
    mean_pll   = mean(pll,   na.rm=TRUE),
    se_pll     = sd(pll,     na.rm=TRUE)/sqrt(.N),
    mean_brier = mean(brier, na.rm=TRUE),
    se_brier   = sd(brier,   na.rm=TRUE)/sqrt(.N)
  ), by=.(scenario,p_fit)]

  macro_auc5 <- auc_long[theta==1.0 & tau==5, .(
    mean_auc5 = mean(mean_auc,na.rm=TRUE),
    se_auc5   = sqrt(sum(se_auc^2,na.rm=TRUE))/.N
  ), by=.(scenario,p_fit)]

  tbl_s5 <- merge(
    merge(macro_pred, macro_auc5, by=c("scenario","p_fit")),
    sel_agg[level=="main_effect",
            .(scenario,p_fit,mean_tpr,se_tpr,mean_fpr,se_fpr,mean_auc,se_auc)],
    by=c("scenario","p_fit"), all.x=TRUE
  )[scenario %in% c("A1","A2","A3","C"), .(
    Scenario            = scenario,
    `$P^{\\text{fit}}$` = p_fit,
    `Pois. LL`          = sprintf("%.4f (%.4f)", mean_pll,   se_pll),
    `Brier`             = sprintf("%.4f (%.4f)", mean_brier, se_brier),
    `tdAUC (5yr)`       = sprintf("%.3f (%.3f)", mean_auc5,  se_auc5),
    `TPR (main)`        = sprintf("%.3f (%.3f)", mean_tpr,   se_tpr),
    `FPR (main)`        = sprintf("%.3f (%.3f)", mean_fpr,   se_fpr),
    `Sel. AUC (main)`   = sprintf("%.3f (%.3f)", mean_auc,   se_auc)
  )][order(Scenario)]
  save_sim_tex(xtable(tbl_s5,
    caption="Master simulation summary ($\\theta=1.0$; $R=100$; macro-averaged).
             Predictive metrics from $K=5$ CV; selection from full-data fit.",
    label="tab:sim_master"),
    "tab_s5_master_summary", sanitize.text.function=identity)
})

# --------------------------------------------------------------------------- #
# 10.  Final summary
# --------------------------------------------------------------------------- #

section("10. Final manifest and summary", {
  cat("\n")
  results_log("Results pipeline complete.")
  cat(sprintf("\n  All outputs written.\n"))
  cat(sprintf("  Figures : sim_results/figures/\n"))
  cat(sprintf("  Tables  : sim_results/tables/\n"))
  cat(sprintf("  Log     : %s\n", RESULTS_LOG))
  cat(sprintf("  Manifest: %s\n", MANIFEST))
  cat("\n  Output manifest:\n")
  cat(readLines(MANIFEST), sep="\n")
})
