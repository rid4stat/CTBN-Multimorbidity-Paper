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

# --------------------------------------------------------------------------- #
# 2.  Load data and CV results
# --------------------------------------------------------------------------- #

#################################### Run for all targets

DT_wide3 <- fread("final_ukb_H40.csv")
#cv_raw   <- fread("~/cv_results_cv_thetas.csv")

fxcovs <- c(
  "(Intercept)", "sexMale", "immigrantNot immigrant",
  "crp", "wbc", "platelet", "neutrophill", "lymphocyte", "glyca", "albumin",
  "ethnicityAsian", "ethnicityDo not know", "ethnicityEuropean", "ethnicityMixed/Other",
  "qualSchool", "qualVocational",
  "cur_emp_stsHomemaker", "cur_emp_stsUnemployed", "own_rent_accoRent",
  "bmi_g>30", "bmi_g25-30", "smokerPrevious", "smokerCurrent",
  "DII_quintileQ2", "DII_quintileQ3", "DII_quintileQ4", "DII_quintileQ5",
  "townsend_quintileQ2", "townsend_quintileQ3", "townsend_quintileQ4", "townsend_quintileQ5"
)
tvcovs <- c("age")


# ---------------------------------------------------------------------------
# Build fit_fns programmatically across the full grid
# ---------------------------------------------------------------------------

priors     <- c("spike_slab", "horseshoe", "lasso", "structured")
max_orders <- c(1L, 2L, 3L)
thetas     <- c(0.5, 1.0, 2.0)

# Short label helpers
prior_label <- c(spike_slab = "spike", horseshoe = "horse",
                 lasso = "lasso", structured = "struct")
theta_label <- function(th) gsub("\\.", "", formatC(th, format = "f", digits = 1))
 theta_label(0.5) -> "05"; theta_label(1.0) -> "10"; theta_label(2.0) -> "20"

# Factory function: returns a fit_fn closure for one grid cell.
# All hyperparameters are baked into the closure at creation time so
# parallel workers receive fully self-contained functions.
make_fit_fn <- function(prior, max_order, theta) {

  # Capture all values into the closure explicitly — avoids late-binding
  # issues where loop variables change after the closure is created.
  .prior     <- prior
  .max_order <- max_order
  .theta     <- theta
  .fxcovs    <- fxcovs[-1]   # drop "(Intercept)"; ctbn_map_fast adds it internally
  .tvcovs    <- tvcovs

  function(dt, ...) {
    ctbn_map_fast(
      DT_wide           = dt,
      prior             = .prior,
      fixed_covs        = .fxcovs,
      time_varying_covs = .tvcovs,
      max_order         = .max_order,
      theta             = .theta,
      variable_select   = FALSE,   # selection via PIP threshold post-fit
      lbfgs_maxit       = 500L,
      compute_se        = TRUE,    # Laplace SEs needed for CIs and covariate section
      verbose           = FALSE,
      ...   # receives target_conditions from ctbn_map_cv_compare
    )
  }
}

# Build the named list of 36 fit functions
fit_fns <- list()

for (pr in priors) {
  for (mo in max_orders) {
    for (th in thetas) {
      nm <- sprintf("%s_m%d_t%s", prior_label[pr], mo, theta_label(th))
      fit_fns[[nm]] <- make_fit_fn(prior = pr, max_order = mo, theta = th)
    }
  }
}

message(sprintf("Grid: %d models (%d priors × %d max_orders × %d thetas)",
                length(fit_fns),
                length(priors), length(max_orders), length(thetas)))
message("Models: ", paste(names(fit_fns), collapse = ", "))




library(future)
plan(multisession, workers = parallel::detectCores() - 1L)
library(progressr)
handlers(global = TRUE)
start_time <- Sys.time()
with_progress(cv_thetas <- ctbn_map_cv_compare(
  DT_wide = DT_wide3[,-c(1,16)],        # not [,-1] — keep eid
  fit_fns = fit_fns,
  targets = NULL,   # this becomes target_conditions inside
  fixed_covs = fxcovs[-1], time_varying_covs = tvcovs,
  k_folds = 5, eval_times = c(1, 2, 5, 7, 10),
  .progress = TRUE
))
end_time <- Sys.time()
message(sprintf("CV completed in %.1f minutes",
                as.numeric(difftime(end_time, start_time, units = "mins"))))
plan(sequential)   # always restore after

write.csv(cv_thetas, "cv_results_cv_thetas.csv")


DT_wide3 <- fread("final_ukb_H40.csv")
cv_raw   <- fread("~/cv_results_cv_thetas.csv")

fxcovs <- c(
  "(Intercept)", "sexMale", "immigrantNot immigrant",
  "crp", "wbc", "platelet", "neutrophill", "lymphocyte", "glyca", "albumin",
  "ethnicityAsian", "ethnicityDo not know", "ethnicityEuropean", "ethnicityMixed/Other",
  "qualSchool", "qualVocational",
  "cur_emp_stsHomemaker", "cur_emp_stsUnemployed", "own_rent_accoRent",
  "bmi_g>30", "bmi_g25-30", "smokerPrevious", "smokerCurrent",
  "DII_quintileQ2", "DII_quintileQ3", "DII_quintileQ4", "DII_quintileQ5",
  "townsend_quintileQ2", "townsend_quintileQ3", "townsend_quintileQ4", "townsend_quintileQ5"
)
tvcovs <- c("age")

# --------------------------------------------------------------------------- #
# 3.  Parse model names
# --------------------------------------------------------------------------- #
cv <- copy(cv_raw)
cv[, prior      := sub("_m.*", "", model)]
cv[, max_order  := as.integer(sub(".*_m([0-9]+)_.*", "\\1", model))]
cv[, theta      := as.numeric(paste0(
  sub(".*_t([0-9])([0-9])$", "\\1", model), ".",
  sub(".*_t([0-9])([0-9])$", "\\2", model)))]
prior_labels <- c(spike = "Spike-and-slab", horse = "Horseshoe",
                  lasso = "LASSO",          struct = "Structured")
cv[, prior_label  := prior_labels[prior]]
cv[, target_full  := full_label(target)]
cv[, target_fac   := cond_factor(target_full)]
cv[, target_short := short_label(target)]

# =============================================================================
# 5.2.1  Effect of Interaction Order
# =============================================================================
cat("\n--- 5.2.1 ---\n")


ggplot(cvv, aes(y=poisson_ll, x = theta, col = prior_label)) +
  facet_wrap(~target_full)

order_auc <- cv |> filter(theta == 1)[, .(
  mean_tdAUC = mean(tdauc,  na.rm = TRUE), se_tdAUC = sd(tdauc,  na.rm = TRUE) / sqrt(.N),
  mean_brier = mean(brier,  na.rm = TRUE), se_brier = sd(brier,  na.rm = TRUE) / sqrt(.N),
  mean_poissonll = mean(poisson_ll,  na.rm = TRUE), se_poissonll = sd(poisson_ll,  na.rm = TRUE) / sqrt(.N)
), by = .(prior_label, max_order, eval_time)]

fig_521a <- ggplot(order_auc,
                   aes(x = eval_time, y = mean_tdAUC,
                       colour = factor(max_order), group = factor(max_order))) +
  geom_ribbon(aes(ymin = mean_tdAUC - 1.96 * se_tdAUC,
                  ymax = mean_tdAUC + 1.96 * se_tdAUC,
                  fill = factor(max_order)), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  facet_wrap(~prior_label, nrow = 2) +
  scale_colour_brewer(palette = "Dark2", name = "Interaction order P") +
  scale_fill_brewer(  palette = "Dark2", name = "Interaction order P") +
  scale_x_continuous(breaks = c(1, 2, 5, 7, 10)) +
  labs(#title    = "Figure 5.2.1a \u2014 Time-dependent AUC by interaction order",
       #subtitle = "Mean \u00b1 95% CI across all conditions and CV folds; \u03b8 = 0.5",
       x = "Evaluation horizon (years)", y = "Mean time-dependent AUC") +
  theme_ctbn()
save_fig(fig_521a, "fig_521a_order_tdauc", w = 10, h = 7)

fig_521b <- ggplot(order_auc,
                   aes(x = eval_time, y = mean_brier,
                       colour = factor(max_order), group = factor(max_order))) +
  geom_ribbon(aes(ymin = mean_brier - 1.96 * se_brier,
                  ymax = mean_brier + 1.96 * se_brier,
                  fill = factor(max_order)), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  facet_wrap(~prior_label, nrow = 2) +
  scale_colour_brewer(palette = "Dark2", name = "Interaction order P") +
  scale_fill_brewer(  palette = "Dark2", name = "Interaction order P") +
  scale_x_continuous(breaks = c(1, 2, 5, 7, 10)) +
  labs(#title    = "Figure 5.2.1b \u2014 Brier score by interaction order",
       #subtitle = "Lower is better; mean \u00b1 95% CI; \u03b8 = 0.5",
       x = "Evaluation horizon (years)", y = "Mean Brier score") +
  theme_ctbn()
save_fig(fig_521b, "fig_521b_order_brier", w = 10, h = 7)

tbl_521 <- order_auc[eval_time == 2, .(
  Prior = prior_label, `Order (P)` = max_order,
  `tdAUC (mean)` = sprintf("%.3f", mean_tdAUC), `tdAUC (SE)` = sprintf("%.4f", se_tdAUC),
  `Brier (mean)` = sprintf("%.4f", mean_brier),  `Brier (SE)` = sprintf("%.5f", se_brier),
  `Poisson LL (mean)` = sprintf("%.4f", mean_poissonll),  `Poisson LL (SE)` = sprintf("%.5f", se_poissonll)
)]
setorder(tbl_521, Prior, `Order (P)`)
save_tex(xtable(tbl_521,
                caption = "Mean tdAUC, Brier score and Poisson LL at 2-year horizon by prior
                           and interaction order. Averaged over all conditions and folds.",
                label   = "tab:order_performance"),
         "tab_521_order_performance", sanitize.text.function = identity)


# =============================================================================
# 5.2.2  Sensitivity to theta
# =============================================================================
cat("\n--- 5.2.2 ---\n")

theta_auc <- cv[, .(
  mean_tdAUC = mean(tdauc, na.rm = TRUE), se_tdAUC = sd(tdauc, na.rm = TRUE) / sqrt(.N)
), by = .(prior_label, theta, max_order, eval_time)]

fig_522 <- ggplot(theta_auc,
                  aes(x = eval_time, y = mean_tdAUC,
                      colour = factor(theta), group = factor(theta))) +
  geom_ribbon(aes(ymin = mean_tdAUC - 1.96 * se_tdAUC,
                  ymax = mean_tdAUC + 1.96 * se_tdAUC,
                  fill = factor(theta)), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  facet_grid(prior_label ~ max_order,
             labeller = labeller(max_order = function(x) paste0("P = ", x))) +
  scale_colour_brewer(palette = "Set1", name = "\u03b8") +
  scale_fill_brewer(  palette = "Set1", name = "\u03b8") +
  scale_x_continuous(breaks = c(1, 2, 5, 7, 10)) +
  labs(#title    = "Figure 5.2.2 \u2014 Sensitivity of tdAUC to order penalty \u03b8",
       #subtitle = "Rows = prior; columns = interaction order P",
       x = "Evaluation horizon (years)", y = "Mean time-dependent AUC") +
  theme_ctbn()
save_fig(fig_522, "fig_522_theta_sensitivity", w = 14, h = 10)

tbl_522 <- theta_auc[eval_time == 5, .(
  Prior      = prior_label,
  `P`        = max_order,
  `θ`  = theta,
  `tdAUC`    = sprintf("%.3f", mean_tdAUC),
  `SE`       = sprintf("%.4f", se_tdAUC)
)]
setorder(tbl_522, Prior, P, `θ`)
save_tex(xtable(tbl_522,
                caption = "Mean 5-year tdAUC for each (prior, $P$, $\\theta$) combination.",
                label   = "tab:theta_sensitivity"),
         "tab_522_theta_sensitivity", sanitize.text.function = identity)


# =============================================================================
# 5.2.3  Best model
# =============================================================================
cat("\n--- 5.2.3 ---\n")

best_cfg <- cv[eval_time == 5, .(
  mean_tdAUC = mean(tdauc, na.rm = TRUE),
  se_tdAUC   = sd(tdauc, na.rm = TRUE) / sqrt(.N)
), by = .(prior_label, max_order, theta)][order(-mean_tdAUC)]

tbl_523 <- best_cfg[1:10, .(
  Rank       = seq_len(.N),
  Prior      = prior_label,
  `P`        = max_order,
  `θ`  = theta,
  `tdAUC`    = sprintf("%.3f", mean_tdAUC),
  `SE`       = sprintf("%.4f", se_tdAUC)
)]
save_tex(xtable(tbl_523,
                caption = "Top 10 model configurations ranked by mean 5-year tdAUC.",
                label   = "tab:best_model"),
         "tab_523_best_model", sanitize.text.function = identity)

# Per-condition profile of best model
best_row   <- best_cfg[1]
best_model <- cv[prior_label == best_row$prior_label &
                   max_order   == best_row$max_order   &
                   theta       == best_row$theta]

per_cond <- best_model[, .(
  mean_tdAUC = mean(tdauc, na.rm = TRUE),
  se_tdAUC   = sd(tdauc, na.rm = TRUE) / sqrt(.N)
), by = .(target_full, target_fac, eval_time)]
per_cond[, is_infl := target_full %in% INFL_LABELS]

fig_523 <- ggplot(per_cond,
                  aes(x = eval_time, y = mean_tdAUC,
                      colour = target_fac, linetype = is_infl, group = target_fac)) +
  geom_ribbon(aes(ymin = mean_tdAUC - 1.96 * se_tdAUC,
                  ymax = mean_tdAUC + 1.96 * se_tdAUC,
                  fill = target_fac), alpha = 0.10, colour = NA) +
  geom_line(linewidth = 0.9) +
  scale_colour_brewer(palette = "Paired", name = "Condition") +
  scale_fill_brewer(  palette = "Paired", name = "Condition") +
  scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "dashed"),
                        labels = c("No", "Yes"), name = "Inflammatory [Infl]") +
  scale_x_continuous(breaks = c(1, 2, 5, 7, 10)) +
  labs(#title    = sprintf("Figure 5.2.3 \u2014 Per-condition tdAUC: %s, P=%d, \u03b8=%.1f",
       #                   best_row$prior_label, best_row$max_order, best_row$theta),
       #subtitle = "Dashed = [Infl] inflammatory conditions | Ribbon = \u00b11.96 SE",
       x = "Evaluation horizon (years)", y = "Time-dependent AUC") +
  theme_ctbn() +
  theme(legend.position = "right")
save_fig(fig_523, "fig_523_best_model_profile", w = 12, h = 6)


# =============================================================================
# 5.2.4  Network: load fit object
# =============================================================================
cat("\n--- 5.2.4 ---\n")


# Fit spike-slab with MAP at best (theta, P) for final inference.
# If spike-slab is not the best prior, still produce the spike-slab network
# for interpretability — it is the only prior with explicit PIPs.
cat("\n", strrep("=", 70), "\n")
cat("5.2.4  Inferred Multimorbidity Network\n")
cat(strrep("=", 70), "\n\n")

# Configuration for network: use spike-slab at best (theta, P)
net_theta     <- 1
net_max_order <- 2

cat(sprintf("Fitting spike-slab (MCMC) at theta=%.1f, max_order=%d for network...\n",
            net_theta, net_max_order))
start_time <- Sys.time()
fit_network <-   ctbn_map_fast(
  DT_wide           = DT_wide3[,-c(1,16)],
  prior             = "spike_slab",
  fixed_covs        = fxcovs[-1],
  time_varying_covs = tvcovs,
  target_conditions = NULL,
  max_order         = net_max_order,
  theta             = net_theta,
  variable_select   = FALSE,   # selection via PIP threshold post-fit
  lbfgs_maxit       = 500L,
  compute_se        = TRUE,    # Laplace SEs needed for CIs and covariate section
  verbose           = FALSE
)
end_time <- Sys.time()
message(sprintf("CV completed in %.1f minutes",
                as.numeric(difftime(end_time, start_time, units = "mins"))))
saveRDS(fit_network, "net_rds_map.rds")




# fit_network must be loaded externally, e.g.:
# fit_network <- readRDS("fit_spike_p2_t1.rds")
# Assumed available in environment from here onward.

conds       <- rownames(fit_network$pip_matrix)
PIP_THRESH  <- 0.50
pip_net     <- fit_network$pip_matrix
beta_net    <- fit_network$beta_matrix

# ── Figure 5.2.4a: PIP heatmap ────────────────────────────────────────────── #
pip_long <- setDT(melt(pip_net, varnames = c("influencer","target"), value.name = "pip"))
pip_long[, inf_lab := full_label(influencer)]
pip_long[, tgt_lab := full_label(target)]
pip_long[, inf_fac := cond_factor(inf_lab)]
pip_long[, tgt_fac := cond_factor(tgt_lab)]
pip_long[influencer == target, pip := NA_real_]

fig_524a <- ggplot(pip_long, aes(x = tgt_fac, y = inf_fac, fill = pip)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(data = pip_long[!is.na(pip) & pip >= PIP_THRESH],
            aes(label = sprintf("%.2f", pip)),
            size = 2.5, colour = "white", fontface = "bold") +
  scale_fill_gradient2(low = "white", mid = "#9ecae1", high = "#08519c",
                       midpoint = 0.3, limits = c(0, 1), name = "PIP",
                       na.value = "grey90") +
  scale_x_discrete(position = "top") +
  labs(#title    = "Figure 5.2.4a \u2014 Posterior inclusion probability matrix",
       #subtitle = "Bold white: active edge (PIP \u2265 0.50) | Amber bar = [Infl] cluster",
       x = "Target condition", y = "Influencer condition",
       caption = "[Infl] = RS, DE, CLRD") +
  theme_ctbn(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 8),
        axis.text.y = element_text(size = 8),
        panel.grid  = element_blank(),
        plot.margin = margin(5, 5, 5, 35))
fig_524a <- add_infl_strips(fig_524a)
save_fig(fig_524a, "fig_524a_pip_heatmap", w = 12, h = 10)

# ── Figure 5.2.4b: RR heatmap ─────────────────────────────────────────────── #
rr_long <- setDT(melt(exp(beta_net), varnames = c("influencer","target"), value.name = "rr"))
rr_long <- merge(rr_long, pip_long[, .(influencer, target, pip)],
                 by = c("influencer","target"), all.x = TRUE)
rr_long[is.na(pip) | pip < PIP_THRESH | influencer == target, rr := NA_real_]
rr_long[, inf_fac := cond_factor(full_label(influencer))]
rr_long[, tgt_fac := cond_factor(full_label(target))]

fig_524b <- ggplot(rr_long, aes(x = tgt_fac, y = inf_fac, fill = rr)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(data = rr_long[!is.na(rr)],
            aes(label = sprintf("%.2f", rr)), size = 2.5, colour = "black") +
  rr_fill_scale(name = "RR") +
  scale_x_discrete(position = "top") +
  labs(#title    = "Figure 5.2.4b \u2014 Posterior mean rate ratio matrix",
       #subtitle = "Blue = RR > 1 | Red = RR < 1 | Grey = inactive (PIP < 0.50) | Amber bar = [Infl]",
       x = "Target condition", y = "Influencer condition",
       caption = "[Infl] = RS, DE, CLRD") +
  theme_ctbn(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 8),
        axis.text.y = element_text(size = 8),
        panel.grid  = element_blank(),
        plot.margin = margin(5, 5, 5, 35))
fig_524b <- add_infl_strips(fig_524b)
save_fig(fig_524b, "fig_524b_rr_heatmap", w = 12, h = 10)

# ── Baseline intensity helper ─────────────────────────────────────────────── #
compute_baseline_q1000 <- function(fit, DT) {
  ca       <- fit$call_args
  all_covs <- c(ca$fixed_covs, ca$time_varying_covs)
  results  <- list()

  for (target in names(fit$map_fits)) {
    mfit <- fit$map_fits[[target]]
    if (is.null(mfit)) next

    # MAP gamma estimates and Laplace SEs
    gamma_hat <- mfit$gamma_hat
    se_gamma  <- mfit$se_gamma %||% rep(0, length(gamma_hat))

    # Covariate names from the MAP fit (z_cols slot)
    cov_names <- if (!is.null(mfit$z_cols)) mfit$z_cols else all_covs
    if (length(gamma_hat) != length(cov_names)) next

    # Build z_vec at population means (intercept = 1, continuous at mean)
    z_vec                              <- rep(0, length(cov_names))
    names(z_vec)                       <- cov_names
    z_vec[cov_names == "(Intercept)"]  <- 1
    for (cv in setdiff(cov_names, "(Intercept)")) {
      v <- DT[[cv]]
      if (!is.null(v) && is.numeric(v)) z_vec[cv] <- mean(v, na.rm = TRUE)
    }

    # Point estimate: exp(z' gamma_hat) * 1000
    lp_hat   <- sum(z_vec * gamma_hat)
    q_mean   <- exp(lp_hat) * 1000

    # Propagated Laplace SE for CI: Var(lp) = z' diag(se^2) z = sum((z*se)^2)
    lp_se    <- sqrt(sum((z_vec * se_gamma)^2))
    q_lo     <- exp(lp_hat - 1.96 * lp_se) * 1000
    q_hi     <- exp(lp_hat + 1.96 * lp_se) * 1000

    results[[target]] <- data.table(
      target      = target,
      q_mean_1000 = q_mean,
      q_lo_1000   = q_lo,
      q_hi_1000   = q_hi
    )
  }
  rbindlist(results)
}

q_baseline <- compute_baseline_q1000(fit_network, DT_wide3[, -c(1,16)])

# ── Edge list ─────────────────────────────────────────────────────────────── #
edge_list <- rbindlist(lapply(conds, function(inf) {
  rbindlist(lapply(conds, function(tgt) {
    if (inf == tgt) return(NULL)
    pv <- pip_net[inf, tgt]; bv <- beta_net[inf, tgt]
    if (is.na(pv) || pv < PIP_THRESH) return(NULL)
    q0     <- q_baseline[target == tgt, q_mean_1000]
    q0     <- if (length(q0) == 0) NA_real_ else q0
    q_cond <- q0 * exp(bv)
    data.table(from = full_label(inf), to = full_label(tgt),
               from_short = short_label(inf), to_short = short_label(tgt),
               pip = round(pv, 3), beta = round(bv, 4), rr = round(exp(bv), 3),
               q_baseline = round(q0, 2), q_cond = round(q_cond, 2),
               edge_label = sprintf("%.1f/1k", q_cond),
               direction  = ifelse(bv > 0, "excitatory", "inhibitory"))
  }))
}))

# ── Figure 5.2.4c: Network graph ──────────────────────────────────────────── #
if (nrow(edge_list) > 0) {
  node_df <- data.frame(
    name      = COND_ORDER,
    abbr      = sub(".*\\((.+)\\)$", "\\1", COND_ORDER),
    is_infl   = COND_ORDER %in% INFL_LABELS,
    stringsAsFactors = FALSE
  )
  node_df$name_wrap <- str_wrap(node_df$name, width = 14)
  set.seed(123) # reproducible network
  g_net   <- graph_from_data_frame(edge_list, directed = TRUE, vertices = node_df)
  set.seed(123) # reproducible network
  fig_524c <- ggraph(g_net, layout = "fr") +
    geom_edge_arc(aes(colour = direction, width = pip, alpha = pip, label = edge_label),
                  arrow = arrow(length = unit(3,"mm"), type = "closed"),
                  end_cap = circle(7,"mm"), strength = 0.25,
                  angle_calc = "along", label_dodge = unit(3,"mm"),
                  label_size = 2.8, label_colour = "grey20", show.legend = TRUE) +
    geom_node_point(aes(colour = is_infl), size = 18, alpha = 0.88) +
    geom_node_text(aes(label = abbr), size = 3.0, colour = "white", fontface = "bold") +
    geom_node_label(aes(label = name_wrap), size = 2.4, colour = "grey15",
                    fill = "white", alpha = 0.7, label.size = 0, nudge_y = -0.55) +
    scale_edge_colour_manual(values = c(excitatory = "#2166ac", inhibitory = "#d73027"),
                             name = "Direction") +
    scale_edge_width(range = c(0.6, 2.8), name = "PIP") +
    scale_edge_alpha(range = c(0.45, 1),  name = "PIP") +
    scale_colour_manual(values = c(`FALSE` = "#4d4d4d", `TRUE` = INFL_COL),
                        labels = c("Non-inflammatory","Inflammatory [Infl]"),
                        name = "Condition type") +
    labs(#title    = "Figure 5.2.4c \u2014 Inferred directed multimorbidity network",
         #subtitle = "Edge label = conditional intensity (events/1,000 person-years)",
         caption  = "[Infl] = RS, DE, CLRD") +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey40", size = 9),
          legend.position = "right")
  save_fig(fig_524c, "fig_524c_network_graph", w = 13, h = 10)
}

# ── Table 5.2.4 ───────────────────────────────────────────────────────────── #
tbl_524 <- edge_list[order(-pip), .(
  Influencer = from_short, Target = to_short,
  `[Infl]*`  = ifelse(grepl("^\\[Infl\\]", from_short) |
                        grepl("^\\[Infl\\]", to_short), "*", ""),
  PIP = sprintf("%.3f", pip), `log RR` = sprintf("%.4f", beta),
  RR  = sprintf("%.3f", rr), `q0 (/1kpy)` = sprintf("%.2f", q_baseline)
)]
save_tex(xtable(tbl_524,
                caption = "Active edges (PIP $\\geq 0.50$) in the multimorbidity network.
                           $q_0$ = baseline intensity (events/1,000 person-years).
                           * = at least one node in the [Infl] cluster.",
                label   = "tab:network_edges"),
         "tab_524_network_edges", sanitize.text.function = identity)


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
  
  mfit <- fit$map_fits[[target]]
  if (is.null(mfit)) { warning("No MAP fit for: ", target); return(NA_real_) }

  # MAP point estimates directly — no MCMC draws needed
  beta_mean  <- mfit$beta_hat
  gamma_mean <- mfit$gamma_hat

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
# 5.2.5  Conditional Risk: Main Effects
#         Q1: What is the risk of a SECOND LTC given ONE existing condition?
#
#  For each pair (j, m), with all other conditions absent and covariates at
#  population means:
#    F_m(tau | X_j = 0)  — baseline risk
#    F_m(tau | X_j = 1)  — conditional risk given condition j present
#    excess risk          — F1 - F0
#    rate ratio           — exp(beta^m_{j})
# =============================================================================
cat("\n--- 5.2.5 ---\n")

conditional_risk_table <- function(fit, DT_wide, tau = 5) {
  ca        <- fit$call_args
  all_covs  <- c(ca$fixed_covs, ca$time_varying_covs)
  reserved  <- c("eid","time_to_event","dt", all_covs)
  all_conds <- setdiff(names(DT_wide), reserved)
  
  zbar <- sapply(all_covs, function(cv) {
    v <- DT_wide[[cv]]
    if (is.numeric(v)) mean(v, na.rm = TRUE)
    else as.numeric(names(which.max(table(v))))
  })
  names(zbar) <- all_covs
  
  results <- list()
  for (target in all_conds) {
    mfit <- fit$map_fits[[target]]; if (is.null(mfit)) next
    influencers <- setdiff(all_conds, target)
    # MAP beta point estimates — first n_main entries are the main effects
    bm <- mfit$beta_hat[seq_along(influencers)]
    if (length(bm) < length(influencers)) next
    names(bm) <- influencers
    
    make_ref <- function(x_val) {
      ref <- as.data.table(as.list(rep(0, length(all_conds))))
      setnames(ref, all_conds)
      ref[, (existing) := as.numeric(x_val)]
      ref[, (target)   := 0]
      ref[, eid           := 1L]
      ref[, time_to_event := 0]
      ref[, dt            := tau]
      for (cv in all_covs) ref[, (cv) := zbar[[cv]]]
      ref
    }
    
    for (existing in influencers) {
      F0 <- tryCatch(compute_F_m(fit, target, make_ref(0), tau), error = function(e) NA_real_)
      F1 <- tryCatch(compute_F_m(fit, target, make_ref(1), tau), error = function(e) NA_real_)
      results[[length(results) + 1]] <- data.table(
        existing_condition = existing, target_condition = target,
        baseline_risk      = round(F0, 4), conditional_risk = round(F1, 4),
        excess_risk        = round(F1 - F0, 4),
        rate_ratio         = round(exp(bm[existing]), 3),
        log_RR             = round(bm[existing], 3), tau_years = tau
      )
    }
  }
  rbindlist(results)
}

risk_rds <- "risk_table.rds"
if (file.exists(risk_rds)) {
  risk_tbl <- readRDS(risk_rds); cat("Loaded risk_table.\n")
} else {
  risk5    <- conditional_risk_table(fit_network, DT_wide3[, -c(1,16)], tau = 5)
  risk10   <- conditional_risk_table(fit_network, DT_wide3[, -c(1,16)], tau = 10)
  risk_tbl <- rbind(risk5, risk10)
  saveRDS(risk_tbl, risk_rds)
}

risk_tbl[, existing_lab  := full_label(existing_condition)]
risk_tbl[, target_lab    := full_label(target_condition)]
risk_tbl[, existing_fac  := cond_factor(existing_lab)]
risk_tbl[, target_fac    := cond_factor(target_lab)]
risk_tbl[, infl_flag     := fifelse(
  existing_condition %in% INFL_CONDS | target_condition %in% INFL_CONDS, "*", "")]

# ── Figure 5.2.5a: Excess risk heatmap ────────────────────────────────────── #
fig_525a <- ggplot(risk_tbl[tau_years == 5],
                   aes(x = target_fac, y = existing_fac, fill = excess_risk)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(is.na(excess_risk), "",
                               sprintf("%.3f", excess_risk))),
            size = 2.5, colour = "black") +
  scale_fill_gradient2(low = "#d73027", mid = "#f7f7f7", high = "#2166ac",
                       midpoint = 0, na.value = "grey90",
                       name = "Excess\nrisk",
                       labels = percent_format(accuracy = 0.1)) +
  scale_x_discrete(position = "top") +
  labs(#title    = "Figure 5.2.5a \u2014 5-year excess cumulative incidence (main effects)",
       #subtitle = paste0("F_m(5 | X_j=1, Z=Zbar) \u2212 F_m(5 | X_j=0, Z=Zbar)",
       #                  " | Red = increased risk | Amber bar = [Infl] cluster"),
       x = "Target condition (to develop)", y = "Existing condition (already present)",
       caption = "Covariates at population means. [Infl] = RS, DE, CLRD") +
  theme_ctbn(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 8),
        axis.text.y = element_text(size = 8),
        panel.grid  = element_blank(),
        plot.margin = margin(5, 5, 5, 35))
fig_525a <- add_infl_strips(fig_525a)
save_fig(fig_525a, "fig_525a_excess_risk_heatmap", w = 12, h = 10)

# ── Figure 5.2.5b: RR heatmap ─────────────────────────────────────────────── #
fig_525b <- ggplot(risk_tbl[tau_years == 5],
                   aes(x = target_fac, y = existing_fac, fill = rate_ratio)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(is.na(rate_ratio), "",
                               sprintf("%.2f", rate_ratio))),
            size = 2.5, colour = "black") +
  rr_fill_scale(name = "RR") +
  scale_x_discrete(position = "top") +
  labs(#title    = "Figure 5.2.5b \u2014 Rate ratio matrix: risk of second LTC",
       #subtitle = "RR = exp(\u03b2\u1d50\u2c7c) | Blue = RR > 1 | Red = RR < 1 | Amber bar = [Infl]",
       x = "Target condition", y = "Existing condition",
       caption = "[Infl] = RS, DE, CLRD") +
  theme_ctbn(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 8),
        axis.text.y = element_text(size = 8),
        panel.grid  = element_blank(),
        plot.margin = margin(5, 5, 5, 35))
fig_525b <- add_infl_strips(fig_525b)
save_fig(fig_525b, "fig_525b_rr_heatmap", w = 12, h = 10)

# ── Table 5.2.5 ───────────────────────────────────────────────────────────── #
top15 <- risk_tbl[tau_years == 5 & existing_condition != target_condition & !is.na(excess_risk)]
setorder(top15, -excess_risk)
save_tex(xtable(top15[1:min(15, .N), .(
  `Existing`          = existing_lab, `Target` = target_lab, `[Infl]*` = infl_flag,
  `Baseline risk`     = sprintf("%.1f%%", baseline_risk    * 100),
  `Conditional risk`  = sprintf("%.1f%%", conditional_risk * 100),
  `Excess risk`       = sprintf("+%.1f%%", excess_risk     * 100),
  RR                  = sprintf("%.2f", rate_ratio))],
  caption = "Top 15 condition pairs by 5-year excess cumulative incidence (main effects).
             * = at least one condition in the [Infl] cluster.",
  label   = "tab:conditional_risk"),
  "tab_525_conditional_risk", sanitize.text.function = identity)


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
    mfit <- fit$map_fits[[target]]; if (is.null(mfit)) next
    influencers <- setdiff(all_conds, target)

    # Build column names in the exact order used during model fitting
    tmp_data <- as.data.frame(matrix(0, 1, length(influencers)))
    colnames(tmp_data) <- influencers
    ir       <- build_interaction_cols(tmp_data, influencers, ca$max_order)
    int_cols <- ir$cols   # e.g. "DM_x_HTN", "DM_x_IHD", ...

    if (length(int_cols) == 0) next

    n_main <- length(influencers)
    n_int  <- length(int_cols)

    # Guard: MAP fit must have at least n_main + n_int beta coefficients
    if (length(mfit$beta_hat) < n_main + n_int) next

    # Main effect MAP estimates (positions 1..n_main in beta_hat)
    beta_main_mean <- mfit$beta_hat[seq_len(n_main)]
    names(beta_main_mean) <- influencers

    # Interaction MAP estimates (positions n_main+1 .. n_main+n_int)
    beta_int_mean <- mfit$beta_hat[n_main + seq_len(n_int)]
    names(beta_int_mean) <- int_cols

    # Interaction PIPs / pseudo-PIPs from the MAP fit:
    #   spike_slab / structured / lasso  -> pip_hat  (position-aligned with beta_hat)
    #   horseshoe                        -> 1 - kappa_hat  (shrinkage complement)
    # Both pip_hat and kappa_hat cover the full parameter vector (main + interactions).
    pip_int_mean <- if (fit$prior == "horseshoe" && !is.null(mfit$kappa_hat) &&
                        length(mfit$kappa_hat) >= n_main + n_int) {
      pv <- 1 - mfit$kappa_hat[n_main + seq_len(n_int)]
      setNames(pv, int_cols)
    } else if (!is.null(mfit$pip_hat) &&
               length(mfit$pip_hat) >= n_main + n_int) {
      pv <- mfit$pip_hat[n_main + seq_len(n_int)]
      setNames(pv, int_cols)
    } else {
      setNames(rep(NA_real_, n_int), int_cols)
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


# ── Run and cache interaction effects ─────────────────────────────────────── #
int_rds <- "interaction_effects.rds"
if (file.exists(int_rds)) {
  int_tbl <- readRDS(int_rds); cat("Loaded interaction_effects.\n")
} else {
  int_tbl <- compute_interaction_effects(fit_network, DT_wide3[, -c(1,16)],
                                         tau = 5, pip_thresh = 0.50)
  saveRDS(int_tbl, int_rds)
}

cat(sprintf("Interaction table: %d significant triplets\n", nrow(int_tbl)))

if (nrow(int_tbl) > 0) {
  
  # Add display labels
  int_tbl[, target_lab := full_label(target_condition)]
  int_tbl[, j_lab      := full_label(condition_j)]
  int_tbl[, k_lab      := full_label(condition_k)]
  int_tbl[, target_fac := cond_factor(target_lab)]
  int_tbl[, infl_flag  := fifelse(
    target_condition %in% INFL_CONDS |
      condition_j      %in% INFL_CONDS |
      condition_k      %in% INFL_CONDS, "*", "")]
  int_tbl[, triplet_label := paste0(short_label(condition_j), " + ",
                                    short_label(condition_k),
                                    " \u2192 ", short_label(target_condition))]
  
  setorder(int_tbl, -delta_F)
  
  # ── Table 5.2.6a: Full interaction decomposition ────────────────────────── #
  save_tex(xtable(int_tbl[, .(
    `Target`       = target_lab,
    `Cond. j`      = j_lab,
    `Cond. k`      = k_lab,
    `[Infl]*`      = infl_flag,
    `RR_j`         = sprintf("%.3f", RR_j),
    `RR_k`         = sprintf("%.3f", RR_k),
    `Syn. mult.`   = sprintf("%.3f", syn_multiplier),
    `Joint RR`     = sprintf("%.3f", joint_RR),
    `$\\Delta F_m(5)$` = sprintf("+%.1f%%", delta_F * 100),
    `PIP_jk`       = sprintf("%.3f", pip_jk))],
    caption = "Pairwise interaction decomposition for all significant triplets
             ($\\text{PIP}_{jk} \\geq 0.50$, 5-year horizon).
             RR$_j$ and RR$_k$ are main-effect rate ratios.
             Syn.\\ mult.\\ $= \\exp(\\hat{\\beta}^m_{jk})$ is the synergistic
             multiplier beyond individual effects.
             Joint RR $= $ RR$_j \\times$ RR$_k \\times$ Syn.\\ mult.
             $\\Delta F_m(5) = F_m(5 \\mid j,k) - F_m(5 \\mid j) -
             F_m(5 \\mid k) + F_m(5 \\mid \\varnothing)$ is the absolute
             synergistic excess cumulative incidence.
             * = at least one condition in the [Infl] cluster.",
    label   = "tab:interaction_decomp"),
    "tab_526a_interaction_decomp", sanitize.text.function = identity)
  # ── Figure 5.2.6a: Forest plot of synergistic multipliers ─────────────── #
  setorder(int_tbl, -syn_multiplier)
  int_tbl[, direction := ifelse(syn_multiplier > 1, "Super-multiplicative","Sub-multiplicative")]
  
  fig_526a <- ggplot(int_tbl,
                     aes(x = syn_multiplier,
                         y = reorder(triplet_label, syn_multiplier),
                         colour = direction, shape = infl_flag != "")) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60", linewidth = 0.8) +
    geom_point(size = 3.5) +
    scale_colour_manual(
      values = c("Super-multiplicative" = "#2166ac", "Sub-multiplicative" = "#d73027"),
      name   = "Synergy direction") +
    scale_shape_manual(
      values = c(`FALSE` = 16, `TRUE` = 18),
      name   = "[Infl] triplet", labels = c("No","Yes")) +
    scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3),
                  labels = c("0.50","0.75","1.00","1.50","2.00","3.00")) +
    labs(
      #title    = "Figure 5.2.6a \u2014 Synergistic multipliers: exp(\u03b2\u1d50\u2c7c\u2c7c)",
      #subtitle = paste0("Point > 1 (blue): joint effect exceeds product of individual effects\n",
      #                  "Point < 1 (red): sub-multiplicative interaction | ",
      #                  "Diamond (\u25c6) = [Infl] triplet | log x-axis"),
      x        = "Synergistic multiplier exp(\u03b2\u1d50\u2c7c\u2c7c) [log scale]",
      y        = NULL,
      caption  = "Only triplets with PIP_jk \u2265 0.50. [Infl] = RS, DE, CLRD"
    ) +
    theme_ctbn(base_size = 10) +
    theme(legend.position = "right",
          axis.text.y     = element_text(size = 8))
  save_fig(fig_526a, "fig_526a_synergy_forest", w = 12, h = max(5, nrow(int_tbl) * 0.45))
  
  # ── Figure 5.2.6b: Absolute synergistic excess forest plot ────────────── #
  setorder(int_tbl, -delta_F)
  n_show <- min(20L, nrow(int_tbl))
  
  fig_526b <- ggplot(int_tbl[seq_len(n_show)],
                     aes(x = delta_F * 100,
                         y = reorder(triplet_label, delta_F),
                         colour = infl_flag != "",
                         fill   = infl_flag != "")) +
    geom_col(alpha = 0.70, width = 0.65) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.6) +
    geom_text(aes(label = sprintf("+%.1f%%", delta_F * 100),
                  hjust = ifelse(delta_F >= 0, -0.1, 1.1)),
              size = 2.8, colour = "grey20") +
    scale_colour_manual(values = c(`FALSE` = "#4d78aa", `TRUE` = INFL_COL),
                        name = "[Infl] triplet", labels = c("No","Yes")) +
    scale_fill_manual(  values = c(`FALSE` = "#4d78aa", `TRUE` = INFL_COL),
                        name = "[Infl] triplet", labels = c("No","Yes")) +
    labs(
      #title    = "Figure 5.2.6b \u2014 Absolute synergistic excess: \u0394F\u2098(5)",
      #subtitle = paste0("\u0394F\u2098(5) = F\u2098(5|j,k) \u2212 F\u2098(5|j) \u2212 F\u2098(5|k)",
      #                  " + F\u2098(5|\u2205)",
      #                  " | Amber bars = [Infl] cluster triplets"),
      x        = "Synergistic excess risk \u0394F\u2098(5) (percentage points)",
      y        = NULL,
      caption  = "Top 20 triplets by absolute synergistic excess. [Infl] = RS, DE, CLRD"
    ) +
    theme_ctbn(base_size = 10) +
    theme(legend.position = "right",
          axis.text.y     = element_text(size = 8))
  save_fig(fig_526b, "fig_526b_synergy_excess", w = 12, h = 8)
  
  # ── Figure 5.2.6c: Cumulative incidence curves with synergistic excess ribbon ─ #
  # For the top 5 triplets by |delta_F|, plot four F_m(t) curves:
  #   F_baseline, F_j_only, F_k_only, F_joint, and the additive counterfactual
  # The ribbon between F_joint and the additive counterfactual visualises Delta F_m(t).
  
  compute_F_profile_grid <- function(fit, target, j, k, DT_wide,
                                     t_grid = seq(0, 10, by = 0.5), zbar) {
    all_covs  <- c(fit$call_args$fixed_covs, fit$call_args$time_varying_covs)
    all_conds <- setdiff(names(DT_wide), c("eid","time_to_event","dt", all_covs))
    
    make_profile <- function(xj, xk) {
      ref <- as.data.table(as.list(rep(0, length(all_conds))))
      setnames(ref, all_conds)
      ref[, (j)      := as.numeric(xj)]
      ref[, (k)      := as.numeric(xk)]
      ref[, (target) := 0]
      ref[, eid := 1L]; ref[, time_to_event := 0]; ref[, dt := max(t_grid)]
      for (cv in all_covs) ref[, (cv) := zbar[[cv]]]
      ref
    }
    
    profiles <- list(
      list(label = "Baseline (neither)", xj = 0, xk = 0),
      list(label = paste0(short_label(j), " only"), xj = 1, xk = 0),
      list(label = paste0(short_label(k), " only"), xj = 0, xk = 1),
      list(label = paste0("Both ", short_label(j), " & ", short_label(k)), xj = 1, xk = 1)
    )
    
    rbindlist(lapply(profiles, function(p) {
      F_t <- vapply(t_grid, function(tau_t) {
        if (tau_t == 0) return(0)
        tryCatch(
          compute_F_m(fit, target, make_profile(p$xj, p$xk), tau = tau_t),
          error = function(e) NA_real_)
      }, numeric(1))
      data.table(profile = p$label, t = t_grid, F = F_t)
    }))
  }
  
  top5 <- int_tbl[order(-abs(delta_F))][seq_len(min(5L, .N))]
  
  curve_list <- list()
  for (i in seq_len(nrow(top5))) {
    row  <- top5[i]
    zbar <- sapply(c(fxcovs, tvcovs), function(cv) {
      v <- DT_wide3[[cv]]; if (is.null(v)) return(0)
      if (is.numeric(v)) mean(v, na.rm = TRUE)
      else as.numeric(names(which.max(table(v))))
    })
    names(zbar) <- c(fxcovs, tvcovs)
    
    dt_curves <- compute_F_profile_grid(
      fit_network, row$target_condition, row$condition_j, row$condition_k,
      DT_wide3[, -c(1,16)], t_grid = seq(0, 10, by = 0.25), zbar = zbar)
    dt_curves[, triplet := row$triplet_label]
    dt_curves[, target  := row$target_condition]
    curve_list[[i]] <- dt_curves
  }
  curve_dt <- rbindlist(curve_list)
  
  # Additive counterfactual: F_add = F_j + F_k - F_baseline
  library(data.table)
  library(ggplot2)
  library(scales)
  
  # ── 1. Wide format ───────────────────────────────────────────────
  curve_wide <- setDT(dcast(
    curve_dt,
    triplet + target + t ~ profile,
    value.var = "F"
  ))
  
  # ── 2. Extract conditions ────────────────────────────────────────
  curve_wide[, c("cond_j", "cond_k") :=
               tstrsplit(triplet, " \\+ | → ", keep = 1:2)]
  
  # Clean [Infl] prefix ONLY for matching
  curve_wide[, cond_j_clean := gsub("^\\[Infl\\] ", "", cond_j)]
  curve_wide[, cond_k_clean := cond_k]
  
  # ── 3. Compute additive counterfactual ───────────────────────────
  curve_wide[, F_additive := {
    j_col <- grep(paste0(cond_j_clean, " only$"), names(.SD), value = TRUE)
    if (length(j_col) == 0) {
      j_col <- grep(paste0("\\[Infl\\] ", cond_j_clean, " only$"),
                    names(.SD), value = TRUE)
    }
    
    k_col <- grep(paste0(cond_k_clean, " only$"), names(.SD), value = TRUE)
    
    if (length(j_col) == 0 | length(k_col) == 0) return(NA_real_)
    
    .SD[[j_col[1]]] + .SD[[k_col[1]]] - .SD[["Baseline (neither)"]]
  }, by = .(triplet, t)]
  
  # ── 4. Compute joint (Both X & Y) ────────────────────────────────
  curve_wide[, joint := {
    pattern1 <- paste0("^Both ", cond_j_clean, " & ", cond_k_clean, "$")
    pattern2 <- paste0("^Both ", cond_k_clean, " & ", cond_j_clean, "$")
    pattern3 <- paste0("^Both \\[Infl\\] ", cond_j_clean, " & ", cond_k_clean, "$")
    
    both_col <- grep(paste(pattern1, pattern2, pattern3, sep = "|"),
                     names(.SD), value = TRUE)
    
    if (length(both_col) == 0) return(NA_real_)
    
    .SD[[both_col[1]]]
  }, by = .(triplet, t)]
  
  # ── 5. Ribbon bounds ─────────────────────────────────────────────
  curve_wide[, `:=`(
    ymin_val = pmin(F_additive, joint, na.rm = TRUE),
    ymax_val = pmax(F_additive, joint, na.rm = TRUE)
  )]
  
  # ── 6. Prepare long format (ONLY valid profiles) ─────────────────
  profile_cols <- grep("only$|^Both|Baseline", names(curve_wide), value = TRUE)
  
  curve_long <- setDT(melt(
    curve_wide,
    id.vars = c("triplet", "target", "t"),
    measure.vars = profile_cols,
    variable.name = "profile",
    value.name = "F"
  ))
  
  # ── 7. Profile types ─────────────────────────────────────────────
  curve_long[, profile_type := fcase(
    grepl("only$", profile), "main_effect",
    grepl("^Both", profile), "joint",
    profile == "Baseline (neither)", "baseline",
    default = "other"
  )]
  
  # ── 8. Plot ──────────────────────────────────────────────────────
  fig_526c <- ggplot() +
    
    # Ribbon (synergy)
    geom_ribbon(
      data = curve_wide[!is.na(joint)],
      aes(x = t, ymin = ymin_val, ymax = ymax_val, group = triplet),
      fill = "#f4a582", alpha = 0.45
    ) +
    
    # Profile curves
    geom_line(
      data = curve_long,
      aes(x = t, y = F,
          colour = profile,
          linetype = profile_type,
          group = interaction(triplet, profile)),
      linewidth = 0.85
    ) +
    
    # Additive counterfactual
    geom_line(
      data = curve_wide,
      aes(x = t, y = F_additive, group = triplet),
      colour = "grey50", linetype = "dotted", linewidth = 0.7
    ) +
    
    facet_wrap(~triplet, scales = "free_y", ncol = 2) +
    
    scale_colour_brewer(palette = "Set1", name = "Profile") +
    
    scale_linetype_manual(
      values = c(
        baseline = "solid",
        main_effect = "solid",
        joint = "solid",
        other = "dashed"
      ),
      guide = "none"
    ) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_x_continuous(breaks = c(0, 2, 5, 10)) +
    
    labs(
      #title    = "Figure 5.2.6c \u2014 Cumulative incidence curves with synergistic excess",
      #subtitle = paste0("Shaded ribbon = synergistic excess \u0394F\u2098(t) = joint \u2212 additive counterfactual\n",
      #                  "Dotted line = additive counterfactual F\u2098(j) + F\u2098(k) \u2212 F\u2098(\u2205)"),
      x = "Years from cohort entry", y = "Cumulative incidence F\u2098(t)",
      caption = "Top 5 triplets by |delta_F|. [Infl] = RS, DE, CLRD"
    ) +
    theme_ctbn(base_size = 10) +
    theme(legend.position = "bottom",
          strip.text      = element_text(face = "bold", size = 8))
  
  # ── 9. Save ─────────────────────────────────────────────────────
  save_fig(fig_526c, "fig_526c_synergy_curves", w = 14, h = 12)
  
} else {
  message("No significant pairwise interactions found at PIP >= 0.50. Skipping Section 5.2.6 figures.")
}


# =============================================================================
# 5.2.7  Population-Averaged Risk Stratified by Disease Burden
# =============================================================================
cat("\n--- 5.2.7 ---\n")

stratified_cumulative_incidence <- function(fit, DT_wide,
                                            t_grid = seq(0, 10, by = 0.5),
                                            n_sample = 500, seed = 42) {
  set.seed(seed)
  ca        <- fit$call_args
  all_covs  <- c(ca$fixed_covs, ca$time_varying_covs)
  reserved  <- c("eid","time_to_event","dt", all_covs)
  all_conds <- setdiff(names(DT_wide), reserved)
  setorder(DT_wide, eid, time_to_event)
  results <- list()
  
  for (target in all_conds) {
    cat(sprintf("  [strat] target: %s\n", target))
    influencers <- setdiff(all_conds, target)
    atr         <- DT_wide[get(target) == 0 & dt > 0]
    if (nrow(atr) == 0) next
    atr[, disease_burden := rowSums(.SD == 1L), .SDcols = influencers]
    atr[, burden_stratum := fcase(
      disease_burden == 0, "0 conditions",
      disease_burden == 1, "1 condition",
      disease_burden == 2, "2 conditions",
      disease_burden >= 3, "3+ conditions")]
    atr[, burden_stratum := factor(burden_stratum,
                                   levels = c("0 conditions","1 condition","2 conditions","3+ conditions"))]
    
    for (stratum in levels(atr$burden_stratum)) {
      pats  <- unique(atr[burden_stratum == stratum]$eid); if (length(pats) == 0) next
      psamp <- if (length(pats) > n_sample) sample(pats, n_sample) else pats
      F_mat <- matrix(NA_real_, length(psamp), length(t_grid))
      for (pi in seq_along(psamp)) {
        pdc <- copy(atr[eid == psamp[pi]])
        for (ti in seq_along(t_grid)) {
          tau_t <- t_grid[ti]; if (tau_t == 0) { F_mat[pi, ti] <- 0; next }
          F_mat[pi, ti] <- tryCatch(
            compute_F_m(fit, target,
                        pdc[, c("disease_burden","burden_stratum") := NULL],
                        tau = tau_t, zbar = NULL, unit = 10),
            error = function(e) NA_real_)
        }
      }
      results[[length(results) + 1]] <- data.table(
        target = target, burden_stratum = stratum, n_patients = length(psamp),
        t = t_grid,
        mean_F  = colMeans(F_mat, na.rm = TRUE),
        lower_F = apply(F_mat, 2, quantile, 0.025, na.rm = TRUE),
        upper_F = apply(F_mat, 2, quantile, 0.975, na.rm = TRUE))
    }
  }
  rbindlist(results)
}

strat_rds <- "strat_cumincidence.rds"
if (file.exists(strat_rds)) {
  strat_dt <- readRDS(strat_rds); cat("Loaded strat_cumincidence.\n")
} else {
  strat_dt <- stratified_cumulative_incidence(fit_network, DT_wide3[, -c(1,16)])
  saveRDS(strat_dt, strat_rds)
}

strat_dt[, target_lab     := full_label(target)]
strat_dt[, is_infl        := target %in% INFL_CONDS]
strat_dt[, target_fac     := cond_factor(target_lab)]
strat_dt[, burden_stratum := factor(burden_stratum,
                                    levels = c("0 conditions","1 condition","2 conditions","3+ conditions"))]

burden_pal <- c("0 conditions"  = "#2166AC", "1 condition"   = "#4DAC26",
                "2 conditions"  = "#F4A582", "3+ conditions" = "#D6604D")

fig_527_base <- ggplot(strat_dt,
                       aes(x = t, y = mean_F, colour = burden_stratum,
                           fill = burden_stratum, group = burden_stratum)) +
  geom_ribbon(aes(ymin = lower_F, ymax = upper_F), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~target_fac, scales = "free_y", ncol = 4) +
  scale_colour_manual(values = burden_pal, name = "Existing disease burden") +
  scale_fill_manual(  values = burden_pal, name = "Existing disease burden") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = c(0, 2, 5, 10)) +
  labs(title    = "Figure 5.2.7 \u2014 Cumulative incidence by disease burden at entry",
       subtitle = "Ribbon = 2.5th\u201397.5th percentile | Amber-shaded facet headers = [Infl]",
       x = "Years from cohort entry", y = "Cumulative incidence F\u2098(t)",
       caption = "[Infl] = RS, DE, CLRD") +
  theme_ctbn(base_size = 9) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold", size = 7))
save_fig_gt(colour_infl_strips(fig_527_base), "fig_527_stratified_curves", w = 16, h = 10)

tbl_527 <- dcast(
  strat_dt[abs(t - 10) < 0.01,
           .(target_lab, target_fac, is_infl, burden_stratum, mean_F)],
  target_lab + target_fac + is_infl ~ burden_stratum, value.var = "mean_F")
setorder(tbl_527, target_fac)
for (col in c("0 conditions","1 condition","2 conditions","3+ conditions"))
  if (col %in% names(tbl_527)) tbl_527[[col]] <- sprintf("%.1f%%", tbl_527[[col]] * 100)
tbl_527[, `[Infl]` := ifelse(is_infl, "*", "")]
save_tex(xtable(tbl_527,
                caption = "Ten-year cumulative incidence $F_m(10)$ by disease burden stratum.
                           * = [Infl] cluster.",
                label   = "tab:stratified_cumincidence"),
         "tab_527_stratified_cumincidence", sanitize.text.function = identity)


# =============================================================================
# 5.2.8  Survival and Cumulative Incidence Functions
# =============================================================================
cat("\n--- 5.2.8 ---\n")

overall_dt <- strat_dt[, .(
  mean_F  = weighted.mean(mean_F,  n_patients, na.rm = TRUE),
  lower_F = weighted.mean(lower_F, n_patients, na.rm = TRUE),
  upper_F = weighted.mean(upper_F, n_patients, na.rm = TRUE)
), by = .(target_lab, target_fac, is_infl, t)]

fig_528a <- ggplot(overall_dt,
                   aes(x = t, y = mean_F, colour = target_lab, fill = target_lab,
                       linetype = is_infl, group = target_lab)) +
  geom_ribbon(aes(ymin = lower_F, ymax = upper_F), alpha = 0.10, colour = NA) +
  geom_line(linewidth = 0.9) +
  scale_colour_brewer(palette = "Paired", name = "Condition") +
  scale_fill_brewer(  palette = "Paired", name = "Condition") +
  scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "dashed"),
                        name = "Inflammatory [Infl]", labels = c("No","Yes")) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = c(0, 2, 5, 7, 10)) +
  labs(title    = "Figure 5.2.8a \u2014 Population-averaged cumulative incidence F\u2098(t)",
       subtitle = "Dashed = [Infl] | Ribbon = 2.5th\u201397.5th percentile",
       x = "Years from cohort entry", y = "Cumulative incidence F\u2098(t)",
       caption = "[Infl] = RS, DE, CLRD") +
  theme_ctbn() + guides(colour = guide_legend(ncol = 2), fill = guide_legend(ncol = 2)) +
  theme(legend.position = "right")
save_fig(fig_528a, "fig_528a_cumincidence_overall", w = 12, h = 6)

surv_dt <- copy(overall_dt)
surv_dt[, `:=`(mean_S = 1 - mean_F, lower_S = 1 - upper_F, upper_S = 1 - lower_F)]

fig_528b <- ggplot(surv_dt,
                   aes(x = t, y = mean_S, colour = target_lab, fill = target_lab,
                       linetype = is_infl, group = target_lab)) +
  geom_ribbon(aes(ymin = lower_S, ymax = upper_S), alpha = 0.10, colour = NA) +
  geom_line(linewidth = 0.9) +
  scale_colour_brewer(palette = "Paired", name = "Condition") +
  scale_fill_brewer(  palette = "Paired", name = "Condition") +
  scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "dashed"),
                        name = "Inflammatory [Infl]", labels = c("No","Yes")) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_x_continuous(breaks = c(0, 2, 5, 7, 10)) +
  labs(title    = "Figure 5.2.8b \u2014 Survival functions S\u2098(t) = 1 \u2212 F\u2098(t)",
       subtitle = "Dashed = [Infl] | Ribbon = 2.5th\u201397.5th percentile",
       x = "Years from cohort entry", y = "Survival probability S\u2098(t)",
       caption = "[Infl] = RS, DE, CLRD") +
  theme_ctbn() + guides(colour = guide_legend(ncol = 2), fill = guide_legend(ncol = 2)) +
  theme(legend.position = "right")
save_fig(fig_528b, "fig_528b_survival_overall", w = 12, h = 6)

med_dt <- overall_dt[, {
  idx <- which(mean_F >= 0.50)[1]
  list(median_t = if (is.na(idx)) NA_real_ else t[idx], is_infl = is_infl[1])
}, by = .(target_lab, target_fac)]
setorder(med_dt, target_fac)
save_tex(xtable(med_dt[, .(
  Condition = target_lab, `[Infl]` = ifelse(is_infl, "*", ""),
  `Median time (yr)` = ifelse(is.na(median_t), ">10", sprintf("%.1f", median_t)))],
  caption = "Approximate median time (years) at which $F_m(t)$ reaches 50\\%.
             $>$10 = not reached within follow-up. * = [Infl] cluster.",
  label   = "tab:median_time"),
  "tab_528_median_time", sanitize.text.function = identity)


# =============================================================================
# 5.2.9  Covariate Associations
# =============================================================================
cat("\n--- 5.2.9 ---\n")

param_names <- names(fit_network$pip_cov_list[[1]])
targets     <- names(fit_network$pip_cov_list)
n_cov       <- length(param_names) 

pip_cov_mat  <- matrix(0, n_cov, length(targets), dimnames = list(param_names, targets))
gamma_mat    <- matrix(0, n_cov, length(targets), dimnames = list(param_names, targets))
lci_gamma_mat    <- matrix(0, n_cov, length(targets), dimnames = list(param_names, targets))
uci_gamma_mat    <- matrix(0, n_cov, length(targets), dimnames = list(param_names, targets))

expgamma_mat <- matrix(0, n_cov, length(targets), dimnames = list(param_names, targets))
lci_expgamma_mat <- matrix(0, n_cov, length(targets), dimnames = list(param_names, targets))
uci_expgamma_mat <- matrix(0, n_cov, length(targets), dimnames = list(param_names, targets))

for (t in targets) {
  pip_cov_mat[, t] <- fit_network$pip_cov_list[[t]] #[-1]
  mfit_t           <- fit_network$map_fits[[t]]
  if (is.null(mfit_t)) next
  # gamma_hat[1] = intercept; [-1] drops it to match param_names[-1]
  gamma_hat_t          <- mfit_t$gamma_hat #[-1]
  se_gamma_hat_t          <- mfit_t$se_gamma #[-1]
  lci_gamma_hat_t          <- gamma_hat_t - 1.96*se_gamma_hat_t
  uci_gamma_hat_t          <- gamma_hat_t + 1.96*se_gamma_hat_t
  gamma_mat[, t]       <- gamma_hat_t
  lci_gamma_mat[, t]       <- lci_gamma_hat_t
  uci_gamma_mat[, t]       <- uci_gamma_hat_t
  expgamma_mat[, t]    <- exp(gamma_hat_t)
  lci_expgamma_mat[, t]    <- exp(lci_gamma_hat_t)
  uci_expgamma_mat[, t]    <- exp(uci_gamma_hat_t)
}

mk_cov_long <- function(mat, val) {
  dt <- setDT(melt(as.data.table(mat, keep.rownames = "covariate"),
                   id.vars = "covariate", variable.name = "target", value.name = val))
  dt[, target_lab := full_label(target)]
  dt[, target_fac := cond_factor(target_lab)]
  dt[, is_infl    := target %in% INFL_CONDS]; dt
}
pip_cov_long  <- mk_cov_long(pip_cov_mat,  "pip")
gamma_long    <- mk_cov_long(gamma_mat,    "gamma")
expgamma_long <- mk_cov_long(expgamma_mat, "expgamma")
lci_expgamma_long <- mk_cov_long(lci_expgamma_mat, "lci_expgamma")
uci_expgamma_long <- mk_cov_long(uci_expgamma_mat, "uci_expgamma")

gamma_long <- merge(gamma_long,
                    pip_cov_long[,  .(covariate, target, pip)],    by = c("covariate","target"))
gamma_long <- merge(gamma_long,
                    expgamma_long[, .(covariate, target, expgamma)], by = c("covariate","target"))

gamma_long <- merge(gamma_long,
                    lci_expgamma_long[, .(covariate, target, lci_expgamma)], by = c("covariate","target"))

gamma_long <- merge(gamma_long,
                    uci_expgamma_long[, .(covariate, target, uci_expgamma)], by = c("covariate","target"))

gamma_long[pip < 0.5, `:=`(gamma = NA_real_, expgamma = NA_real_)]

fig_529a <- ggplot(pip_cov_long, aes(x = target_fac, y = covariate, fill = pip)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(data = pip_cov_long[pip >= 0.5],
            aes(label = sprintf("%.2f", pip)), size = 2.3, colour = "white", fontface = "bold") +
  scale_fill_gradient2(low = "white", mid = "#9ecae1", high = "#08519c",
                       midpoint = 0.3, limits = c(0,1), name = "PIP") +
  scale_x_discrete(position = "top") +
  labs(#title    = "Figure 5.2.9a \u2014 Covariate posterior inclusion probabilities",
       #subtitle = "Bold text: PIP \u2265 0.50 | Amber bar = [Infl] columns",
       x = "Target condition", y = "Covariate",
       caption = "[Infl] = RS, DE, CLRD") +
  theme_ctbn(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 8),
        plot.margin = margin(5, 5, 30, 5))
fig_529a <- add_infl_strip_x_only(fig_529a, n_rows = n_cov)
save_fig(fig_529a, "fig_529a_pip_covariates", w = 12, h = 11)

fig_529b <- ggplot(gamma_long, aes(x = target_fac, y = covariate, fill = expgamma)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(data = gamma_long[!is.na(expgamma)],
            aes(label = sprintf("%.2f", expgamma)), size = 2.3, colour = "black") +
  rr_fill_scale(name = "RR") +
  scale_x_discrete(position = "top") +
  labs(#title    = "Figure 5.2.9b \u2014 Covariate posterior mean relative risk",
       #subtitle = "RR = exp(\u03b3\u0302) | Blue = RR > 1 | Red = RR < 1 | Grey = PIP < 0.50",
       x = "Target condition", y = "Covariate",
       caption = "[Infl] = RS, DE, CLRD") +
  theme_ctbn(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 8),
        plot.margin = margin(5, 5, 30, 5))
fig_529b <- add_infl_strip_x_only(fig_529b, n_rows = n_cov)
save_fig(fig_529b, "fig_529b_rr_covariates", w = 12, h = 11)

forest_dt <- gamma_long[!is.na(expgamma)][order(-expgamma)]
forest_dt[, label := paste0(covariate, " \u2192 ", target_lab)]


fig_529c <- ggplot(forest_dt[seq_len(min(30, .N))],
                   aes(x = expgamma, y = reorder(label, expgamma),
                       colour = ifelse(expgamma > 1,"RR > 1","RR < 1"), shape = is_infl)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_point(size = 3.2) +
  scale_colour_manual(values = c("RR > 1" = "#2166ac","RR < 1" = "#d73027"), name = "Direction") +
  scale_shape_manual( values = c(`FALSE` = 16, `TRUE` = 18),
                      name = "[Infl] target", labels = c("No","Yes")) +
  scale_x_log10() +
  labs(#title    = "Figure 5.2.9c \u2014 Top 30 covariate effects (log RR scale)",
    #subtitle = "PIP \u2265 0.50 | Diamond = [Infl] target | Blue = RR > 1 | Red = RR < 1",
    x = "Posterior mean RR = exp(\u03b3\u0302)", y = NULL,
    caption = "[Infl] = RS, DE, CLRD") +
  theme_ctbn(base_size = 9) +
  theme(legend.position = "right", axis.text.y = element_text(size = 8))





fig_529c <- ggplot(forest_dt[seq_len(min(30, .N))],
                   aes(x = expgamma, y = reorder(label, expgamma),
                       colour = ifelse(expgamma > 1,"RR > 1","RR < 1"), shape = is_infl)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_point(size = 3.2) +
  geom_errorbar(aes(xmin = lci_expgamma, xmax = uci_expgamma),
                 height = 0.2, size = 0.6, alpha = 0.8) +
  scale_colour_manual(values = c("RR > 1" = "#2166ac","RR < 1" = "#d73027"), name = "Direction") +
  scale_shape_manual( values = c(`FALSE` = 16, `TRUE` = 18),
                      name = "[Infl] target", labels = c("No","Yes")) +
  scale_x_log10() +
  labs(#title    = "Figure 5.2.9c \u2014 Top 30 covariate effects (log RR scale)",
       #subtitle = "PIP \u2265 0.50 | Diamond = [Infl] target | Blue = RR > 1 | Red = RR < 1",
       x = "Posterior mean RR = exp(\u03b3\u0302)", y = NULL,
       caption = "[Infl] = RS, DE, CLRD") +
  theme_ctbn(base_size = 9) +
  theme(legend.position = "right", axis.text.y = element_text(size = 8))
save_fig(fig_529c, "fig_529c_covariate_forest", w = 12, h = 9)

tbl_529 <- gamma_long[!is.na(expgamma), .(
  Covariate = covariate, Target = target_lab, `[Infl]` = ifelse(is_infl, "*", ""),
  PIP = sprintf("%.3f", pip), `log RR` = sprintf("%.3f", gamma),
  RR  = sprintf("%.3f", expgamma))]
setorder(tbl_529, -RR)
save_tex(xtable(tbl_529,
                caption = "Covariate effects with PIP $\\geq 0.50$.
                           * = inflammatory target ([Infl] cluster).
                           Spike-and-slab; $\\theta=1$; $P=2$.",
                label   = "tab:covariate_effects"),
         "tab_529_covariate_effects", sanitize.text.function = identity)

# =============================================================================
# DONE
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("Figures  ->  ./figures/\n")
cat("Tables   ->  ./tables/\n")
cat(strrep("=", 70), "\n\n", sep = "")
