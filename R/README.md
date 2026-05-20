# `R/` — Core library code

This folder contains the **reusable library code** that the entry-point
scripts in `../scripts/` source. Each file is self-contained and can be
used either via `source()` (stand-alone mode) or through
`devtools::load_all()` once the package skeleton is wired up.

| File | What it provides |
|---|---|
| `ctbn_map_fast.R` | The fast L-BFGS-B MAP + Laplace estimator, all four order-dependent shrinkage priors, `build_design_matrix()`, `compute_selection_stats()`, S3 methods for `"ctbn_map"`, and `ctbn_map_parallel()`. |
| `ctbn_stan_fit.R` | The full Stan / NUTS estimator (`ctbn_fit()`), embedded Stan model strings (also available in stand-alone form in `inst/stan/`), patient-level cross-validation (`ctbn_cv_compare()`), per-fold metric functions, `compute_F_m()`, `compute_interaction_effects()`, and plotting helpers. |
| `data.R` | Lazy loader (`load_toy_DT_wide()`) for the bundled toy dataset. |
| `zzz.R` | `.onLoad()` defaults and the package start-up banner. |

The two estimator files have **interchangeable interfaces**: a script that
calls `ctbn_map_fast(DT_wide, ...)` will work unchanged after swapping
the call to `ctbn_fit(DT_wide, ...)`. The output objects share the same
slots (`beta_hat`, `gamma_hat`, `pip_hat`, `se_beta`, …), so the
downstream `compute_F_m()`, `compute_interaction_effects()`, and metric
functions in `ctbn_stan_fit.R` dispatch transparently across both.
