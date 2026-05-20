# CTBN-Multimorbidity — NEWS

## Version 1.0.0 (2026-05-20)

First public release accompanying the manuscript

> Olaniran, O. R., Paria, S. S., Khondoker, M., MacGregor, A., and Lewin, A.
> *Continuous-Time Bayesian Networks with Structured Shrinkage Priors for
> Modelling Multimorbidity Trajectories in Large-Scale Electronic Health
> Records.* Manuscript (under review), 2026.

### Features

* **Stan / NUTS MCMC implementation** (`R/ctbn_stan_fit.R`):
  full Bayesian inference under four order-dependent shrinkage priors —
  Hierarchical Structured (HS), Bayesian LASSO, Regularised Horseshoe,
  and (continuous) Spike-and-Slab. Posterior inclusion probabilities (PIPs)
  derived analytically per MCMC draw under the spike-and-slab; pseudo-PIPs
  derived from posterior shrinkage factors for the continuous priors.
* **Fast L-BFGS-B MAP + Laplace approximation** (`R/ctbn_map_fast.R`):
  50–200× faster than the full Stan fit, preserving the full interface
  (`ctbn_map_fast()` is a drop-in replacement for `ctbn_fit()`),
  used for the simulation study and for scalability testing on the
  UK Biobank cohort.
* `ctbn_cv_parallel()` for patient-level *K*-fold cross-validation with
  per-fold parallel execution and warm-start support.
* `compute_F_m()` for posterior-mean and plug-in cumulative incidence
  functions; verified equivalence within 5 × 10⁻⁴ in the UK Biobank
  application.
* `compute_interaction_effects()` and `compute_synergistic_excess()` for
  the new §5.4 interaction / synergy analysis.
* Per-replicate caching, resumability, live progress dashboard via R6
  (`SimProgress`), and per-scenario log files in
  `sim_results/logs/scenario_X.log`.

### Reproducibility artefacts

* All 33 main-paper figures (`paper/figures/`) and 7 main-paper LaTeX
  tables (`paper/tables/`).
* All 22 simulation figures (`paper/sim_results/figures/`) and 6 simulation
  tables (`paper/sim_results/tables/`).
* CONSORT diagram for the InflAim cohort
  (`paper/figures/inflAim_consort_diagram_v2.png`).
* Full LaTeX source for the main manuscript and supplementary material
  (under `paper/`, OUP authoring-template format).

### Known limitations

* The UK Biobank-scale Stan fit requires ~16 GB RAM and ~10 hours on
  8 cores for *M* = 10 conditions, *P* = 2. The MAP estimator is
  the recommended path for routine fitting and for the planned
  extension to the full *M* = 60 InflAim condition set.
* This release does *not* ship the UK Biobank data; see the
  [Data availability](README.md#data-availability) section of the README.
