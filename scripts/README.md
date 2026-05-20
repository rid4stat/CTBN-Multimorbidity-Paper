# `scripts/` — entry-point pipelines

This folder holds the **runnable** scripts that reproduce the simulation
study and UK Biobank analysis. The reusable library code that they
`source()` lives in `../R/`.

Run order is implied by the leading digits in each file name:

| File | Purpose | Calls library code from |
|---|---|---|
| `simulation/01_map_simulation.R`       | Generates synthetic CTBN data and runs the fast MAP-based simulation across priors, P ∈ {1,2,3}, θ ∈ {0.5,1,2}. Caches each replicate to `results/sim_results/replicates/`. | `R/ctbn_map_fast.R` |
| `simulation/02_map_simulation_results.R` | Post-processing of MAP simulation: recovery, selection, and prediction metrics; plots and LaTeX tables. | — |
| `simulation/03_stan_simulation.R`       | Full Stan / NUTS simulation study reported in the paper (hours). Same design as the MAP version, but with full posterior inference. | `R/ctbn_stan_fit.R` |
| `simulation/04_stan_simulation_results.R` | Post-processing of the Stan simulation results. | — |
| `ukb_analysis/03_map_results.R`         | UK Biobank empirical results pipeline (fast MAP), all sections 5.2.1–5.2.9. | `R/ctbn_map_fast.R` |
| `ukb_analysis/04_stan_results.R`        | **Paper pipeline** — UK Biobank empirical results with full Stan. | `R/ctbn_stan_fit.R` |
| `ukb_analysis/05_stan_results_alt.R`    | An alternative Stan-results pipeline kept for traceability; not required for reproduction. | `R/ctbn_stan_fit.R` |

The MAP and Stan paths produce **identical interfaces** (same output
slots, same metric functions) — the MAP path is recommended for routine
work and the Stan path is run once for the paper.
