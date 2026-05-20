# `paper/` — Manuscript artefacts

This folder contains the LaTeX source of the manuscript plus every figure
and table that appears in the paper or supplementary material.

## Contents

| Path | Description |
|---|---|
| `main.tex`                       | Full manuscript source (OUP authoring template). |
| `supplementary.tex`              | Supplementary material (prior derivations, additional results). |
| `reference.bib`                  | BibTeX bibliography. |
| `manifest.txt`                   | OUP manifest. |
| `figures/`                       | 33 main-paper figures (PNG + PDF) plus the InflAim CONSORT diagram. |
| `tables/`                        | 7 main-paper LaTeX tables. |
| `sim_results/figures/`           | 22 simulation figures (PNG + PDF). |
| `sim_results/tables/`            | 6 simulation tables. |

## How to (re)build the paper

```bash
cd paper
latexmk -pdf main.tex
latexmk -pdf supplementary.tex
```

## How figures and tables map back to scripts

* All `fig_52X_*.{png,pdf}` and `tab_52X_*.tex` files in `figures/` /
  `tables/` are produced by
  [`scripts/ukb_analysis/04_stan_results.R`](../scripts/ukb_analysis/04_stan_results.R)
  (or its MAP-based counterpart
  [`03_map_results.R`](../scripts/ukb_analysis/03_map_results.R)).
* All `fig*` and `tab_s*` files in `sim_results/` are produced by
  [`scripts/simulation/04_stan_simulation_results.R`](../scripts/simulation/04_stan_simulation_results.R)
  (Stan path, used for the paper) or
  [`02_map_simulation_results.R`](../scripts/simulation/02_map_simulation_results.R)
  (MAP path, faster).
* The CONSORT diagram `figures/inflAim_consort_diagram_v2.png` is
  produced from the cohort-derivation logic described in §4 of the
  paper; the diagram itself was authored in a separate tool and is
  shipped as a static image.

## Licence

The LaTeX artefacts in this folder are based on the OUP authoring template
and follow the OUP authoring-template licence terms. They are provided
for the purposes of academic review and reproduction of the manuscript
results, not for general redistribution. See the project [`LICENSE`](../LICENSE)
for the full licence policy.
