# Contributing to CTBN-Multimorbidity

Thank you for your interest in this repository. Because it accompanies
a peer-reviewed manuscript, contributions are very welcome but follow a
slightly more conservative policy than a typical open-source R package.

## Scope of accepted contributions

We actively welcome:

* **Bug fixes** in the R code, with a reproducer in `tests/testthat/`.
* **Performance improvements** that do not change estimator behaviour
  (please include a short benchmark in the pull request description).
* **Documentation improvements** — typos, clarifications, additional
  worked examples in `vignettes/`.
* **New unit tests** for existing functionality.
* **Portability fixes** (Windows / macOS / Linux, different R or Stan
  versions).

Please raise an issue *before* opening a pull request for:

* New estimators or new prior families. The four priors compared in the
  paper were selected after careful methodological consideration; a fifth
  alternative is a methodological contribution rather than a maintenance
  change.
* Refactoring of `R/ctbn_map_fast.R` or `R/ctbn_stan_fit.R` that changes
  function signatures or output slots — downstream scripts depend on the
  current interface.
* Anything affecting the paper artefacts (`paper/`).

## How to contribute code

1. **Fork** the repository and create a feature branch from `main`:
   ```bash
   git checkout -b fix/short-description
   ```
2. **Run the test suite** before pushing:
   ```r
   devtools::test()
   ```
3. **Style**: please follow the [tidyverse style guide](https://style.tidyverse.org/)
   except where the existing code disagrees — in that case, match the
   surrounding style.
4. **Commit messages**: use the imperative mood, e.g. "Fix off-by-one in
   `build_design_matrix()` for P = 3" rather than "Fixed off-by-one...".
5. **Pull request**: target `main`, describe what changed and why, and
   reference any related issue. CI must pass before review.

## How to contribute a bug report

Please include in the issue:

* R version (`sessionInfo()`)
* Operating system
* `rstan::stan_version()` if the bug is in the Stan pipeline
* A **minimal reproducible example** — ideally a small synthetic dataset
  from `scripts/simulation/01_map_simulation.R` rather than UK Biobank
  data
* The full error message and traceback

## Data policy

This repository must never contain UK Biobank participant-level data.
If you submit a PR that adds data files, please confirm in the PR
description that the data are either:

* Fully synthetic (e.g. from the simulation generator), **or**
* Aggregate summary statistics that have been approved for release by
  the UK Biobank Access Committee.

## Code of conduct

This project follows the [Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md).
By participating you agree to abide by its terms.
