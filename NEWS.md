# peppwR 0.1.0

Initial CRAN release.

## Core Features

* `fit_distributions()` fits candidate distributions (gamma, lognormal, normal,
  inverse Gaussian, inverse gamma) to peptide abundance data
* `power_analysis()` performs power analysis in two modes:

  - **Aggregate mode**: Specify distribution and parameters directly
  - **Per-peptide mode**: Use fitted distributions from pilot data

* Three analysis questions supported via `find` parameter:
  - `find = "sample_size"`: What N do I need for target power?
  - `find = "power"`: What's my power at given N?
  - `find = "effect_size"`: What's the minimum detectable effect?

## Statistical Tests

* Wilcoxon rank-sum test (`test = "wilcoxon"`, default)
* Bootstrap-t test (`test = "bootstrap_t"`)
* Bayes factor t-test (`test = "bayes_t"`)

## Missing Data Handling

* `compute_missingness()` calculates NA rates and MNAR scores per peptide
* MNAR (Missing Not At Random) detection for low-abundance dropouts
* `simulate_with_missingness()` incorporates missing data patterns in power simulations

## FDR-Aware Mode

* `apply_fdr = TRUE` in per-peptide mode simulates whole-peptidome experiments
  with Benjamini-Hochberg correction
* Configurable `prop_null` for expected proportion of true nulls

## Diagnostic Plots

* `plot_density_overlay()`: Observed histogram with fitted density curve
* `plot_qq()`: QQ plots for goodness-of-fit assessment
* `plot_power_heatmap()`: N x effect size power lookup grid
* `plot_power_vs_effect()`: Power sensitivity at fixed N
* `plot_param_distribution()`: Distribution of fit quality across peptidome
* `plot_missingness()`: NA rate and MNAR score distributions

## Empirical Bootstrap

* `on_fit_failure = "empirical"` option uses bootstrap resampling when
  parametric fitting fails
