# Package index

## Core Functions

Main functions for distribution fitting and power analysis

- [`fit_distributions()`](https://danmaclean.github.io/peppwR/reference/fit_distributions.md)
  : Fit distributions to peptide abundance data
- [`power_analysis()`](https://danmaclean.github.io/peppwR/reference/power_analysis.md)
  : Power analysis for peptide abundance data
- [`power_analysis(`*`<character>`*`)`](https://danmaclean.github.io/peppwR/reference/power_analysis.character.md)
  [`power_analysis(`*`<default>`*`)`](https://danmaclean.github.io/peppwR/reference/power_analysis.character.md)
  : Power analysis with specified distribution (aggregate mode)
- [`power_analysis(`*`<peppwr_fits>`*`)`](https://danmaclean.github.io/peppwR/reference/power_analysis.peppwr_fits.md)
  : Power analysis for per-peptide mode using fitted distributions

## Diagnostic Plots

Visualization functions for fits and power results

- [`plot_density_overlay()`](https://danmaclean.github.io/peppwR/reference/plot_density_overlay.md)
  : Plot density overlay: observed histogram with fitted density curve
- [`plot_missingness()`](https://danmaclean.github.io/peppwR/reference/plot_missingness.md)
  : Plot missingness statistics
- [`plot_param_distribution()`](https://danmaclean.github.io/peppwR/reference/plot_param_distribution.md)
  : Plot distribution of fitted parameters across peptidome
- [`plot_power_heatmap()`](https://danmaclean.github.io/peppwR/reference/plot_power_heatmap.md)
  : Plot power heatmap: N x effect size grid
- [`plot_power_vs_effect()`](https://danmaclean.github.io/peppwR/reference/plot_power_vs_effect.md)
  : Plot power vs effect size at fixed N
- [`plot_qq()`](https://danmaclean.github.io/peppwR/reference/plot_qq.md)
  : Plot QQ plots for goodness-of-fit
- [`plot(`*`<peppwr_fits>`*`)`](https://danmaclean.github.io/peppwR/reference/plot.peppwr_fits.md)
  : Plot method for peppwr_fits
- [`plot(`*`<peppwr_power>`*`)`](https://danmaclean.github.io/peppwR/reference/plot.peppwr_power.md)
  : Plot method for peppwr_power

## Simulation Functions

Functions for running power simulations

- [`run_power_sim()`](https://danmaclean.github.io/peppwR/reference/run_power_sim.md)
  : Run power simulation
- [`run_power_sim_empirical()`](https://danmaclean.github.io/peppwR/reference/run_power_sim_empirical.md)
  : Run power simulation using empirical bootstrap
- [`run_power_sim_fdr()`](https://danmaclean.github.io/peppwR/reference/run_power_sim_fdr.md)
  : Run FDR-aware power simulation for whole peptidome
- [`run_power_sim_with_missingness()`](https://danmaclean.github.io/peppwR/reference/run_power_sim_with_missingness.md)
  : Run power simulation with missingness
- [`simulate_empirical()`](https://danmaclean.github.io/peppwR/reference/simulate_empirical.md)
  : Simulate experiment using empirical bootstrap
- [`simulate_experiment()`](https://danmaclean.github.io/peppwR/reference/simulate_experiment.md)
  : Simulate an experiment with control and treatment groups
- [`simulate_with_missingness()`](https://danmaclean.github.io/peppwR/reference/simulate_with_missingness.md)
  : Simulate experiment with realistic missingness

## Missingness

Functions for handling missing data

- [`compute_missingness()`](https://danmaclean.github.io/peppwR/reference/compute_missingness.md)
  : Compute missingness statistics for a vector of values
- [`get_mnar_peptides()`](https://danmaclean.github.io/peppwR/reference/get_mnar_peptides.md)
  : Get peptides showing MNAR pattern

## Statistical Tests

Test functions used in power simulations

- [`test_bayes_t()`](https://danmaclean.github.io/peppwR/reference/test_bayes_t.md)
  : Bayes factor t-test
- [`test_bootstrap_t()`](https://danmaclean.github.io/peppwR/reference/test_bootstrap_t.md)
  : Bootstrap-t test
- [`test_rankprod()`](https://danmaclean.github.io/peppwR/reference/test_rankprod.md)
  : Rank Products test
- [`test_wilcoxon()`](https://danmaclean.github.io/peppwR/reference/test_wilcoxon.md)
  : Wilcoxon rank-sum test

## Classes and Methods

S3 classes and their methods

- [`new_peppwr_fits()`](https://danmaclean.github.io/peppwR/reference/new_peppwr_fits.md)
  : Create a new peppwr_fits object
- [`new_peppwr_power()`](https://danmaclean.github.io/peppwR/reference/new_peppwr_power.md)
  : Create a new peppwr_power object
- [`print(`*`<peppwr_fits>`*`)`](https://danmaclean.github.io/peppwR/reference/print.peppwr_fits.md)
  : Print method for peppwr_fits
- [`print(`*`<peppwr_power>`*`)`](https://danmaclean.github.io/peppwR/reference/print.peppwr_power.md)
  : Print method for peppwr_power
- [`print(`*`<summary.peppwr_fits>`*`)`](https://danmaclean.github.io/peppwR/reference/print.summary.peppwr_fits.md)
  : Print method for summary.peppwr_fits
- [`print(`*`<summary.peppwr_power>`*`)`](https://danmaclean.github.io/peppwR/reference/print.summary.peppwr_power.md)
  : Print method for summary.peppwr_power
- [`summary(`*`<peppwr_fits>`*`)`](https://danmaclean.github.io/peppwR/reference/summary.peppwr_fits.md)
  : Summary method for peppwr_fits
- [`summary(`*`<peppwr_power>`*`)`](https://danmaclean.github.io/peppwR/reference/summary.peppwr_power.md)
  : Summary method for peppwr_power
- [`validate_peppwr_fits()`](https://danmaclean.github.io/peppwR/reference/validate_peppwr_fits.md)
  : Validate a peppwr_fits object
- [`validate_peppwr_power()`](https://danmaclean.github.io/peppwR/reference/validate_peppwr_power.md)
  : Validate a peppwr_power object

## Utilities

Helper functions

- [`get_distribution_preset()`](https://danmaclean.github.io/peppwR/reference/get_distribution_preset.md)
  : Get distributions for a preset
- [`is_count_data()`](https://danmaclean.github.io/peppwR/reference/is_count_data.md)
  : Check if data appears to be count data (non-negative integers)
