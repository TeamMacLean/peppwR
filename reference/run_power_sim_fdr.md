# Run FDR-aware power simulation for whole peptidome

Simulates an entire peptidome experiment with a mix of true nulls and
true alternatives, then applies Benjamini-Hochberg correction to compute
FDR-adjusted power.

## Usage

``` r
run_power_sim_fdr(
  fits,
  effect_size,
  n_per_group,
  prop_null = 0.9,
  fdr_threshold = 0.05,
  alpha = 0.05,
  test = "wilcoxon",
  n_sim = 1000
)
```

## Arguments

- fits:

  A peppwr_fits object

- effect_size:

  Fold change for treatment

- n_per_group:

  Number of samples per group

- prop_null:

  Proportion of peptides that are true nulls (no effect). Default 0.9
  (90% unchanged).

- fdr_threshold:

  FDR threshold for calling discoveries. Default 0.05.

- alpha:

  Nominal significance level (used for simulation). Default 0.05.

- test:

  Statistical test to use. Default "wilcoxon".

- n_sim:

  Number of simulation iterations. Default 1000.

## Value

FDR-adjusted power estimate (proportion of true alternatives detected)
