# Run power simulation using empirical bootstrap

Estimates power by repeatedly resampling from observed data rather than
sampling from fitted distributions.

## Usage

``` r
run_power_sim_empirical(
  raw_data,
  n_per_group,
  effect_size,
  alpha = 0.05,
  test = "wilcoxon",
  n_sim = 1000
)
```

## Arguments

- raw_data:

  Numeric vector of observed values

- n_per_group:

  Number of samples per group

- effect_size:

  Fold change for treatment

- alpha:

  Significance level

- test:

  Statistical test to use

- n_sim:

  Number of simulations

## Value

Power estimate (proportion of significant results)
