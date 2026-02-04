# Run power simulation with missingness

Estimates power accounting for realistic missing data patterns.

## Usage

``` r
run_power_sim_with_missingness(
  distribution,
  params,
  n_per_group,
  effect_size,
  na_rate = 0,
  mnar_score = 0,
  alpha = 0.05,
  test = "wilcoxon",
  n_sim = 1000
)
```

## Arguments

- distribution:

  Distribution name

- params:

  List of distribution parameters

- n_per_group:

  Number of samples per group

- effect_size:

  Fold change for treatment

- na_rate:

  Proportion of values that will be NA

- mnar_score:

  MNAR intensity (0 = MCAR)

- alpha:

  Significance level

- test:

  Statistical test to use

- n_sim:

  Number of simulations

## Value

Power estimate (proportion of significant results)
