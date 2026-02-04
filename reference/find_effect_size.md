# Find minimum detectable effect size

Find minimum detectable effect size

## Usage

``` r
find_effect_size(
  distribution,
  params,
  n_per_group,
  target_power,
  alpha,
  test,
  n_sim
)
```

## Arguments

- distribution:

  Distribution name

- params:

  Distribution parameters

- n_per_group:

  Sample size per group

- target_power:

  Target power level

- alpha:

  Significance level

- test:

  Statistical test

- n_sim:

  Number of simulations per effect size

## Value

List with effect_size and effect_curve
