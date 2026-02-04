# Find required sample size for target power

Find required sample size for target power

## Usage

``` r
find_sample_size(
  distribution,
  params,
  effect_size,
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

- effect_size:

  Effect size to detect

- target_power:

  Target power level

- alpha:

  Significance level

- test:

  Statistical test

- n_sim:

  Number of simulations per sample size

## Value

List with n and power_curve
