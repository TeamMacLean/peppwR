# Simulate experiment using empirical bootstrap

Resamples from observed data instead of using parametric distributions.
Useful when distribution fitting fails or for non-parametric analysis.

## Usage

``` r
simulate_empirical(raw_data, n_per_group, effect_size)
```

## Arguments

- raw_data:

  Numeric vector of observed values

- n_per_group:

  Number of samples per group

- effect_size:

  Fold change for treatment (multiplicative effect)

## Value

List with control and treatment vectors
