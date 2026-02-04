# Simulate experiment with realistic missingness

Generates control and treatment samples from a distribution, then
introduces missing values according to the specified rate and MNAR
pattern.

## Usage

``` r
simulate_with_missingness(
  distribution,
  params,
  n_per_group,
  effect_size,
  na_rate = 0,
  mnar_score = 0
)
```

## Arguments

- distribution:

  Distribution name (e.g., "norm", "gamma", "lnorm")

- params:

  List of distribution parameters

- n_per_group:

  Number of samples per group

- effect_size:

  Fold change for treatment (multiplicative effect)

- na_rate:

  Proportion of values to make NA (0-1)

- mnar_score:

  MNAR intensity: 0 = MCAR, positive = low values more likely to be
  missing. Typical values: 0-3.

## Value

List with control and treatment vectors (may contain NAs)
