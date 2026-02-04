# Simulate an experiment with control and treatment groups

Simulate an experiment with control and treatment groups

## Usage

``` r
simulate_experiment(distribution, params, n_per_group, effect_size)
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

## Value

List with control and treatment vectors
