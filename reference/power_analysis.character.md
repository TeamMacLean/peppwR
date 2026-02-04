# Power analysis with specified distribution (aggregate mode)

Power analysis with specified distribution (aggregate mode)

Default power analysis method (aggregate mode with defaults)

## Usage

``` r
# S3 method for class 'character'
power_analysis(
  distribution,
  params,
  effect_size = NULL,
  n_per_group = NULL,
  target_power = NULL,
  alpha = 0.05,
  test = "wilcoxon",
  find = "power",
  n_sim = 1000,
  ...
)

# Default S3 method
power_analysis(distribution, ...)
```

## Arguments

- distribution:

  Distribution name (e.g., "norm", "gamma", "lnorm")

- params:

  List of distribution parameters

- effect_size:

  Fold change to detect

- n_per_group:

  Sample size per group (required for find="power")

- target_power:

  Target power (required for find="sample_size" or find="effect_size")

- alpha:

  Significance level (default 0.05)

- test:

  Statistical test to use (default "wilcoxon")

- find:

  What to solve for: "power", "sample_size", or "effect_size"

- n_sim:

  Number of simulations (default 1000)

- ...:

  Additional arguments (ignored)

## Value

A peppwr_power object
