# Run power simulation

Run power simulation

## Usage

``` r
run_power_sim(
  distribution,
  params,
  n_per_group,
  effect_size,
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

- alpha:

  Significance level

- test:

  Statistical test to use ("wilcoxon", "bootstrap_t", "bayes_t",
  "rankprod")

- n_sim:

  Number of simulations

## Value

Power estimate (proportion of significant results)
