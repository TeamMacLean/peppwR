# Plot power heatmap: N x effect size grid

Plot power heatmap: N x effect size grid

## Usage

``` r
plot_power_heatmap(
  distribution,
  params,
  n_range,
  effect_range,
  n_steps = 6,
  n_sim = 100,
  test = "wilcoxon",
  alpha = 0.05
)
```

## Arguments

- distribution:

  Distribution name

- params:

  Distribution parameters

- n_range:

  Range of sample sizes (vector of length 2)

- effect_range:

  Range of effect sizes (vector of length 2)

- n_steps:

  Number of grid points per dimension

- n_sim:

  Number of simulations per grid cell

- test:

  Statistical test to use

- alpha:

  Significance level

## Value

A ggplot object
