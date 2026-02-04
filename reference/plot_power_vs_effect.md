# Plot power vs effect size at fixed N

Plot power vs effect size at fixed N

## Usage

``` r
plot_power_vs_effect(power_result, effect_range, n_steps = 10, n_sim = NULL)
```

## Arguments

- power_result:

  A peppwr_power object

- effect_range:

  Range of effect sizes to explore

- n_steps:

  Number of effect size values to compute

- n_sim:

  Number of simulations per point

## Value

A ggplot object
