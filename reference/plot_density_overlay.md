# Plot density overlay: observed histogram with fitted density curve

Plot density overlay: observed histogram with fitted density curve

## Usage

``` r
plot_density_overlay(fits, peptide_id = NULL, n_overlay = 6)
```

## Arguments

- fits:

  A peppwr_fits object

- peptide_id:

  Specific peptide to plot (NULL for multiple)

- n_overlay:

  Number of peptides to overlay when peptide_id is NULL

## Value

A ggplot object
