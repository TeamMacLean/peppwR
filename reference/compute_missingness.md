# Compute missingness statistics for a vector of values

Calculates the number and proportion of missing (NA) values.

## Usage

``` r
compute_missingness(values)
```

## Arguments

- values:

  Numeric vector (may contain NAs)

## Value

List with:

- n_total: Total number of values

- n_missing: Number of NA values

- na_rate: Proportion missing (0-1)

## See also

[`compute_dataset_mnar()`](https://teammaclean.github.io/peppwR/reference/compute_dataset_mnar.md)
for dataset-level MNAR detection
