# Compute missingness statistics for a vector of values

Calculates NA rate and MNAR (Missing Not At Random) score. MNAR
detection uses the observation that under MCAR, the mean rank of
observed values should be (n+1)/2. Under MNAR with low values missing,
the mean rank will be higher.

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

- mnar_score: Z-score measuring MNAR pattern. Positive values indicate
  low values are more likely to be missing. Values \> 2 suggest MNAR.
