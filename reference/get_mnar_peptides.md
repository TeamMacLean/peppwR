# Get peptides showing MNAR pattern

Returns peptides where the MNAR score exceeds a threshold, indicating
that low values are systematically missing (Missing Not At Random).

## Usage

``` r
get_mnar_peptides(fits, threshold = 2)
```

## Arguments

- fits:

  A peppwr_fits object

- threshold:

  MNAR score threshold (default 2, corresponds to ~95% confidence)

## Value

A tibble with columns:

- peptide_id: Peptide identifier

- condition: Group/condition

- na_rate: Proportion of missing values

- mnar_score: Z-score for MNAR pattern

- mean_abundance: Mean of observed values Sorted by mnar_score
  descending. Returns empty tibble if no peptides exceed threshold.
