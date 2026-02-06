# Fit distributions to peptide abundance data

Fits candidate distributions to each peptide's abundance values and
selects the best fit by AIC. Also computes missingness statistics
including dataset-level MNAR detection.

## Usage

``` r
fit_distributions(data, id, group, value, distributions = "continuous")
```

## Arguments

- data:

  A data frame containing peptide abundance data

- id:

  Column name for peptide identifier

- group:

  Column name for group/condition

- value:

  Column name for abundance values

- distributions:

  Which distributions to fit: "continuous" (default), "counts", "all",
  or a character vector of distribution names

## Value

A peppwr_fits object containing:

- `$data`: Nested tibble with original data and fit results

- `$best`: Best-fitting distribution for each peptide

- `$missingness`: Per-peptide missingness statistics

- `$dataset_mnar`: Dataset-level MNAR correlation and interpretation

## Missingness Tracking

The returned object includes:

- Per-peptide NA rates (in `$missingness`)

- Dataset-level MNAR correlation (in `$dataset_mnar`)

The dataset-level MNAR metric correlates log(mean_abundance) with NA
rate across peptides. A negative correlation indicates low-abundance
peptides have more missing values - typical of detection-limit-driven
MNAR.

Print the result to see both metrics. For small sample sizes (N \< 15),
the dataset-level correlation is more reliable than per-peptide scores.

## See also

[`compute_dataset_mnar()`](https://teammaclean.github.io/peppwR/reference/compute_dataset_mnar.md)
for dataset-level MNAR details,
[`plot_missingness()`](https://teammaclean.github.io/peppwR/reference/plot_missingness.md)
to visualize missingness patterns
