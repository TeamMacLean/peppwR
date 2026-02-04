# Fit distributions to peptide abundance data

Fit distributions to peptide abundance data

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

A peppwr_fits object
