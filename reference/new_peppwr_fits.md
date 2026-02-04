# Create a new peppwr_fits object

Create a new peppwr_fits object

## Usage

``` r
new_peppwr_fits(data, fits, best, call, missingness = NULL)
```

## Arguments

- data:

  Original data (nested tibble)

- fits:

  Fit results per peptide (list of tibbles with dist, loglik, aic)

- best:

  Best-fitting distribution per peptide (character vector)

- call:

  Original function call

- missingness:

  Tibble with missingness statistics per peptide (optional)

## Value

A peppwr_fits object
