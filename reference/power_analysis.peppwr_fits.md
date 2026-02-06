# Power analysis for per-peptide mode using fitted distributions

Power analysis for per-peptide mode using fitted distributions

## Usage

``` r
# S3 method for class 'peppwr_fits'
power_analysis(
  distribution,
  effect_size = NULL,
  n_per_group = NULL,
  target_power = NULL,
  alpha = 0.05,
  test = "wilcoxon",
  find = "power",
  n_sim = 1000,
  on_fit_failure = "exclude",
  proportion_threshold = 0.5,
  include_missingness = FALSE,
  apply_fdr = FALSE,
  prop_null = 0.9,
  fdr_threshold = 0.05,
  ...
)
```

## Arguments

- distribution:

  A peppwr_fits object from fit_distributions()

- effect_size:

  Fold change to detect

- n_per_group:

  Sample size per group (required for find="power")

- target_power:

  Target power (required for find="sample_size")

- alpha:

  Significance level (default 0.05)

- test:

  Statistical test to use (default "wilcoxon")

- find:

  What to solve for: "power" or "sample_size"

- n_sim:

  Number of simulations per peptide (default 1000)

- on_fit_failure:

  How to handle failed fits: "exclude", "empirical", or "lognormal"

- proportion_threshold:

  Proportion of peptides that must reach target_power (default 0.5)

- include_missingness:

  If TRUE, incorporate peptide-specific NA rates into simulations

- apply_fdr:

  If TRUE, use FDR-aware simulation with Benjamini-Hochberg correction.
  Note: not compatible with `test = "bayes_t"` (Bayes factors cannot be
  converted to p-values)

- prop_null:

  Proportion of true null peptides (default 0.9 = 90% unchanged)

- fdr_threshold:

  FDR threshold for calling discoveries (default 0.05)

- ...:

  Additional arguments (ignored)

## Value

A peppwr_power object
