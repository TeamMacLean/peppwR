# Compute dataset-level MNAR metric

Calculates Spearman correlation between log(mean_abundance) and NA rate
across peptides to detect whether low-abundance peptides have more
missing values than high-abundance peptides.

## Usage

``` r
compute_dataset_mnar(missingness)
```

## Arguments

- missingness:

  A tibble with columns na_rate and mean_abundance (as produced by
  fit_distributions)

## Value

A list with:

- correlation: Spearman correlation coefficient (negative = MNAR
  pattern)

- p_value: p-value for the correlation

- n_peptides: Number of peptides with missing data used in calculation

- interpretation: Human-readable interpretation string

## MNAR Detection

MNAR (Missing Not At Random) in mass spectrometry typically manifests as
low-abundance peptides having higher rates of missing values due to
detection limits. This function detects this pattern by correlating mean
abundance with NA rate across all peptides.

A negative correlation indicates that low-abundance peptides have more
missing values - the hallmark of detection-limit-driven MNAR.

## Interpretation

|                   |                         |
|-------------------|-------------------------|
| Correlation (r)   | Interpretation          |
| r \> -0.1         | No evidence of MNAR     |
| -0.3 \< r \< -0.1 | Weak evidence           |
| -0.5 \< r \< -0.3 | Moderate evidence       |
| r \< -0.5         | Strong evidence of MNAR |
