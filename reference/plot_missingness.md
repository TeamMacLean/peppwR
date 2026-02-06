# Plot missingness statistics

Creates a visualization of missing data patterns with two panels:

## Usage

``` r
plot_missingness(fits)
```

## Arguments

- fits:

  A peppwr_fits object

## Value

A ggplot object or gtable (combined panels)

## Details

- **Panel 1**: Distribution of NA rates across peptides

- **Panel 2**: Mean abundance vs NA rate scatter plot showing the
  **dataset-level MNAR correlation**. The subtitle displays the Spearman
  correlation coefficient and p-value. A negative correlation indicates
  that low-abundance peptides have more missing values - the hallmark of
  detection-limit-driven MNAR.

## MNAR Detection

MNAR (Missing Not At Random) in mass spectrometry typically manifests as
low-abundance peptides having higher rates of missing values due to
detection limits. Panel 2 visualizes this relationship and reports the
correlation coefficient.

## See also

[`compute_dataset_mnar()`](https://teammaclean.github.io/peppwR/reference/compute_dataset_mnar.md)
for the underlying correlation calculation
