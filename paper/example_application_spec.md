# Example Application Rmd Specification

Supplementary material for Bioinformatics Application Note.

**File**: `paper/supplementary_example.Rmd`

**Purpose**: Fully reproducible code generating all results and Figure 1 panels B-D.

---

## YAML Header

```yaml
title: "peppwR Example Application: Supplementary Material"
author: "Dan MacLean"
output:
  html_document:
    toc: true
    toc_float: true
    code_folding: show
  pdf_document:
    toc: true
---
```

---

## Section 1: Setup

```r
library(peppwR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

set.seed(42)  # Reproducibility
```

**Output**: None (setup only)

---

## Section 2: Simulate Phosphoproteomics Pilot Data

### 2.1 Define peptide parameters

- 500 peptides
- Shape: uniform(1.5, 5) - captures low-to-high abundance peptides
- Rate: uniform(0.05, 0.2) - captures variance heterogeneity
- Results in mean abundances ranging ~7.5 to 100 (realistic intensity range after log-transform)

```r
n_peptides <- 500
n_pilot <- 6  # samples per condition in pilot

peptide_params <- tibble(

  peptide = paste0("pep_", sprintf("%04d", 1:n_peptides)),
  shape = runif(n_peptides, min = 1.5, max = 5),
  rate = runif(n_peptides, min = 0.05, max = 0.2)
) |>
  mutate(
    theoretical_mean = shape / rate,
    theoretical_cv = 1 / sqrt(shape)  # CV for gamma

)
```

**Output**: Table showing parameter ranges and summary statistics

### 2.2 Generate abundance values

```r
pilot_data <- peptide_params |>
  rowwise() |>
  mutate(
    abundance = list(rgamma(n_pilot, shape = shape, rate = rate))
  ) |>
  ungroup() |>
  select(peptide, abundance) |>
  unnest(abundance) |>
  mutate(
    sample = rep(paste0("S", 1:n_pilot), n_peptides),
    condition = "control"
  )
```

**Output**: Dimensions of pilot_data, head()

### 2.3 Introduce MNAR missingness

- Target ~15% overall missingness
- Probability of NA inversely related to abundance
- Use logistic model: P(NA) = 1 / (1 + exp(k * (abundance - threshold)))

```r
# Calculate per-peptide mean abundance
peptide_means <- pilot_data |>
  group_by(peptide) |>
  summarise(mean_abundance = mean(abundance))

# Set missingness threshold at 25th percentile of means
miss_threshold <- quantile(peptide_means$mean_abundance, 0.25)
k <- 0.15  # steepness of logistic curve

pilot_data <- pilot_data |>
  left_join(peptide_means, by = "peptide") |>
  mutate(
    p_missing = 1 / (1 + exp(k * (mean_abundance - miss_threshold))),
    is_missing = runif(n()) < p_missing,
    abundance = if_else(is_missing, NA_real_, abundance)
  ) |>
  select(peptide, sample, condition, abundance)
```

**Output**:
- Overall missingness rate (target ~15%)
- Correlation between mean abundance and missingness rate
- Plot: missingness rate vs mean abundance (showing MNAR pattern)

---

## Section 3: Distribution Fitting

```r
fits <- fit_distributions(
  pilot_data,
  id = "peptide",
  group = "condition",
  value = "abundance",
  distributions = "continuous"
)

print(fits)
summary(fits)
```

**Output**:
- Fit summary (distribution counts, failure rate)
- Missingness diagnostics from fits object

---

## Section 4: Primary Power Analysis

### 4.1 Find required sample size

```r
result_primary <- power_analysis(

  fits,
  effect_size = 2,
  target_power = 0.8,
  proportion_threshold = 0.8,  # 80% of peptides

  find = "sample_size",
  test = "wilcoxon",
  n_sim = 1000
)

print(result_primary)
```

**Output**: "N = X per group required for 80% of peptides to achieve 80% power to detect 2-fold changes"

**Expected**: N ~ 8

### 4.2 Extract detailed results for plotting

```r
# Power at each sample size for all peptides
power_by_n <- result_primary$simulations
```

---

## Section 5: Comparison 1 - Per-Peptide vs Homogeneous

### 5.1 Homogeneous (aggregate) approach

Calculate pooled parameters:

```r
# Pool all non-NA values
pooled_data <- pilot_data |>
  filter(!is.na(abundance))

# Fit single gamma to pooled data
pooled_fit <- fitdistrplus::fitdist(pooled_data$abundance, "gamma")
pooled_params <- as.list(pooled_fit$estimate)

# Power analysis with homogeneous assumption
result_homogeneous <- power_analysis(
  "gamma",
  params = pooled_params,
  effect_size = 2,
  target_power = 0.8,
  find = "sample_size",
  test = "wilcoxon",
  n_sim = 1000
)

print(result_homogeneous)
```

**Output**: "N = X per group" (expected ~5)

### 5.2 Comparison summary

```r
comparison_1 <- tibble(

  approach = c("Homogeneous (pooled)", "Per-peptide (peppwR)"),
  required_n = c(result_homogeneous$answer, result_primary$answer),
  description = c(
    "Single distribution fitted to pooled data",
    "Distribution fitted per peptide, 80% threshold"
  )
)

knitr::kable(comparison_1, caption = "Sample size estimates by approach")
```

**Output**: Table showing per-peptide requires ~60% more samples than homogeneous

---

## Section 5.5: Comparison 1b - With vs Without MNAR

peppwR can incorporate observed missingness patterns into simulations.

```r
# With MNAR-aware simulation
result_with_mnar <- power_analysis(
  fits,
  effect_size = 2,
  target_power = 0.8,
  proportion_threshold = 0.8,
  find = "sample_size",
  test = "wilcoxon",
  n_sim = 1000,
  include_missingness = TRUE  # <-- key parameter
)

print(result_with_mnar)
```

**Expected finding:**
- Without MNAR: N = 8
- With MNAR: N = 9-10 (missingness reduces effective sample size)

**Message**: Accounting for realistic missingness patterns shows that more samples are needed to maintain power when some observations will be missing.

---

## Section 6: Comparison 2 - Statistical Test Choice

### 6.1 Bayes factor t-test

```r
result_bayes <- power_analysis(
  fits,
  effect_size = 2,
  target_power = 0.8,
  proportion_threshold = 0.8,
  find = "sample_size",
  test = "bayes_t",
  n_sim = 1000
)

print(result_bayes)
```

**Output**: "N = X per group" (expected ~6-7)

### 6.2 Comparison summary

```r
comparison_2 <- tibble(
  test = c("Wilcoxon rank-sum", "Bayes factor t-test"),
  required_n = c(result_primary$answer, result_bayes$answer),
  characteristics = c(
    "Non-parametric, conservative, no assumptions",
    "Parametric Bayesian, more efficient when assumptions hold"
  )
)

knitr::kable(comparison_2, caption = "Sample size estimates by statistical test")
```

**Output**: Table showing Bayes ~15-20% more efficient

---

## Section 7: Figure 1 Generation

### 7.1 Panel A: Workflow diagram (conceptual)

Not generated in R - create separately as vector graphic or use DiagrammeR:

```r
# Optional: simple DiagrammeR version
library(DiagrammeR)

grViz("
  digraph workflow {
    rankdir=LR
    node [shape=box, style=rounded]

    pilot [label='Pilot Data']
    fit [label='fit_distributions()']
    power [label='power_analysis()']
    result [label='Result']

    pilot -> fit -> power -> result
  }
")
```

**Note**: For publication, create clean vector version in Inkscape/Illustrator.

### 7.2 Panel B: Power heterogeneity histogram

```r
# NOTE: find="sample_size" doesn't store per-peptide power
# Need separate call with find="power" at fixed N=6
result_n6 <- power_analysis(
  fits,
  effect_size = 2,
  n_per_group = 6,
  find = "power",
  test = "wilcoxon",
  n_sim = 1000
)

# Extract per-peptide power vector
power_at_n6 <- data.frame(power = result_n6$simulations$peptide_power)

panel_b <- ggplot(power_at_n6, aes(x = power * 100)) +

geom_histogram(
    binwidth = 5,
    fill = "#FEB24C",  # YlOrRd mid-tone
    color = "white",
    boundary = 0
  ) +
  geom_vline(xintercept = 80, linetype = "dashed", color = "#BD0026", linewidth = 0.8) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    x = "Power (%)",
    y = "Number of peptides",
    title = "B) Per-peptide power at N=6"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank()
  )
```

**Output**: Histogram showing heterogeneity; bimodal expected

### 7.3 Panel C: Power curve with test comparison

```r
# Extract power curves from find="sample_size" results
# These contain proportion_powered at each N

curve_data <- bind_rows(
  result_primary$simulations$power_curve |>
    mutate(pct_powered = proportion_powered * 100, test = "Wilcoxon"),
  result_bayes$simulations$power_curve |>
    mutate(pct_powered = proportion_powered * 100, test = "Bayes factor")
)

panel_c <- ggplot(curve_data, aes(x = n_per_group, y = pct_powered,
                                   color = test, linetype = test)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 80, linetype = "dotted", color = "gray40") +
  scale_color_manual(values = c("Wilcoxon" = "#BD0026", "Bayes factor" = "#31A354")) +
  scale_linetype_manual(values = c("Wilcoxon" = "solid", "Bayes factor" = "dashed")) +
  scale_x_continuous(breaks = seq(3, 15, 2)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    x = "Sample size per group",
    y = "% peptides at 80% power",
    title = "C) Power curve by statistical test",
    color = "Test",
    linetype = "Test"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = c(0.75, 0.25),
    panel.grid.minor = element_blank()
  )
```

**Output**: Two curves showing Bayes reaches threshold at lower N

### 7.4 Panel D: Power heatmap

```r
# Generate heatmap data across effect sizes
effect_sizes <- c(1.5, 2, 2.5, 3)
sample_sizes <- 4:12

# Run power analysis for each combination (computationally intensive)
heatmap_data <- expand_grid(
  effect_size = effect_sizes,
  n_per_group = sample_sizes
) |>
  rowwise() |>
  mutate(
    result = list(
      power_analysis(fits, effect_size = effect_size, n_per_group = n_per_group,
                     find = "power", test = "wilcoxon", n_sim = 500)
    ),
    pct_powered = mean(result$simulations$power >= 0.8) * 100
  ) |>
  ungroup() |>
  select(effect_size, n_per_group, pct_powered)

panel_d <- ggplot(heatmap_data, aes(x = factor(effect_size), y = factor(n_per_group),
                                     fill = pct_powered)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(pct_powered)), size = 3, color = "black") +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, limits = c(0, 100)) +
  labs(
    x = "Effect size (fold-change)",
    y = "Sample size per group",
    title = "D) Power lookup table",
    fill = "% powered"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  )
```

**Output**: Heatmap showing trade-off between N and effect size

### 7.5 Combine into Figure 1

```r
# Panel A placeholder (or import from file)
panel_a <- ggdraw() +
  draw_label("Panel A: Workflow Diagram\n(see vector graphic)",
             size = 10, color = "gray50")

figure_1 <- plot_grid(
panel_a, panel_b,
  panel_c, panel_d,
  nrow = 2,
  labels = c("", "", "", ""),  # labels already in titles
  align = "hv"
)

ggsave("figure_1.png", figure_1, width = 8, height = 7, dpi = 300)
ggsave("figure_1.pdf", figure_1, width = 8, height = 7)
```

**Output**: Combined figure for publication

---

## Section 8: Session Info

```r
sessionInfo()
```

---

## Expected Results Summary

| Analysis | Result |
|----------|--------|
| Per-peptide, Wilcoxon | N = 8 |
| Per-peptide, Wilcoxon + MNAR | N = 9-10 |
| Homogeneous, Wilcoxon | N = 5 |
| Per-peptide, Bayes | N = 6-7 |
| More samples needed (per-peptide vs homogeneous) | ~60% |
| Efficiency gain (Bayes vs Wilcoxon) | ~15-20% |

---

## Computational Notes

- `n_sim = 1000` for primary analyses (publication quality)
- `n_sim = 500` for heatmap grid (36 combinations)
- Expected runtime: ~5-10 minutes on modern laptop
- Set `n_sim = 100` for draft/testing

---

## File Outputs

| File | Description |
|------|-------------|
| `figure_1.png` | Combined figure (300 dpi) |
| `figure_1.pdf` | Combined figure (vector) |
| `supplementary_example.html` | Rendered supplement |

