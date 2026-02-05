# Real-World Analysis Examples Implementation Plan (v3)

## Overview

Create two static analysis documents demonstrating peppwR on real
phosphoproteomics data. These are pre-rendered examples (not vignettes)
stored in `inst/examples/` with committed HTML output.

**Key changes from v2:** - Fixed incorrect assumptions about results -
text must match rendered figures - Wilcoxon test has essentially zero
power at N=3 with 2-fold effect - restructure to discover this early -
Use Bayes factor test as primary test after comparison reveals
Wilcoxon’s limitations - Distribution fitting shows Pareto/Skew Normal
dominance (artifact of small sample size) - explain why - Remove/rework
minimum detectable effect section (empty plots with Wilcoxon) - Improve
PRM pedagogical tone to match DDA document

------------------------------------------------------------------------

## Critical Learnings from v2 Rendering

### Distribution Fitting Results

- **Observed**: Pareto and Skew Normal are best fits, not
  lognormal/gamma
- **Reason**: With only N=3 per condition (6 observations per peptide),
  there’s insufficient data for reliable distribution selection.
  Gamma/lognormal would be expected with more data.
- **Action**: Explain this is an artifact of small sample size, not a
  statement about true underlying distributions

### Wilcoxon Test Power

- **Observed**: Essentially zero power at N=3 with 2-fold effect
- **Reason**: Non-parametric tests like Wilcoxon are conservative and
  need larger samples
- **Action**: Compare tests EARLY in the analysis to discover this, then
  switch to Bayes factor

### Minimum Detectable Effect

- **Observed**: Empty plots showing nothing detectable even at 0% power
  threshold
- **Reason**: Wilcoxon has no power, so MDE analysis is meaningless
- **Action**: Either remove this section or run it with Bayes factor
  test

### Power Analysis Plots

- **Observed**: Empty histograms showing no peptides achieve any power
- **Action**: After discovering Wilcoxon’s limitations, redo with Bayes
  factor to show meaningful results

------------------------------------------------------------------------

## Data Summary

### DDA Data (`sample_data/dda_data.csv`)

- **Source**: Arabidopsis phosphoproteomics
- **Peptides**: 2,228 unique phosphopeptides
- **Comparison**: Timepoint 0 vs 600 (earliest vs latest)
- **Replicates**: 3 biological replicates per timepoint (6 observations
  per peptide total)
- **Missing**: 0%
- **Key columns**: `new_annotation` (peptide ID), `genotype`, `bio_rep`,
  `timepoints`, `timepoints_values` (abundance)

### PRM Data (`sample_data/prm_data.csv`)

- **Source**: Fungal phosphoproteomics (rice blast)
- **Peptides**: 285 unique phosphopeptides
- **Comparison**: Timepoint 0 vs 6 (earliest vs latest)
- **Replicates**: 3 bio × 2 tech per condition
- **Missing**: 17% (good for missingness demo)
- **Key columns**: `peptide_modified_sequence` (peptide ID), `genotype`,
  `timepoint`, `bio_rep`, `tech_rep`, `total_area` (abundance)

------------------------------------------------------------------------

## Document 1: DDA Time-Course Power Analysis

**File**: `inst/examples/dda-time-course-power.Rmd`

### Revised Outline

``` markdown
# Power Analysis for Time-Course Phosphoproteomics (DDA)

## Introduction
- What is statistical power and why does it matter?
- The challenge: phosphoproteomics experiments are expensive, sample sizes are small
- Goal: Determine what effects we can detect and plan future experiments
- Overview of the peppwR approach: fit distributions → simulate experiments → estimate power

## About the Data
- Arabidopsis phosphoproteomics time-course experiment
- DDA (Data-Dependent Acquisition) mass spectrometry
- Comparing early (t=0) vs late (t=600) timepoints
- 3 biological replicates per condition (important limitation!)
- 2,228 phosphopeptides quantified

## Data Preparation
- Load and filter data
- Reshape for peppwR (peptide_id, condition, abundance)
- Quick EDA: abundance distribution across conditions

## Distribution Fitting
- Why fit distributions? (To simulate realistic experiments)
- Run fit_distributions()
- **Interpret results honestly**: Pareto and Skew Normal dominate
- **Explain why**: With only 6 observations per peptide, distribution selection is unreliable
- This is an artifact of small sample size, not the "true" underlying distributions
- With more replicates, gamma/lognormal would likely fit better
- plot_param_distribution() with sample size context

## Power Analysis

### Choosing the Right Statistical Test

**IMPORTANT: Start with test comparison before doing other analyses**

With small samples (N=3), the choice of statistical test dramatically affects power:

| Test | Type | Characteristics |
|------|------|-----------------|
| Wilcoxon rank-sum | Non-parametric | Conservative, needs larger N |
| Bootstrap-t | Resampling | Handles non-normality |
| Bayes factor | Bayesian | More powerful at small N |

Compare all three tests at N=3 with 2-fold effect:

```r
# Run all three tests
power_wilcox <- power_analysis(fits, effect_size = 2, n_per_group = 3, test = "wilcoxon")
power_boot <- power_analysis(fits, effect_size = 2, n_per_group = 3, test = "bootstrap_t")
power_bayes <- power_analysis(fits, effect_size = 2, n_per_group = 3, test = "bayes_t")
```

Present results as a clear table:

| Statistical Test  | Median Power | % Peptides \> 50% Power | % Peptides \> 80% Power |
|-------------------|--------------|-------------------------|-------------------------|
| Wilcoxon rank-sum | X.XX         | X%                      | X%                      |
| Bootstrap-t       | X.XX         | X%                      | X%                      |
| Bayes factor      | X.XX         | X%                      | X%                      |

**Key finding**: Wilcoxon has essentially zero power at N=3. This is
expected - non-parametric tests are conservative and need larger
samples. The Bayes factor test provides usable power estimates.

**For remaining analyses, use Bayes factor test.**

### The Three Questions (Using Bayes Factor)

Power analysis can answer three related questions: 1. Given sample size
and effect size → what is my power? 2. Given target power and effect
size → what sample size do I need? 3. Given sample size and target power
→ what’s the minimum detectable effect?

### Question 1: Current Power (N=3, Bayes Factor)

- With 3 replicates, what power do we have to detect a 2-fold change?
- Show actual distribution of power across peptides
- Interpret: “X% of peptides achieve \>80% power”

### Question 2: Sample Size for Target Power

- What N do we need for 80% power at 2-fold change? (use realistic
  effect)
- Understanding the power curve
- Show % peptides reaching target at each N

### Question 3: Minimum Detectable Effect

**KNOWN ISSUE:** `find = "effect_size"` is NOT YET IMPLEMENTED for
per-peptide mode. See CLAUDE.md “Known Issues” section for details and
fix plan.

**Current workaround options:** 1. **Remove this section** from
per-peptide examples until implemented 2. **Use aggregate mode** with
representative parameters derived from fits: \`\`\`r \# Derive
representative lognormal params from fitted data \# (most peptides are
approximately lognormal on log scale) log_abundances \<-
log(pilot$abundance\lbrack pilot$abundance \> 0\]) representative_params
\<- list( meanlog = median(log_abundances), sdlog = sd(log_abundances) )

\# Run MDE in aggregate mode min_effect \<- power_analysis( distribution
= “lnorm”, params = representative_params, n_per_group = 3, target_power
= 0.8, find = “effect_size”, test = “bayes_t”, n_sim = 100 ) \`\`\` 3.
**Wait for v2.2** implementation of per-peptide MDE

**Recommendation:** For now, remove Question 3 from per-peptide examples
OR use aggregate mode with clear explanation that this represents a
“typical” peptide, not the full peptidome distribution.

## Power Heatmap

### Why Aggregate Mode for Heatmaps?

The
[`plot_power_heatmap()`](https://teammaclean.github.io/peppwR/reference/plot_power_heatmap.md)
function uses **aggregate mode** (single distribution) rather than
per-peptide mode because:

1.  **Computational cost**: Per-peptide heatmaps would require
    `n_peptides × n_values × effect_values × n_sim` simulations (e.g.,
    2228 × 10 × 10 × 100 = 22M simulations)
2.  **Visualization**: A single heatmap is more interpretable than 2000+
    individual heatmaps
3.  **Purpose**: Heatmaps answer “what’s the general tradeoff?” not
    “what’s the power for peptide X?”

### Deriving Representative Parameters

To make the heatmap relevant to your data, derive parameters from the
fitted distributions:

``` r
# Option 1: Use median log-abundance parameters
log_abundances <- log(pilot$abundance[pilot$abundance > 0])
params <- list(
  meanlog = median(log_abundances),
  sdlog = mad(log_abundances, constant = 1)  # robust SD estimate
)

# Option 2: For gamma-like data, use moment matching
abundances <- pilot$abundance[pilot$abundance > 0]
m <- mean(abundances)
v <- var(abundances)
params <- list(
  shape = m^2 / v,
  rate = m / v
)
```

### Heatmap Code with Explanation

``` r
# Derive representative parameters from the data
log_abundances <- log(pilot$abundance[pilot$abundance > 0])
derived_params <- list(
  meanlog = median(log_abundances),
  sdlog = mad(log_abundances, constant = 1)
)

cat("Using representative lognormal parameters:\n")
cat("  meanlog =", round(derived_params$meanlog, 2), "\n")
cat("  sdlog =", round(derived_params$sdlog, 2), "\n")

# Generate heatmap
plot_power_heatmap(
  distribution = "lnorm",
  params = derived_params,
  n_range = c(3, 12),
  effect_range = c(1.5, 4),
  test = "wilcoxon"  # heatmap uses aggregate mode
)
```

### Interpreting the Heatmap

The heatmap shows power for a **representative peptide** with typical
abundance characteristics. Individual peptides will vary: -
Low-variability peptides will have higher power than shown -
High-variability peptides will have lower power than shown - Use the
per-peptide `find = "power"` results to understand the distribution
across peptides

## Summary and Recommendations

- Key findings:
  - Distribution fitting is unreliable at N=3 (Pareto/Skew Normal
    artifact)
  - Wilcoxon test has no power at N=3 - use Bayes factor instead
  - Actual power estimates with Bayes factor
- Practical recommendations:
  - Minimum sample size for reliable detection
  - Effect sizes that are detectable
- Caveats and limitations

## Session Info

    ### Key Code Changes

    ```r
    # Test comparison - do this FIRST
    power_wilcox <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                                    find = "power", test = "wilcoxon", n_sim = 100)
    power_boot <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                                  find = "power", test = "bootstrap_t", n_sim = 100)
    power_bayes <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                                   find = "power", test = "bayes_t", n_sim = 100)

    # Create comparison table
    comparison <- tibble::tibble(
      Test = c("Wilcoxon rank-sum", "Bootstrap-t", "Bayes factor"),
      `Median Power` = c(
        median(power_wilcox$simulations$peptide_power, na.rm = TRUE),
        median(power_boot$simulations$peptide_power, na.rm = TRUE),
        median(power_bayes$simulations$peptide_power, na.rm = TRUE)
      ),
      `% > 50% Power` = c(
        mean(power_wilcox$simulations$peptide_power > 0.5, na.rm = TRUE) * 100,
        mean(power_boot$simulations$peptide_power > 0.5, na.rm = TRUE) * 100,
        mean(power_bayes$simulations$peptide_power > 0.5, na.rm = TRUE) * 100
      ),
      `% > 80% Power` = c(
        mean(power_wilcox$simulations$peptide_power > 0.8, na.rm = TRUE) * 100,
        mean(power_boot$simulations$peptide_power > 0.8, na.rm = TRUE) * 100,
        mean(power_bayes$simulations$peptide_power > 0.8, na.rm = TRUE) * 100
      )
    )
    knitr::kable(comparison, digits = 2)

    # Subsequent analyses use Bayes factor
    power_n3 <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                               find = "power", test = "bayes_t", n_sim = 100)

    sample_size <- power_analysis(fits, effect_size = 2, target_power = 0.8,
                                  find = "sample_size", test = "bayes_t", n_sim = 100)

    # Power heatmap - derive representative parameters from data
    log_abundances <- log(pilot$abundance[pilot$abundance > 0])
    derived_params <- list(
      meanlog = median(log_abundances),
      sdlog = mad(log_abundances, constant = 1)
    )

    plot_power_heatmap(
      distribution = "lnorm",
      params = derived_params,
      n_range = c(3, 12),
      effect_range = c(1.5, 4),
      test = "wilcoxon"  # heatmap uses aggregate mode
    )

------------------------------------------------------------------------

## Document 2: PRM Analysis with Missing Data

**File**: `inst/examples/prm-genotype-power.Rmd`

### Revised Outline (Improved Pedagogy)

``` markdown
# Power Analysis for Targeted Proteomics with Missing Data (PRM)

## Introduction

### What Makes Targeted Proteomics Different?

Mass spectrometry-based proteomics operates in two fundamentally different modes:

**Discovery proteomics (DDA)** casts a wide net:
- Measures thousands of proteins/peptides
- Identifies whatever is abundant enough to detect
- Great for hypothesis generation
- High multiple testing burden

**Targeted proteomics (PRM/SRM)** takes aim at specific targets:
- Measures a pre-defined panel of peptides
- Higher sensitivity and reproducibility
- Better for hypothesis testing
- Lower multiple testing burden

This document analyzes a PRM dataset, demonstrating how peppwR handles the unique challenges of targeted proteomics.

### The Missing Data Challenge

Even with targeted methods, missing values are ubiquitous in mass spectrometry. But not all missing data is created equal:

**MCAR (Missing Completely At Random)**
- Values missing due to random technical failures
- No relationship to abundance
- Reduces sample size but doesn't bias results

**MNAR (Missing Not At Random)**
- Low-abundance peptides systematically missing
- Below detection limit = no measurement
- Can bias abundance estimates upward
- Common in mass spectrometry

peppwR tracks both the rate of missingness AND evidence for MNAR patterns, enabling more realistic power estimates.

### Goals of This Analysis

1. Characterize missingness patterns in the data
2. Identify peptides with MNAR evidence
3. Compare power with/without accounting for missingness
4. Determine realistic sample size requirements
5. Understand FDR impact on power

## About the Data

This analysis uses fungal phosphoproteomics data from a targeted PRM experiment:

| Property | Value |
|----------|-------|
| Organism | *Magnaporthe oryzae* (rice blast fungus) |
| Method | PRM (Parallel Reaction Monitoring) |
| Comparison | Early (t=0) vs late (t=6) timepoints |
| Biological replicates | 3 per condition |
| Technical replicates | 2 per biological replicate |
| Peptides | 285 phosphopeptides |
| Missing data | ~17% |

## Data Preparation

[Same as before but with better explanatory text]

### Why Average Technical Replicates?

Technical replicates measure the same biological sample multiple times. They capture *measurement* variability, not *biological* variability.

For power analysis, biological replicates are the true unit of replication - they represent independent observations of the biological phenomenon we're studying. We average technical replicates to get one value per biological replicate.

## Missingness Analysis

### Examining the Extent of Missingness

Before diving into power analysis, we need to understand our missing data:

```r
fits <- fit_distributions(pilot, "peptide_id", "condition", "abundance")
print(fits)  # Shows missingness summary
```

### Visualizing Missingness Patterns

The
[`plot_missingness()`](https://teammaclean.github.io/peppwR/reference/plot_missingness.md)
function provides three complementary views:

1.  **NA Rate Distribution**: What fraction of observations are missing
    for each peptide?
2.  **MNAR Score Distribution**: Is there evidence that missingness
    correlates with abundance?
3.  **Abundance vs NA Rate**: Do low-abundance peptides have more
    missing values?

``` r
plot_missingness(fits)
```

### Interpreting MNAR Scores

The MNAR score is a z-statistic testing whether missing observations
have systematically different (lower) abundance than observed values.

| MNAR Score | Interpretation           |
|------------|--------------------------|
| \< 1       | Little evidence for MNAR |
| 1-2        | Weak evidence            |
| 2-3        | Moderate evidence        |
| \> 3       | Strong evidence for MNAR |

``` r
mnar_peptides <- get_mnar_peptides(fits, threshold = 2)
```

### What to Do with MNAR Peptides

Peptides with high MNAR scores require careful handling: - Their
abundance estimates may be biased upward - Power calculations may be
optimistic - Consider robust statistical methods - Report separately in
publications

## Distribution Fitting

### Why Fit Distributions?

\[Same pedagogical content as DDA document\]

### Results: A Cautionary Note on Small Samples

``` r
plot(fits)
```

**Important**: Like the DDA dataset, we see Pareto and Skew Normal as
dominant best-fit distributions. This is an **artifact of small sample
size** (only 6 observations per peptide when averaged across tech reps),
not a statement about true underlying distributions.

With more biological replicates, we would expect gamma and lognormal
distributions to fit better - these are the typical distributions for
mass spectrometry abundance data.

## Power Analysis

### Choosing the Right Statistical Test

**Critical first step**: Compare statistical tests before committing to
one.

With small samples (N=3), test choice dramatically affects power:

``` r
power_wilcox <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                                test = "wilcoxon", n_sim = 100)
power_boot <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                              test = "bootstrap_t", n_sim = 100)
power_bayes <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                               test = "bayes_t", n_sim = 100)
```

\[Present as table like DDA document\]

**Key finding**: \[Describe actual results - likely Wilcoxon has no
power, Bayes factor is best\]

**For remaining analyses, use Bayes factor test.**

### What Effect Size Can We Detect?

Before asking “what’s my power at 2-fold?”, ask “what effect size IS
detectable?”

``` r
# First, understand what's detectable
min_effect <- power_analysis(fits, n_per_group = 3, target_power = 0.8,
                             find = "effect_size", test = "bayes_t", n_sim = 100)
```

This informs whether a 2-fold threshold is realistic or if we need to
consider larger effects.

### Power at Realistic Effect Sizes

\[Use effect size informed by min_effect analysis\]

### Impact of Missingness on Power

Compare power estimates with and without accounting for missingness:

``` r
power_no_miss <- power_analysis(fits, effect_size = X, n_per_group = 3,
                                 test = "bayes_t", include_missingness = FALSE)
power_with_miss <- power_analysis(fits, effect_size = X, n_per_group = 3,
                                   test = "bayes_t", include_missingness = TRUE)
```

Missingness reduces effective sample size. The difference between these
estimates shows the “cost” of missing data.

### FDR-Aware Power

With 285 peptides, multiple testing correction matters:

``` r
power_fdr <- power_analysis(fits, effect_size = X, n_per_group = 3,
                            test = "bayes_t", apply_fdr = TRUE,
                            prop_null = 0.8, n_sim = 100)
```

**Understanding `prop_null`**: This parameter specifies the assumed
proportion of true null hypotheses. With `prop_null = 0.8`, we assume
80% of peptides have no true effect. This affects how stringent the FDR
correction is.

## Summary and Recommendations

### Key Findings

1.  **Missingness**: ~17% missing with \[X\] peptides showing MNAR
    evidence
2.  **Test selection**: Wilcoxon has no power at N=3; Bayes factor is
    usable
3.  **Distribution fitting**: Pareto/Skew Normal dominance is artifact
    of small N
4.  **Power estimates**: \[Actual results\]
5.  **FDR impact**: \[Actual results\]

### Practical Recommendations

\[Based on actual results\]

### Caveats

- Small sample size limits distribution fitting reliability
- MNAR models are approximations
- prop_null requires assumptions about true effect rates

## Session Info

\`\`\`

------------------------------------------------------------------------

## Implementation Steps

1.  Update DDA Rmd following revised outline

    - Move test comparison to EARLY in analysis
    - Fix distribution fitting interpretation
    - Use Bayes factor for main analyses
    - Create proper comparison table
    - Remove/fix minimum detectable effect section

2.  Render and **carefully verify** DDA HTML

    - Check each figure matches its accompanying text
    - Verify table values are populated correctly

3.  Update PRM Rmd following revised outline

    - Improve pedagogical tone throughout
    - Same test comparison approach as DDA
    - Better MNAR explanation
    - Use effect size informed by data, not assumed 2-fold

4.  Render and **carefully verify** PRM HTML

    - Same verification as DDA

5.  Run devtools::check(vignettes = FALSE)

6.  Commit changes

------------------------------------------------------------------------

## Writing Guidelines

### Text Must Match Figures

**NEVER** write explanatory text that assumes a result. Instead: - Run
the analysis first - Look at the actual output/figure - Write text that
accurately describes what you see

Bad: “The results show lognormal distributions dominate…” Good: “The
results show \[describe actual bars in figure\]…”

### Explain Unexpected Results

When results don’t match intuition, explain why: - Pareto/Skew Normal
dominance → small sample size artifact - Zero Wilcoxon power →
non-parametric tests need larger N - Empty MDE plots → using wrong test
for the data

### Use Tables for Comparisons

When comparing multiple options (tests, conditions), use a clear table:

| Option A | Option B | Option C |
|----------|----------|----------|
| value    | value    | value    |

### Pedagogical Consistency

Both documents should have similar tone and structure: - Explain “why”
before “how” - Use consistent terminology - Similar section headers
where appropriate

------------------------------------------------------------------------

## Verification Checklist

After implementation: - \[ \] DDA Rmd renders without error - \[ \] DDA
figures match their explanatory text - \[ \] Test comparison table shows
actual computed values - \[ \] Bayes factor used for main analyses (not
Wilcoxon) - \[ \] Distribution fitting text explains Pareto/Skew Normal
artifact - \[ \] PRM Rmd renders without error - \[ \] PRM figures match
their explanatory text - \[ \] PRM has same pedagogical quality as DDA -
\[ \] Effect size threshold informed by data (not assumed) - \[ \]
devtools::check(vignettes = FALSE) passes with 0 errors, 0 warnings - \[
\] Changes committed
