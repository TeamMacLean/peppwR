# Real-World Analysis Examples Implementation Plan

## Overview

Create two static analysis documents demonstrating peppwR on real phosphoproteomics data. These are pre-rendered examples (not vignettes) stored in `inst/examples/` with committed HTML output.

## Data Summary

### DDA Data (`sample_data/dda_data.csv`)
- **Source**: Arabidopsis phosphoproteomics
- **Peptides**: 2,228 unique phosphopeptides
- **Comparison**: Timepoint 0 vs 600 (earliest vs latest)
- **Replicates**: 3 biological replicates per timepoint
- **Missing**: 0%
- **Key columns**: `new_annotation` (peptide ID), `genotype`, `bio_rep`, `timepoints`, `timepoints_values` (abundance)

### PRM Data (`sample_data/prm_data.csv`)
- **Source**: Fungal phosphoproteomics (rice blast)
- **Peptides**: 285 unique phosphopeptides
- **Comparison**: Timepoint 0 vs 6 (earliest vs latest), both genotypes combined OR Guy11 only
- **Replicates**: 3 bio × 2 tech per condition
- **Missing**: 17% (good for missingness demo)
- **Key columns**: `peptide_modified_sequence` (peptide ID), `genotype`, `timepoint`, `bio_rep`, `tech_rep`, `total_area` (abundance)

---

## File Structure

```
inst/
  examples/
    dda-time-course-power.Rmd
    dda-time-course-power.html  (pre-rendered)
    prm-genotype-power.Rmd
    prm-genotype-power.html     (pre-rendered)
```

---

## Document 1: DDA Time-Course Power Analysis

**File**: `inst/examples/dda-time-course-power.Rmd`

**Title**: "Power Analysis for Time-Course Phosphoproteomics (DDA)"

### Outline

```markdown
# Power Analysis for Time-Course Phosphoproteomics

## Introduction
- Brief context: Arabidopsis phosphoproteomics time course
- Question: What power do we have to detect changes between early and late timepoints?

## Data Preparation
- Load data from CSV
- Filter to timepoints 0 and 600
- Reshape for peppwR: peptide_id, condition (timepoint), abundance
- Quick EDA: abundance distribution, peptide count

## Distribution Fitting
- fit_distributions() on pilot data
- plot(fits) - show best distribution counts
- Discuss: which distributions fit phosphoproteomics data?

## Diagnostic Plots
- plot_density_overlay() - sample of peptides
- plot_qq() - goodness of fit check
- plot_param_distribution() - AIC across peptidome

## Power Analysis
- Question 1: With N=3, what's our power to detect 2-fold change?
- Question 2: What sample size needed for 80% power at 1.5x effect?
- Question 3: What's the minimum detectable effect at N=3?

## Visualizations
- plot(power_result) - power curves
- plot_power_heatmap() - N vs effect size grid

## Recommendations
- Summary of findings
- Sample size recommendations for future experiments
- Caveats (single genotype, specific timepoints)
```

### Key Code Blocks

```r
# Data prep
dda <- read.csv("../../sample_data/dda_data.csv")
pilot <- dda |>
  filter(timepoints %in% c(0, 600)) |>
  transmute(
    peptide_id = new_annotation,
    condition = paste0("t", timepoints),
    abundance = timepoints_values
  )

# Fitting
fits <- fit_distributions(pilot, "peptide_id", "condition", "abundance")

# Power analyses
power_n3 <- power_analysis(fits, effect_size = 2, n_per_group = 3, find = "power")
sample_size <- power_analysis(fits, effect_size = 1.5, target_power = 0.8, find = "sample_size")
min_effect <- power_analysis(fits, n_per_group = 3, target_power = 0.8, find = "effect_size")
```

---

## Document 2: PRM Genotype Comparison with Missing Data

**File**: `inst/examples/prm-genotype-power.Rmd`

**Title**: "Power Analysis for Targeted Proteomics with Missing Data (PRM)"

### Outline

```markdown
# Power Analysis for Targeted Proteomics with Missing Data

## Introduction
- Context: PRM panel for fungal phosphoproteomics
- Comparison: Early (t=0) vs late (t=6) timepoints
- Highlight: Real-world missing data patterns

## Data Preparation
- Load data from CSV
- Filter to timepoints 0 and 6
- Handle technical replicates: average or keep separate?
- Reshape for peppwR

## Missingness Analysis (Key Section)
- Examine NA patterns in the data
- fit_distributions() with missingness tracking
- print(fits) - show missingness summary
- plot_missingness() - 3-panel visualization
- get_mnar_peptides() - identify problematic peptides
- Discussion: MNAR in targeted proteomics

## Distribution Fitting
- Best fit distributions for PRM data
- Compare to DDA results (normalized vs raw)

## Power Analysis
- Standard power analysis
- Power with include_missingness = TRUE
- Compare: how much does missingness reduce power?

## FDR Considerations
- With 285 peptides, multiple testing matters
- FDR-aware power analysis
- Compare nominal vs FDR-adjusted power

## Recommendations
- Sample size recommendations
- Discussion of targeted vs discovery tradeoffs
- Guidance on handling missing data
```

### Key Code Blocks

```r
# Data prep - average technical replicates
prm <- read.csv("../../sample_data/prm_data.csv")
pilot <- prm |>
  filter(timepoint %in% c(0, 6)) |>
  group_by(peptide_modified_sequence, genotype, timepoint, bio_rep) |>
  summarise(abundance = mean(total_area, na.rm = TRUE), .groups = "drop") |>
  transmute(
    peptide_id = peptide_modified_sequence,
    condition = paste0("t", timepoint),
    abundance = abundance
  )

# Fitting with missingness
fits <- fit_distributions(pilot, "peptide_id", "condition", "abundance")

# Missingness analysis
plot_missingness(fits)
mnar_peps <- get_mnar_peptides(fits, threshold = 2)

# Power with and without missingness
power_standard <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                                  find = "power", include_missingness = FALSE)
power_with_na <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                                 find = "power", include_missingness = TRUE)

# FDR-aware
power_fdr <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                            find = "power", apply_fdr = TRUE, prop_null = 0.8)
```

---

## Implementation Steps

1. **Create directory structure**
   ```bash
   mkdir -p inst/examples
   ```

2. **Create DDA analysis Rmd**
   - Write full document following outline above
   - Use `output: html_document` with self_contained: true

3. **Render DDA analysis**
   ```r
   rmarkdown::render("inst/examples/dda-time-course-power.Rmd")
   ```

4. **Create PRM analysis Rmd**
   - Write full document following outline above
   - Emphasize missingness features

5. **Render PRM analysis**
   ```r
   rmarkdown::render("inst/examples/prm-genotype-power.Rmd")
   ```

6. **Verify outputs**
   - Check HTML renders correctly
   - Verify data paths work from inst/examples/

7. **Commit**
   - Add inst/examples/*.Rmd and inst/examples/*.html
   - Do NOT add sample_data/ to git (keep private)

---

## Notes

- **Data paths**: Use relative paths `../../sample_data/` from inst/examples/
- **Caching**: Use `cache=TRUE` for slow chunks during development
- **Figures**: Set reasonable fig.width/fig.height for HTML output
- **Tone**: Professional but accessible; this is a showcase of peppwR capabilities
- **Reproducibility**: Set seeds for all random operations

---

## Verification

After implementation, confirm:
- [ ] Both Rmd files render without error
- [ ] All peppwR functions demonstrated work correctly
- [ ] Plots display properly in HTML
- [ ] Recommendations are evidence-based and reasonable
- [ ] No private data paths or sensitive info in committed files
