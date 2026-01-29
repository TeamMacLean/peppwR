# Next Session: Implement v2 Features

## Context

peppwR v1 is complete with:
- Distribution fitting (`fit_distributions()`)
- Power analysis (aggregate and per-peptide modes)
- Four statistical tests (wilcoxon, bootstrap_t, bayes_t, rankprod)
- Core plots and three vignettes
- 140 passing tests

## Task

Implement v2 features from `feature_plan_v2.md`. Follow the TDD workflow.

## Files to Read First

1. **feature_plan_v2.md** - Detailed implementation plan with:
   - Phase A: Diagnostic plots (5 new plot functions)
   - Phase B: Empirical bootstrap (implement stubbed `on_fit_failure = "empirical"`)
   - Phase C: Missing data handling (MNAR detection, incorporate into simulations)
   - Phase D: FDR-aware mode (whole-peptidome simulation with BH correction)

2. **CLAUDE.md** - Package specification (updated with v2 phases)

3. **for_CLAUDE.md** - Development workflow and code style

## Recommended Approach

### Start with Phase A (Diagnostic Plots)

These are independent and will help build familiarity with the codebase.

**TDD workflow for each plot function:**

1. Write failing test in `tests/testthat/test-plots-diagnostic.R`:
```r
test_that("plot_density_overlay returns ggplot", {
  fits <- fit_distributions(test_data, "peptide_id", "condition", "abundance")
  p <- plot_density_overlay(fits)
  expect_s3_class(p, "ggplot")
})
```

2. Implement in `R/plots.R` following existing patterns:
   - Use `theme_minimal()`
   - Use `ggplot2::` namespacing
   - Return ggplot objects

3. Run tests: `devtools::test(filter = "plots")`

### Then Phase B (Empirical Bootstrap)

The stub is at `R/power.R` lines 202-206. Need:
- `simulate_empirical()` in `R/simulation.R`
- `run_power_sim_empirical()` in `R/simulation.R`
- Wire up in `power_analysis.peppwr_fits()`

### Then Phase C (Missing Data)

Key insight: Model missingness, don't impute. Add `missingness` slot to `peppwr_fits`:
```r
missingness = tibble(
  peptide_idx, n_total, n_missing, na_rate,
  mnar_score,    # Positive = low values more likely missing
  mnar_pvalue
)
```

### Finally Phase D (FDR)

Depends on Phase C. Simulate all peptides together, apply `p.adjust(method = "BH")`.

## Verification

After each phase:
```r
devtools::test()
devtools::check(vignettes = FALSE)  # Skip vignettes for speed
```

## Quick Start Command

```
/ralph-loop "Implement Phase A diagnostic plots for peppwR.

## Spec
See feature_plan_v2.md Phase A section.

## Files
- R/plots.R (add new functions)
- tests/testthat/test-plots-diagnostic.R (new file)

## Verification
Rscript -e 'devtools::test(filter = \"plots\")'

## Success
All plot tests pass, functions return valid ggplot objects." --max-iterations 15
```