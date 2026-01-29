# peppwR Feature Addition Plan v2

## Overview

Four feature sets to add to peppwR:
1. **Diagnostic Plots** - Visual fit assessment and power exploration
2. **Empirical Bootstrap** - `on_fit_failure = "empirical"` implementation
3. **Missing Data Handling** - Track NA rates, detect MNAR, incorporate into simulations
4. **FDR-Aware Mode** - Multiple testing correction in power analysis

## Design Decisions

### Missing Data Philosophy
- **NO imputation** - ever
- Track and report NA rates per peptide
- Detect MNAR (Missing Not At Random) patterns - when low values are systematically missing
- Incorporate missingness rates into power simulations
- Always try to fit, report % failing due to NAs

### FDR Mode
- User-configurable `prop_null` parameter (default 0.9 = 90% true nulls)
- Simulate full peptidome together, apply Benjamini-Hochberg correction
- Report power after FDR adjustment

---

## Implementation Phases

### Phase A: Diagnostic Plots

**Files**: `R/plots.R`, new `tests/testthat/test-plots-diagnostic.R`

| Task | Function | Description |
|------|----------|-------------|
| A1 | `plot_density_overlay(fits, peptide_id=NULL, n_overlay=6)` | Observed histogram + fitted density curve |
| A2 | `plot_qq(fits, peptide_id=NULL, n_plots=6)` | QQ plots for goodness-of-fit |
| A3 | `plot_power_heatmap(distribution, params, n_range, effect_range)` | N × effect size grid |
| A4 | `plot_power_vs_effect(power_result, effect_range)` | Sensitivity at fixed N |
| A5 | `plot_param_distribution(fits)` | Fitted params across peptidome |

### Phase B: Empirical Bootstrap

**Files**: `R/simulation.R`, `R/power.R`

| Task | Function | Description |
|------|----------|-------------|
| B1 | `simulate_empirical(raw_data, n_per_group, effect_size)` | Bootstrap resample from observed |
| B2 | `run_power_sim_empirical(raw_data, n_per_group, ...)` | Power sim using bootstrap |
| B3 | Update `power_analysis.peppwr_fits()` lines 202-206 | Wire up "empirical" option |

### Phase C: Missing Data Handling

**Files**: `R/fits.R`, `R/classes.R`, `R/simulation.R`, `R/plots.R`

| Task | Function | Description |
|------|----------|-------------|
| C1 | `compute_missingness(values)` | Calculate NA rate, MNAR score, MNAR p-value |
| C2 | Extend `fit_distributions()` | Add `missingness` slot to peppwr_fits |
| C3 | Update `print/summary.peppwr_fits()` | Report missingness statistics |
| C4 | `simulate_with_missingness(dist, params, ..., na_rate, mnar_score)` | Introduce realistic missingness |
| C5 | `run_power_sim_with_missingness()` | Power sim incorporating NA rates |
| C6 | Add `include_missingness` param to `power_analysis.peppwr_fits()` | Use peptide-specific NA rates |
| C7 | `plot_missingness(fits)` | NA rate distribution + MNAR score plot |

### Phase D: FDR-Aware Mode

**Files**: `R/simulation.R`, `R/power.R`, `R/classes.R`

| Task | Function | Description |
|------|----------|-------------|
| D1 | `run_power_sim_fdr(fits, effect_size, n_per_group, prop_null, fdr_threshold)` | Whole-peptidome FDR sim |
| D2 | Add `apply_fdr`, `prop_null`, `fdr_threshold` params to `power_analysis.peppwr_fits()` | Enable FDR mode |
| D3 | Update `print.peppwr_power()` | Describe FDR-adjusted results |

---

## Key Data Structures

### Extended peppwr_fits (Phase C)
```r
peppwr_fits <- list(
  data = nested_tibble,
  fits = list(),
  best = character(),
  call = call,
  missingness = tibble(        # NEW
    peptide_idx, n_total, n_missing, na_rate,
    mnar_score,    # Positive = low values more likely missing
    mnar_pvalue    # Test for MNAR pattern
  )
)
```

### MNAR Detection Logic
```r
# Under MCAR: mean rank of observed values = (n+1)/2
# Under MNAR (low values missing): mean rank is higher
# mnar_score = normalized deviation from expected
mnar_score <- (observed_mean_rank - expected_mean_rank) / se
```

### FDR Simulation Structure
```r
# Each iteration:
# 1. Assign peptides to null (no effect) or alternative (has effect)
# 2. Simulate all peptides
# 3. Run tests, collect p-values
# 4. Apply p.adjust(method = "BH")
# 5. Power = proportion of true alternatives detected at fdr_threshold
```

---

## Dependencies Between Phases

```
Phase A (Plots) -----> Independent
Phase B (Empirical) -> Independent
Phase C (Missing) ---> Should complete before D
Phase D (FDR) -------> Depends on C (uses missingness info)
```

Recommended order: **A → B → C → D** (or A and B in parallel)

---

## Critical Files

| File | Changes |
|------|---------|
| `R/simulation.R` | Add `simulate_empirical`, `simulate_with_missingness`, `run_power_sim_fdr` |
| `R/power.R` | Update lines 202-206; add `apply_fdr`, `include_missingness` params |
| `R/fits.R` | Add `compute_missingness`; extend `fit_distributions` |
| `R/plots.R` | Add 6 new diagnostic plot functions |
| `R/classes.R` | Extend `peppwr_fits` structure; update print/summary methods |

---

## Verification

After each phase:
```r
devtools::test()           # All tests pass
devtools::check(vignettes = FALSE)  # 0 errors, 0 warnings
```

After all phases:
```r
# Manual verification
fits <- fit_distributions(pilot_data, "peptide_id", "condition", "abundance")
print(fits)  # Should show missingness summary

# Test FDR mode
result <- power_analysis(fits, effect_size = 2, n_per_group = 6,
                         find = "power", apply_fdr = TRUE, prop_null = 0.9)
print(result)  # Should show FDR-adjusted power

# Test empirical bootstrap
result <- power_analysis(fits, effect_size = 2, n_per_group = 6,
                         find = "power", on_fit_failure = "empirical")

# Test missingness incorporation
result <- power_analysis(fits, effect_size = 2, n_per_group = 6,
                         find = "power", include_missingness = TRUE)
```

---

## Estimated Scope

- **Phase A**: 5 functions, ~200 lines
- **Phase B**: 3 functions, ~80 lines
- **Phase C**: 7 tasks, ~300 lines
- **Phase D**: 3 tasks, ~150 lines

Total: ~730 lines of new code + tests