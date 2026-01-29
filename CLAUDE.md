# peppwR

Power analysis for phosphopeptide abundance hypothesis tests via simulation.

## Purpose

Help researchers answer:
1. "What sample size do I need for 80% power to detect a given effect?"
2. "With N samples per group, what's my power to detect a 2-fold change?"
3. "Given my sample size and target power, what's the minimum detectable effect?"

## Core Workflow

```
pilot_data → fit_distributions() → power_analysis() → results/recommendation
```

Two modes of operation:

### From Scratch (Aggregate Mode)
No pilot data required. User specifies assumed distribution and parameters (or uses sensible proteomics defaults). Simulates a "typical" peptide.

```r
result <- power_analysis(
  distribution = "gamma",
  params = list(shape = 2, rate = 0.1),
  effect_size = 2,
  target_power = 0.8,
  find = "sample_size"
)
```

### With Pilot Data (Per-Peptide Mode)
Uses real pilot data. Fits distributions to each peptide, then simulates using actual fitted parameters. Reports power across the peptidome.

```r
fits <- fit_distributions(pilot_data, id = "peptide", group = "condition", value = "abundance")
result <- power_analysis(fits, effect_size = 2, target_power = 0.8, find = "sample_size")
# Output: "73% of peptides achieve 80% power with N=6 per group"
```

## S3 Classes

### `peppwr_fits`
Distribution fitting results. Contains nested data with per-peptide fit statistics.

**Slots:**
- `data` - Original data (nested)
- `fits` - Fit results per peptide (distribution, parameters, AIC, loglik)
- `best` - Best-fitting distribution per peptide
- `call` - Original function call

**Methods:**
- `print()` - Summary of fits (best distribution counts, failure rate)
- `plot()` - Distribution evaluation dashboard (best fits, failed fits)
- `summary()` - Detailed fit statistics

### `peppwr_power`
Power analysis results.

**Slots:**
- `mode` - "aggregate" or "per_peptide"
- `question` - Which question was asked: "sample_size", "power", or "effect_size"
- `answer` - The computed answer
- `simulations` - Simulation details (for diagnostics)
- `params` - Input parameters used
- `call` - Original function call

**Methods:**
- `print()` - Clear answer to the question asked
- `plot()` - Power curves, sample size recommendations
- `summary()` - Detailed simulation results

## Core Functions

### Distribution Fitting

```r
fit_distributions(data, id, group, value, distributions = "all")
```
- Fits candidate distributions to each peptide's values
- Returns `peppwr_fits` object
- Default distributions: gamma, normal, lognormal, inverse gaussian, etc.

### Power Analysis

```r
power_analysis(x, ...)
```

Generic with methods for:
- `peppwr_fits` - Per-peptide mode using fitted parameters
- `character` - Aggregate mode with specified distribution
- `default` - Aggregate mode with proteomics defaults

**Key parameters:**
- `effect_size` - Fold change (e.g., 2 for 2-fold)
- `n_per_group` - Sample size per group (if known)
- `target_power` - Desired power (default 0.8)
- `alpha` - Significance level (default 0.05)
- `find` - What to solve for: "sample_size", "power", or "effect_size"
- `test` - Statistical test to use (see below)
- `n_sim` - Number of simulations (default 1000)
- `on_fit_failure` - How to handle peptides with no valid fit: "exclude" (default), "empirical", "lognormal"

### Simulation Engine (Internal)

```r
simulate_experiment(distribution, params, n_per_group, effect_size, n_sim)
run_test(control, treatment, test)
```

## Plots

### Core plots (v1)

**1. `plot(peppwr_fits)` - Distribution summary**
Bar chart showing best-fit distribution counts by metric (AIC/LogLik). Quick sanity check that fitting worked.

**2. `plot(peppwr_power)` aggregate mode - Power curve with recommendation**
Power vs sample size curve with:
- Horizontal line at target power (default 80%)
- Vertical annotation at recommended N
- Clear title: "Recommended sample size: N=X per group"

**3. `plot(peppwr_power)` per-peptide mode - Peptide threshold curve**
"% of peptides reaching target power" vs sample size. Directly answers "how many of my peptides are well-powered at each N?"

### Diagnostic plots (future versions)

| Plot | Purpose |
|------|---------|
| Failed fits bar chart | Data quality check - which distributions fail |
| Density overlay | Visual fit assessment - observed vs fitted |
| QQ plots | Goodness-of-fit diagnostics |
| Power vs effect size | Sensitivity analysis at fixed N |
| Power heatmap | N × effect size lookup grid |
| Parameter distributions | Distribution of fitted params across peptidome |

## Statistical Tests

Configurable via `test` parameter:

| Test | ID | Description |
|------|----|-------------|
| Wilcoxon rank-sum | `"wilcoxon"` | Non-parametric, robust |
| Bootstrap-t | `"bootstrap_t"` | Resampling-based, handles non-normality |
| Bayes factor t-test | `"bayes_t"` | Bayesian evidence for effect |
| Rank Products | `"rankprod"` | Designed for omics, handles small samples |

Default: `"wilcoxon"`

## Simulation Approach

1. **Control group**: Draw `n_per_group` samples from fitted distribution
2. **Treatment group**: Draw from same distribution, apply multiplicative effect
3. **Test**: Run specified statistical test
4. **Repeat**: `n_sim` iterations
5. **Power**: Proportion of significant results

For per-peptide mode, repeat across all peptides and report distribution of power values.

## File Structure

```
R/
  fits.R          - Distribution fitting (fit_distributions, single_fit, do_fits)
  power.R         - Power analysis (power_analysis, simulate_experiment)
  tests.R         - Statistical tests (run_test, wilcoxon, bootstrap_t, etc.)
  plots.R         - Visualization (plot methods, evaldist)
  classes.R       - S3 class constructors and validators
  utils.R         - Helpers (distribution utilities, parameter extraction)

tests/testthat/
  test-fits.R     - Distribution fitting tests
  test-power.R    - Power analysis tests
  test-tests.R    - Statistical test implementations
  test-plots.R    - Plot output tests
```

## Project Plan

Implementation order with dependencies:

### Phase 1: Class Foundation
1. **`peppwr_fits` constructor** - `new_peppwr_fits()` with validation
2. **`print.peppwr_fits`** - Summary output (fit counts, failure rate)
3. **`peppwr_power` constructor** - `new_peppwr_power()` with validation
4. **`print.peppwr_power`** - Clear answer to question asked

### Phase 2: Distribution Fitting (refactor existing)
5. **Refactor `fit_distributions()`** - Wraps existing logic, returns `peppwr_fits`
6. **`plot.peppwr_fits`** - Distribution summary bar chart (refactor `evaldist()`)

### Phase 3: Simulation Engine
7. **`simulate_experiment()`** - Draw from distribution, apply effect, return samples
8. **`test_wilcoxon()`** - Wilcoxon rank-sum test wrapper
9. **`run_power_sim()`** - Run n_sim iterations, return power estimate

### Phase 4: Power Analysis - Aggregate Mode
10. **`power_analysis.default`** - Aggregate mode, find = "power"
11. **Extend for find = "sample_size"** - Search over N values
12. **Extend for find = "effect_size"** - Search over effect values

### Phase 5: Power Analysis - Per-Peptide Mode
13. **`power_analysis.peppwr_fits`** - Per-peptide mode using fitted params
14. **`on_fit_failure` handling** - exclude/empirical/lognormal options

### Phase 6: Core Plots
15. **`plot.peppwr_power` aggregate** - Power curve with recommendation annotation
16. **`plot.peppwr_power` per-peptide** - % peptides at threshold vs N

### Phase 7: Additional Statistical Tests
17. **`test_bootstrap_t()`** - Bootstrap-t implementation
18. **`test_bayes_t()`** - Bayes factor t-test (via BayesFactor)
19. **`test_rankprod()`** - Rank products implementation

### Phase 8: Polish (COMPLETE)
20. **`summary()` methods** - Detailed output for both classes
21. **Vignettes** - Three comprehensive vignettes (getting-started, power-analysis-workflow, benchmarking)

---

## v2 Feature Phases

See `feature_plan_v2.md` for detailed specifications.

### Phase A: Diagnostic Plots
- `plot_density_overlay()` - Observed histogram + fitted density curve
- `plot_qq()` - QQ plots for goodness-of-fit
- `plot_power_heatmap()` - N × effect size lookup grid
- `plot_power_vs_effect()` - Sensitivity at fixed N
- `plot_param_distribution()` - Fitted params across peptidome

### Phase B: Empirical Bootstrap
- `simulate_empirical()` - Bootstrap resample from observed data
- `run_power_sim_empirical()` - Power sim using bootstrap
- Implement `on_fit_failure = "empirical"` (currently stubbed)

### Phase C: Missing Data Handling
**Philosophy: Track and model missingness, never impute**

- `compute_missingness()` - Calculate NA rate, MNAR score, MNAR p-value
- Extend `peppwr_fits` with `missingness` slot
- MNAR detection (Missing Not At Random - when low values systematically missing)
- `simulate_with_missingness()` - Incorporate NA rates into simulations
- `plot_missingness()` - NA rate and MNAR score distributions

### Phase D: FDR-Aware Mode
- `run_power_sim_fdr()` - Whole-peptidome simulation with BH correction
- Add `apply_fdr`, `prop_null`, `fdr_threshold` params to `power_analysis.peppwr_fits()`
- User-configurable `prop_null` (default 0.9 = 90% true nulls)

**Dependencies**: Phase C should complete before Phase D

---

Each item is a candidate for the TDD → Ralph Loop workflow.

## Development

See `for_CLAUDE.md` for code style and `semi-autonomous-feature-development.md` for workflow.

**Key principles:**
- Tidyverse style (pipes, dplyr verbs)
- Explicit namespacing (`dplyr::filter()`)
- S3 classes with print/plot methods
- Discuss → TDD → Ralph Loop workflow

## Dependencies

Current:
- `fitdistrplus`, `univariateML` - Distribution fitting
- `dplyr`, `tidyr`, `purrr`, `tibble` - Data manipulation
- `ggplot2`, `cowplot`, `RColorBrewer` - Visualization

Likely additions:
- `BayesFactor` - Bayes factor t-tests
- `RankProd` or custom implementation - Rank products test
- `boot` - Bootstrap methods

## Design Decisions

### Handling fit failures
Configurable via `on_fit_failure` parameter:
- `"exclude"` (default) - Skip peptides with no valid fit, report count in output
- `"empirical"` - Bootstrap resample from observed values
- `"lognormal"` - Fallback to lognormal with moment-matched parameters

### Multiple testing / FDR
**v1: Single-test power only.** Document clearly that reported power is per-peptide at nominal alpha, not FDR-adjusted.

**Future enhancement:** FDR-aware mode that simulates full peptidome (null + changed) and applies correction.

### Performance
**v1: Single-threaded.** Rely on R's vectorization. Add parallelization (via `future`/`furrr`) only if benchmarks show it's needed.