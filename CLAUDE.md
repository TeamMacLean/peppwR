# peppwR

Power analysis for phosphopeptide abundance hypothesis tests via simulation.

**Repository:** https://github.com/TeamMacLean/peppwR
**Documentation:** https://teammaclean.github.io/peppwR/
**Version:** 0.1.0 (CRAN-ready)

## Related Package

**pepdiff** - Companion package for differential abundance analysis of phosphopeptide data. While peppwR handles experimental design (power analysis), pepdiff handles the actual analysis once data is collected.

## Purpose

Help proteomics researchers answer experimental design questions:

1. **"What sample size do I need?"** → `find = "sample_size"`
2. **"What's my power?"** → `find = "power"`
3. **"What's the minimum detectable effect?"** → `find = "effect_size"`

## Core Workflow

```
pilot_data → fit_distributions() → power_analysis() → results + plots
```

### Aggregate Mode (No Pilot Data)

```r
result <- power_analysis(
  distribution = "gamma",
  params = list(shape = 2, rate = 0.1),
  effect_size = 2,
  target_power = 0.8,
  find = "sample_size"
)
```

### Per-Peptide Mode (With Pilot Data)

```r
fits <- fit_distributions(pilot_data, id = "peptide", group = "condition", value = "abundance")
result <- power_analysis(fits, effect_size = 2, target_power = 0.8, find = "sample_size")
# Output: "73% of peptides achieve 80% power with N=6 per group"
```

## Key Features

| Feature | Description |
|---------|-------------|
| **Distribution fitting** | gamma, lognormal, normal, inverse gaussian via `fitdistrplus`/`univariateML` |
| **Statistical tests** | `wilcoxon` (default), `bootstrap_t`, `bayes_t` |
| **Missing data** | Dataset-level MNAR detection (abundance vs NA rate correlation), missingness-aware simulations |
| **FDR-aware mode** | Whole-peptidome simulation with BH correction |
| **Empirical bootstrap** | Fallback when parametric fitting fails |
| **Diagnostic plots** | Density overlay, QQ, power heatmap, missingness |

## S3 Classes

### `peppwr_fits`
Distribution fitting results with `print()`, `plot()`, `summary()` methods.

### `peppwr_power`
Power analysis results with `print()`, `plot()`, `summary()` methods.

## File Structure

```
R/
├── classes.R       # S3 class constructors and validators
├── fits.R          # fit_distributions(), missingness computation, MNAR detection
├── plots.R         # All plot methods and diagnostic plots
├── power.R         # power_analysis() and methods
├── simulation.R    # Simulation engine, run_power_sim*()
├── tests.R         # Statistical test wrappers
└── utils.R         # Distribution utilities, helpers

vignettes/
├── getting-started.Rmd        # CRAN vignette (quick)
├── find-fits.Rmd              # CRAN vignette (quick)
└── articles/                  # pkgdown articles (not in CRAN package)
    ├── benchmarking.Rmd
    ├── power-analysis-workflow.Rmd
    ├── dda-time-course-power.Rmd
    └── prm-genotype-power.Rmd

sample_data/                   # Example datasets (in .Rbuildignore)
├── dda_data.csv               # Anonymized DDA phosphoproteomics
└── prm_data.csv               # Anonymized PRM phosphoproteomics
```

## Distribution & CI

### GitHub Actions
- **R-CMD-check.yaml** - Multi-platform R CMD check (Ubuntu, macOS, Windows)
- **pkgdown.yaml** - Builds and deploys documentation to gh-pages

### CRAN Preparation
- Version 0.1.0, passes `R CMD check --as-cran` with 0 errors, 0 warnings
- Heavy vignettes moved to pkgdown articles to reduce check time
- Package size optimized (<5MB)

### pkgdown Site
- Hosted at https://teammaclean.github.io/peppwR/
- Includes all articles, function reference, and hex logo

## Development Style

- **Tidyverse style** - Pipes, dplyr verbs
- **Explicit namespacing** - `dplyr::filter()`, `ggplot2::ggplot()`
- **S3 OOP** - Classes with print/plot/summary methods
- **TDD** - testthat with 299 tests

## Development Workflow for Future Changes

See `semi-autonomous-feature-development.md` for detailed workflow.

### Discuss → TDD → Ralph Loop

1. **Discuss** - Human describes feature/bug, Claude asks clarifying questions, agree on scope
2. **TDD** - Write failing test first (the test IS the spec)
3. **Ralph Loop** - Claude iterates autonomously until tests pass

### Key Principles

- **Tests are the contract** - No ambiguity about completion
- **No implementation until test fails** - Red → Green → Refactor
- **Clear context before implementation** - Commit test, start fresh session
- **Self-contained prompts** - Reference files, not discussion history

### For Bug Fixes / Features

```
1. Discuss requirements
2. Write failing test in tests/testthat/
3. Commit the test
4. /clear or new session
5. /ralph-loop with verification command: devtools::test(filter = "test-name")
6. Human smoke test
```

## Dependencies

**Imports:**
- `fitdistrplus`, `univariateML` - Distribution fitting
- `dplyr`, `tidyr`, `purrr`, `tibble` - Data manipulation
- `ggplot2`, `cowplot`, `RColorBrewer`, `scales` - Visualization

**Suggests:**
- `testthat`, `knitr`, `rmarkdown`, `pkgdown`, `bench`

## Design Decisions

### Fit Failures
Configurable via `on_fit_failure`:
- `"exclude"` (default) - Skip peptides, report count
- `"empirical"` - Bootstrap resample from observed values
- `"lognormal"` - Fallback with moment-matched parameters

### FDR-Aware Mode
Use `apply_fdr = TRUE` in per-peptide mode for whole-peptidome simulation with Benjamini-Hochberg correction. Configure `prop_null` (default 0.9).

### Performance
Single-threaded with R vectorization. Benchmarks show acceptable performance for typical datasets (1000s of peptides).

## Project Status

**COMPLETE** - All planned features implemented:

- Core power analysis (aggregate and per-peptide modes)
- All three questions: power, sample_size, effect_size
- Distribution fitting with multiple candidates
- Statistical tests: wilcoxon, bootstrap_t, bayes_t
- Diagnostic plots: density overlay, QQ, heatmap, param distribution
- Missingness handling with dataset-level MNAR detection
- FDR-aware power analysis
- Empirical bootstrap fallback
- Real-world examples (DDA and PRM)
- Comprehensive documentation and pkgdown site
- Hex sticker logo (chilli pepper power curve)
