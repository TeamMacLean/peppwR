# peppwR

**Power analysis for phosphopeptide abundance hypothesis tests via
simulation.**

peppwR helps proteomics researchers answer critical experimental design
questions:

- **“What sample size do I need?”** — Find the N required for 80% power
  to detect a given effect
- **“What’s my power?”** — Calculate power for a given sample size and
  effect
- **“What can I detect?”** — Find the minimum detectable effect size for
  your design

## Features

- **Distribution fitting** — Fit gamma, lognormal, normal, and other
  distributions to pilot data
- **Per-peptide analysis** — Get power estimates across your entire
  peptidome
- **Multiple statistical tests** — Wilcoxon, bootstrap-t, and Bayes
  factor t-tests
- **Missing data handling** — Model MNAR (Missing Not At Random)
  patterns common in proteomics
- **FDR-aware mode** — Account for multiple testing correction in power
  calculations
- **Rich visualizations** — Power curves, heatmaps, QQ plots, and
  diagnostic plots

## Installation

Install the development version from GitHub:

``` r
# install.packages("pak")
pak::pak("TeamMacLean/peppwR")
```

## Quick Start

### From Pilot Data (Per-Peptide Mode)

``` r
library(peppwR)

# Fit distributions to your pilot data
fits <- fit_distributions(
  pilot_data,
  id = "peptide",
  group = "condition",
  value = "abundance"
)

# Find required sample size for 80% power to detect 2-fold change
result <- power_analysis(
  fits,
  effect_size = 2,
  target_power = 0.8,
  find = "sample_size"
)

print(result)
#> Per-peptide power analysis
#> Question: sample_size
#> Answer: N = 6 per group
#> 73% of peptides achieve 80% power

plot(result)
```

### From Scratch (Aggregate Mode)

No pilot data? Specify assumed distribution parameters:

``` r
result <- power_analysis(
  distribution = "gamma",
  params = list(shape = 2, rate = 0.1),
  effect_size = 2,
  target_power = 0.8,
  find = "sample_size"
)
```

## Documentation

- **[Getting
  Started](https://teammaclean.github.io/peppwR/articles/getting-started.html)**
  — Introduction and basic usage
- **[Power Analysis
  Workflow](https://teammaclean.github.io/peppwR/articles/power-analysis-workflow.html)**
  — Complete workflow with real data
- **[Benchmarking](https://teammaclean.github.io/peppwR/articles/benchmarking.html)**
  — Performance characteristics
- **[Function
  Reference](https://teammaclean.github.io/peppwR/reference/index.html)**
  — Full API documentation

### Real-World Examples

- [DDA Time Course
  Analysis](https://teammaclean.github.io/peppwR/articles/dda-time-course-power.html)
- [PRM Genotype
  Analysis](https://teammaclean.github.io/peppwR/articles/prm-genotype-power.html)

## Workflow Overview

``` mermaid
flowchart TD
    A[pilot_data] --> B[fit_distributions]
    B --> C[power_analysis]
    C --> D{find = ?}
    D --> E["power"]
    D --> F["sample_size"]
    D --> G["effect_size"]
    E --> H[results + plots]
    F --> H
    G --> H

    style A fill:#FFFFCC,stroke:#BD0026
    style B fill:#FD8D3C,stroke:#BD0026,color:#fff
    style C fill:#FD8D3C,stroke:#BD0026,color:#fff
    style H fill:#FFFFCC,stroke:#BD0026
```

## Citation

If you use peppwR in your research, please cite:

``` R
MacLean, D. (2026). peppwR: Power Analysis for Phosphopeptide Abundance
Hypothesis Tests. R package version 0.1.0.
https://github.com/TeamMacLean/peppwR
```

## Contributing

Contributions are welcome! Please open an
[issue](https://github.com/TeamMacLean/peppwR/issues) or submit a pull
request.

## License

MIT © [Dan MacLean](https://github.com/danmaclean)
