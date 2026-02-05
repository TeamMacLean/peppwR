# peppwR: Simulation-based power analysis for phosphoproteomics experiments

## Target Journal
**Bioinformatics** (Oxford) - Application Note

## Format Requirements
- ~2,000 words + 1 figure (max 4 pages / 2,600 words)
- Software name in title
- Supplementary material allowed

---

## Abstract (~150 words)

**Background:** Phosphoproteomics experiments require careful experimental design, yet no dedicated tools exist for power analysis that account for the unique statistical properties of peptide-level quantification data.

**Results:** We present peppwR, an R package for simulation-based power analysis tailored to phosphoproteomics. peppwR fits distributions to pilot data on a per-peptide basis, capturing the heterogeneity in abundance and variance across the peptidome. It supports multiple statistical tests (Wilcoxon, bootstrap-t, Bayes factor), handles missing-not-at-random (MNAR) patterns common in mass spectrometry data, and provides FDR-aware power calculations. Users can determine required sample sizes, estimate power for planned experiments, or identify minimum detectable effect sizes.

**Availability:** peppwR is freely available from GitHub (https://github.com/TeamMacLean/peppwR). Documentation at https://teammaclean.github.io/peppwR/.

---

## 1. Introduction (~400 words)

### The power analysis gap in proteomics

- Statistical power = probability of detecting true effects
- Under-powered experiments waste resources; over-powered experiments are inefficient
- Power analysis standard in clinical trials, genomics—but neglected in proteomics

### Why generic tools fail for proteomics

- **Heterogeneity**: Peptides vary dramatically in abundance and variance (cite Perseus paper for context on proteomics complexity)
- **Non-normality**: Abundance distributions often skewed (gamma, lognormal)
- **Missing data**: MNAR patterns—low abundance peptides systematically missing
- **Multiple testing**: Thousands of peptides tested simultaneously

### Existing tools and their limitations

| Tool | Limitation |
|------|------------|
| Perseus | No power analysis functionality despite being the most widely-used platform |
| clippda | Designed for SELDI-TOF peak data; assumes homogeneous variance |
| MultiPower | Multi-omics focus; requires multiple data types |
| G*Power | Generic; assumes normality; no proteomics-specific features |

### peppwR addresses these gaps

- Per-peptide distribution fitting captures heterogeneity
- Simulation-based approach handles non-normality
- MNAR-aware simulations
- FDR-aware power across the peptidome

---

## 2. Implementation (~600 words)

### 2.1 Workflow overview

```
pilot_data → fit_distributions() → power_analysis() → results + plots
```

Two modes:
1. **Aggregate mode**: User specifies distribution (no pilot data needed)
2. **Per-peptide mode**: Fits distributions to each peptide from pilot data

### 2.2 Distribution fitting

- Candidate distributions: gamma, lognormal, normal, inverse Gaussian
- Selection via AIC
- Handles fit failures: exclude, empirical bootstrap, or lognormal fallback
- Uses `fitdistrplus` and `univariateML` packages

### 2.3 Simulation engine

For each peptide:
1. Draw n samples from fitted distribution (control group)
2. Apply multiplicative effect (treatment group)
3. Run statistical test
4. Repeat n_sim times
5. Power = proportion significant

### 2.4 Statistical tests

- **Wilcoxon rank-sum** (default): Non-parametric, robust
- **Bootstrap-t**: Handles non-normality
- **Bayes factor t-test**: Bayesian evidence quantification

### 2.5 Missing data handling

- `compute_missingness()`: Calculates NA rate and MNAR score per peptide
- MNAR detection via correlation between abundance and missingness
- Simulations incorporate realistic missingness patterns

### 2.6 FDR-aware power

- Simulates whole-peptidome experiments
- Applies Benjamini-Hochberg correction
- Reports proportion of true positives discovered at target FDR
- User-configurable proportion of true nulls (`prop_null`)

### 2.7 Three questions answered

| Question | Parameter | Output |
|----------|-----------|--------|
| What sample size? | `find = "sample_size"` | N per group for target power |
| What's my power? | `find = "power"` | Power at given N and effect |
| What can I detect? | `find = "effect_size"` | Minimum detectable fold-change |

---

## 3. Example Application (~450 words)

### 3.1 Simulated phosphoproteomics dataset (~80 words)

- Generate realistic data: 500 peptides with heterogeneous parameters
- Gamma-distributed abundances with varying shape (1.5-5) and rate (0.05-0.2)
- 6 pilot samples per peptide (typical pilot study size)
- Introduce MNAR missingness (~15% overall, correlated with abundance)

```r
set.seed(42)
peptide_params <- tibble(
  peptide = paste0("pep_", sprintf("%04d", 1:500)),
  shape = runif(500, 1.5, 5),
  rate = runif(500, 0.05, 0.2)
)
```

### 3.2 Per-peptide power analysis (~80 words)

```r
fits <- fit_distributions(pilot_data, id = "peptide",
                          group = "condition", value = "abundance")
result <- power_analysis(fits, effect_size = 2,
                         target_power = 0.8, find = "sample_size")
```

**Primary result**: "N = 8 per group required for 80% of peptides to achieve 80% power to detect 2-fold changes"

### 3.3 Comparison 1: Per-peptide vs homogeneous assumption (~100 words)

**Key comparison demonstrating why peppwR matters:**

- Homogeneous approach (generic tools): Pool data, assume uniform variance → N = 5 sufficient
- Per-peptide approach (peppwR): Model heterogeneity → N = 8 required
- Per-peptide approach requires ~60% more samples than homogeneous suggests

**Message**: Generic power tools underestimate sample size because they ignore the substantial heterogeneity in variance across the peptidome.

### 3.4 Comparison 1b: MNAR-aware simulations (~60 words)

peppwR can incorporate observed missingness rates into power simulations:

- Without MNAR: N = 8
- With MNAR (`include_missingness = TRUE`): N = 9-10

**Message**: Accounting for realistic missingness patterns requires additional samples to maintain target power.

### 3.5 Comparison 2: Statistical test choice (~80 words)

Compare Wilcoxon rank-sum vs Bayes factor t-test:

- Wilcoxon: N = 8 (conservative, non-parametric, no distributional assumptions)
- Bayes factor: N = 6-7 (more efficient when parametric assumptions hold)

**Message**: Test choice affects power estimates; users should match their power analysis to their planned statistical approach.

### 3.6 Figure panels (~60 words)

Panels B-D driven by this example data (Panel A is workflow diagram):

- **Panel B**: Histogram of per-peptide power at N=6 (shows heterogeneity)
- **Panel C**: Power curve with Wilcoxon vs Bayes factor comparison lines
- **Panel D**: Power heatmap (N × effect size lookup grid)

---

## 4. Discussion (~200 words)

### Key contributions

- First dedicated power analysis tool for phosphoproteomics
- Per-peptide heterogeneity modeling
- MNAR-aware simulations
- Integration with common analysis workflows

### Limitations

- Computational cost scales with peptide count and n_sim
- Requires pilot data for per-peptide mode (or reasonable assumptions)
- Currently single-threaded

### Future directions

- Integration with pepdiff for complete experimental design → analysis workflow
- Extension to TMT/iTRAQ labeling schemes
- Parallelization for large-scale studies

---

## Figure 1: Multi-panel overview

**Panel A**: Workflow diagram (conceptual)
- Simple flow: `pilot_data → fit_distributions() → power_analysis() → result`
- Clean schematic, not code

**Panel B**: Per-peptide power distribution histogram at N=6 (data-driven)
- X-axis: Power (0-100%)
- Y-axis: Number of peptides
- Shows bimodal distribution: high-abundance peptides near 100%, low-abundance near 20-40%
- Message: Heterogeneity matters—peptides differ dramatically in statistical power

**Panel C**: Power curve with test comparison (data-driven)
- X-axis: Sample size per group (N=3 to N=15)
- Y-axis: % of peptides reaching 80% power
- Two lines: Wilcoxon (solid) vs Bayes factor (dashed)
- Horizontal reference at 80% threshold
- Message: Sample size requirements depend on test choice; Bayes ~15% more efficient

**Panel D**: Power heatmap lookup grid (data-driven)
- Rows: Sample size (N=4 to N=12)
- Columns: Effect size (1.5, 2, 2.5, 3-fold)
- Color: % of peptides reaching 80% power (YlOrRd palette to match logo)
- Message: Practical lookup table for experimental planning

---

## Supplementary Material

- Fully reproducible R code for all analyses and Figure 1
- Benchmarking: Runtime vs peptide count and n_sim
- Additional diagnostic plots

---

## Author Contributions

D.M. conceived the project, designed the software architecture, implemented the package, and wrote the manuscript.

## Acknowledgments

[If applicable]

## Funding

[Grant information]

## Conflict of Interest

None declared.

---

## Rendering Pipeline

### Dependencies

```r
# Add to renv for paper rendering (not package dependencies)
renv::install("rticles")
```

### File Structure

```
paper/
├── paper.Rmd                    # Main manuscript (single source)
├── references.bib               # BibTeX references (exists)
├── outline.md                   # This file
├── example_application_spec.md  # Spec for supplementary Rmd
├── supplementary_example.Rmd    # Reproducible code supplement
└── figures/
    ├── figure_1.pdf             # Combined figure for journal
    └── panel_a_workflow.svg     # Workflow diagram (vector)
```

### Two Output Formats

**1. Bioinformatics Journal Submission**
```r
rmarkdown::render("paper.Rmd", output_format = rticles::bioinformatics_article())
```
- Uses OUP template with copyright footer
- Submit to journal portal

**2. bioRxiv Preprint**
```r
rmarkdown::render("paper.Rmd", output_format = rticles::arxiv_article())
```
- Clean two-column layout, no proprietary footer
- bioRxiv accepts PDF uploads directly

### YAML Header (paper.Rmd)

```yaml
---
title: "peppwR: Simulation-based power analysis for phosphoproteomics experiments"
author:
  - name: Dan MacLean
    affiliation: The Sainsbury Laboratory
    email: dan.maclean@tsl.ac.uk
abstract: |
  **Background:** Phosphoproteomics experiments require careful experimental
  design, yet no dedicated tools exist for power analysis...

  **Availability:** https://github.com/TeamMacLean/peppwR
bibliography: references.bib
output:
  rticles::bioinformatics_article:
    keep_tex: true
  rticles::arxiv_article:
    keep_tex: true
---
```

### Workflow

1. Write content in `paper.Rmd` using standard R Markdown
2. Run `supplementary_example.Rmd` to generate Figure 1
3. Render to both formats from same source
4. Submit preprint to bioRxiv first
5. Submit to Bioinformatics after preprint is posted

### Known Issues

- rticles bioinformatics template has OUP copyright footer
- bioRxiv rejects PDFs with proprietary footers
- Solution: Use `arxiv_article()` for preprint, `bioinformatics_article()` for journal
