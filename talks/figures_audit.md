# Figures audit — `talks/slides.qmd`

Read-only audit. For each slide in `slides.qmd` that carries a `[Figure: …]` placeholder, a `[Cartoon: …]`, a `[Three thumbnails: …]` or a `[Small visual: …]` block, this table notes whether a usable image already exists locally or remotely, or whether the slide needs a fresh hand-drawn figure. Slide numbers below follow `talks/slides_outline.md` (the .qmd does not number slides explicitly).

Local paths are repo-relative. Remote URLs point at `https://teammaclean.github.io/pepdiff/` or `https://teammaclean.github.io/peppwR/` pkgdown sites.

## Power analysis section

| Slide # | Slide title | Slot description (compressed) | Recommendation | Source | Notes |
|---|---|---|---|---|---|
| 6 | The four levers | Annotated two-curve diagram with $\delta$, $\sigma$, $n$, $\alpha$ labelled at the geometric feature each controls; inset showing sampling distribution narrowing with $n$ | draw-fresh |  | Pure conceptual anatomy. Two overlapping Gaussians; horizontal arrow for $\delta$ between peaks; half-width brackets for $\sigma$; tail shading for $\alpha$ under the null; inset of sampling-distribution-of-mean shrinking for $n$. No existing peppwR/pepdiff figure does this. |
| 7 | Two curves, low power and high power | Side-by-side two-Gaussian overlap, low-power vs high-power | draw-fresh |  | Same family of diagram as slide 6 but two-panel. Could be drawn in ggplot easily; no candidate in the docs. |
| 9 | What $p$ actually is | Null distribution with shaded tail past observed mean | draw-fresh |  | Single bell curve with a vertical line at the observed mean and tail shading. Conceptual; not in the docs. |
| 12 | The classical power curve | Power vs $n$ for $D = 0.2, 0.5, 0.8$ from `pwr.t.test`; dashed line at $0.8$; crossings near $n \approx 393, 64, 26$ | draw-fresh |  | Quick `pwr` + ggplot script. peppwR's per-peptide curves (e.g. `plot-pp-sample-size-1.png`) have the same shape but the axes/values and the three-curve family are different — would mislead. |
| 14 | Trick one: merge technical replicates honestly | Cartoon: one mouse + four tubes vs four mice + four tubes | draw-fresh |  | Cartoon; no candidate. Could be drawn in Inkscape / Affinity / drawn-by-hand on iPad. |
| 15 | Ordinal scores need a different merge | Small visual: 0–6 score column → rank column → per-biological-rep mean rank | draw-fresh |  | Tiny three-column schematic. Cleanest as hand-set tables in Quarto rather than an image. |
| 18 | What the re-analysis showed | Original vs re-analysed HR-score figure (Figure 4 in prose) | draw-fresh |  | The original published figure is from an anonymous external paper — not in this repo. Will need a fresh re-render from the underlying data, or a pair of stylised mock violins / stripcharts to convey 75/120 → 42/120. Verify with Dan whether real data is available. |
| 19 | When the textbook recipe runs out | Three thumbnails: many overlapping curves; right-skewed histogram with zero-spike; peptide-by-sample tile matrix with grey gaps | adapt-existing | `vignettes/articles/power-analysis-workflow_files/figure-html/density-overlay-1.png` (middle thumbnail) + draw-fresh for the other two | Density-overlay vignette figure is a close match for the right-skewed histogram thumbnail, though it shows multiple peptides as small multiples rather than a single distribution; crop or remake. The "many overlapping curves" and "missingness tile matrix" thumbnails have no direct candidate — draw fresh. |

## peppwR section

| Slide # | Slide title | Slot description (compressed) | Recommendation | Source | Notes |
|---|---|---|---|---|---|
| 21 | Start with the histogram, not the family | Histogram of DDA *Arabidopsis* phosphopeptide intensities (Figure 5 in prose) | use-existing | `https://teammaclean.github.io/peppwR/articles/dda-time-course-power_files/figure-html/eda-plot-1.png` | Vignette EDA plot shows log2-abundance density across timepoints — close fit. Slight wording mismatch: it is a density/violin, the prose calls it a histogram. Verify with Dan whether the density form is acceptable; otherwise rerender as `hist()` from the DDA pilot data. |
| 27 | Missingness, and why it is not random | Mean abundance vs NA rate scatter (Figure 6 in prose) | use-existing | `vignettes/articles/power-analysis-workflow_files/figure-html/plot-missingness-1.png` (bottom panel) OR `https://teammaclean.github.io/peppwR/articles/prm-genotype-power_files/figure-html/plot-missingness-1.png` | Local plot-missingness figure has three stacked panels — the bottom panel (NA rate vs log mean abundance with dashed loess) is what the slide needs. Crop or rerender just the bottom panel. The PRM vignette version may be cleaner. |
| 35 | Aggregate-mode power curve | $n$ per group vs estimated power; dashed at $0.8$; $n^\ast = 25$ marker (Figure 7 in prose) | use-existing | `vignettes/articles/power-analysis-workflow_files/figure-html/plot-sample-size-1.png` | Exact match: title reads "Recommended sample size: N=25 per group", red dashed at 80%, green marker at $n=25$. Same parameters as slide. |
| 43 | Per-peptide power distribution | Histogram of per-peptide power across 2,228 peptides at $n = 3$ (Figure 8 in prose) | adapt-existing | `https://teammaclean.github.io/peppwR/articles/dda-time-course-power_files/figure-html/plot-power-n3-1.png` OR Panel D of `paper/figures/dda_case_study.png` | DDA vignette `plot-power-n3-1.png` is closest — per-peptide power histogram at $n=3$, 2-fold, Bayes factor. Panel D of the JSS DDA figure is the same content but smaller in a six-panel composite; would need to crop. Local vignette equivalent: `vignettes/articles/power-analysis-workflow_files/figure-html/plot-pp-power-1.png` but that one is at $n=6$ — wrong $n$ for this slide. |

## pepdiff section

| Slide # | Slide title | Slot description (compressed) | Recommendation | Source | Notes |
|---|---|---|---|---|---|
| 50 | Exploratory triptych | PCA / per-sample density / per-peptide MNAR scatter (Figure 10 in prose) | use-existing | `https://teammaclean.github.io/pepdiff/articles/basic_workflow_files/figure-html/plot-data-1.png` | Output of `plot(<pepdiff_data>)` — the documented triptych. Exact match to slide description. Alternative if higher-resolution panels are needed: combine `diagnostic_plots_files/figure-html/pca-good-1.png` + `dist-good-1.png` + `miss-good-1.png` from the diagnostic-plots article. |
| 56 | When the GLM does not fit: diagnostics | Four-panel `plot_fit_diagnostics()` output: deviance histogram, deviance vs log2FC, residuals-vs-fitted for three peptides, pooled QQ (Figure 11 in prose) | use-existing | `https://teammaclean.github.io/pepdiff/articles/checking_fit_files/figure-html/show-bad-plot-1.png` | The "bad" variant is the one that motivates the diagnostic — heavier tails, obvious deviance flags. The "good" version (`show-good-plot-1.png`) is also available if Dan prefers the clean case. Both are the canonical `plot_fit_diagnostics()` output. |
| 62 | Reading the output | Three-panel `plot(<pepdiff_results>)`: volcano, p-value histogram, log2FC density (Figure 12 in prose) | use-existing | `https://teammaclean.github.io/pepdiff/articles/basic_workflow_files/figure-html/plot-results-1.png` | Direct `plot(<pepdiff_results>)` output; matches slide description exactly. Alternative split panels: `diagnostic_plots_files/figure-html/volcano-good-1.png`, `pval-good-1.png`, `fc-good-1.png` if higher resolution per panel is needed. |

## Closing

No figure slots in slides 65–67 (Takeaways, One sentence each, Questions). Skipped.

---

## Summary

### Counts

- Slides with a figure slot reviewed: **14**
- `use-existing`: **6** (slides 21, 27, 35, 50, 56, 62)
- `adapt-existing`: **2** (slides 19, 43) — note slide 19 is a hybrid (one existing + two fresh)
- `draw-fresh`: **6** (slides 6, 7, 9, 12, 14, 15, 18) — seven entries, six rows
- `skip`: **0**

(Slide 19 counts once in the `adapt-existing` total but pulls in extra `draw-fresh` work for its other two thumbnails.)

### Figures we already have locally

- `paper/figures/dda_case_study.png` — six-panel JSS case study figure; Panel D (per-peptide power histogram at $n=3$, 2-fold) is the closest local analog of slide 43.
- `paper/figures/prm_case_study.png` — six-panel JSS case study figure including MNAR detection (Panel A) and missingness-on/off comparison (Panel C).
- `vignettes/articles/power-analysis-workflow_files/figure-html/plot-sample-size-1.png` — aggregate-mode power curve with $N=25$ marker (slide 35, exact match).
- `vignettes/articles/power-analysis-workflow_files/figure-html/plot-missingness-1.png` — three-panel missingness diagnostic; bottom panel is the slide-27 candidate.
- `vignettes/articles/power-analysis-workflow_files/figure-html/density-overlay-1.png` — small-multiples observed-data densities; partial candidate for slide-19 right-skewed thumbnail.
- `vignettes/articles/power-analysis-workflow_files/figure-html/plot-pp-power-1.png` — per-peptide power histogram (at $n=6$, not $n=3$); only a partial fit for slide 43.

### Remote figures worth downloading

(Targets to copy into `talks/figures/` in the next pass.)

- `https://teammaclean.github.io/peppwR/articles/dda-time-course-power_files/figure-html/eda-plot-1.png` — slide 21 candidate.
- `https://teammaclean.github.io/peppwR/articles/dda-time-course-power_files/figure-html/plot-power-n3-1.png` — slide 43, correct $n=3$ version.
- `https://teammaclean.github.io/peppwR/articles/prm-genotype-power_files/figure-html/plot-missingness-1.png` — alternative for slide 27 (cleaner PRM version of the missingness diagnostic).
- `https://teammaclean.github.io/pepdiff/articles/basic_workflow_files/figure-html/plot-data-1.png` — slide 50, the `plot(<pepdiff_data>)` triptych.
- `https://teammaclean.github.io/pepdiff/articles/basic_workflow_files/figure-html/plot-results-1.png` — slide 62, the `plot(<pepdiff_results>)` triptych.
- `https://teammaclean.github.io/pepdiff/articles/checking_fit_files/figure-html/show-bad-plot-1.png` — slide 56, `plot_fit_diagnostics()` four-panel.
- `https://teammaclean.github.io/pepdiff/articles/checking_fit_files/figure-html/show-good-plot-1.png` — slide 56 alternative if the clean case is preferred.

### Anomalies / things to verify with Dan

- Slide 21: prose says "histogram" but the closest rendered figure is a density. Either accept the density or rerender as a true histogram from the DDA pilot CSV (`sample_data/dda_data.csv`).
- Slide 18 ("What the re-analysis showed"): the source is an external anonymous paper; no in-repo figure is appropriate. Needs either a fresh re-render from raw HR-score data (if Dan has it) or a stylised mock.
- Slide 43: only the remote DDA vignette PNG matches the $n=3$ wording; the local `plot-pp-power-1.png` is at $n=6$.
- The PRM case-study composite (`paper/figures/prm_case_study.png`) covers slides 27 and 45 between them — could be cropped panel-by-panel to avoid pulling in extra remote figures.
