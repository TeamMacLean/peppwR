# Next Session: Re-render Examples

## Task

Re-render the DDA and PRM example documents after recent updates.

## What Was Done

1.  Fixed `find = "effect_size"` for per-peptide mode (was not working
    because installed package was out of date)
2.  Updated DDA example with data-derived heatmap parameters
3.  Added explanation of the two thresholds (`target_power` vs
    `proportion_threshold`) to both examples
4.  Updated CLAUDE.md with project status and new bug (misleading red
    line)

## What Needs To Be Done

### 1. Re-render Examples

The Rmd files were updated but HTML not yet re-rendered:

``` r
rmarkdown::render("inst/examples/dda-time-course-power.Rmd")
rmarkdown::render("inst/examples/prm-genotype-power.Rmd")
```

Note: Rendering takes ~20-25 minutes per document due to power
simulations.

### 2. Verify Plots

Check that the effect_size plots now show meaningful curves (not empty).

### 3. Commit

``` bash
git add inst/examples/*.Rmd inst/examples/*.html
git commit -m "Re-render examples with threshold explanation"
```

## Known Bug to Fix (Future Session)

**Misleading red line in per-peptide effect_size plot**

The plot draws a red horizontal line at `proportion_threshold` (50%) but
users confuse it with `target_power` (80%).

Fix in `R/plots.R`: - Remove the red line - Replace with text annotation
explaining the answer

See CLAUDE.md “Known Issues” section for details.

## Files Modified (Not Yet Committed)

- `inst/examples/dda-time-course-power.Rmd` - Added threshold
  explanation
- `inst/examples/prm-genotype-power.Rmd` - Added threshold explanation
- `CLAUDE.md` - Updated project status, added new bug
