# peppwR Paper

Bioinformatics Application Note describing the peppwR package.

## Target

- **Journal**: Bioinformatics (Oxford) - Application Note
- **Preprint**: bioRxiv (same content, different template)
- **Format**: ~2,000 words + 1 figure

## Files

| File | Purpose |
|------|---------|
| `outline.md` | Full paper outline with rendering pipeline |
| `references.bib` | BibTeX references (18 entries) |
| `example_application_spec.md` | Detailed spec for supplementary code |
| `supplementary_example.Rmd` | Reproducible code (to be created) |
| `paper.Rmd` | Main manuscript (to be created) |

## Key Messages

1. **Per-peptide heterogeneity matters**: Generic tools underestimate sample size by ~60%
2. **MNAR-aware simulations**: Accounting for missingness requires more samples
3. **Test choice affects power**: Bayes factor ~15-20% more efficient than Wilcoxon

## Dependencies

```r
# For rendering (add to renv, not package deps)
renv::install("rticles")
```

## Rendering

```r
# Preprint (bioRxiv) - no proprietary footer
rmarkdown::render("paper.Rmd", output_format = rticles::arxiv_article())

# Journal (Bioinformatics) - OUP template
rmarkdown::render("paper.Rmd", output_format = rticles::bioinformatics_article())
```

## Status

- [x] Outline complete
- [x] References collected
- [x] Example application spec complete
- [ ] supplementary_example.Rmd
- [ ] paper.Rmd
- [ ] Figure 1
- [ ] Preprint submission
- [ ] Journal submission
