# Fabricated Citations Cleanup Summary

**Date:** 2026-03-27
**Objective:** Remove 5 confirmed fabricated citations and replace with [CITATION NEEDED] placeholders

## Citations Removed from references.bib

1. **@article{wilkinson2019statistical}** - "Statistical methods in proteomics" - DOI fails, fabricated
2. **@article{hultin2005abundance}** - "Abundance and dynamics of proteins in the human plasma proteome" - DOI 404, doesn't exist
3. **@article{cox2011principles}** - "Principles and practice of quantitative mass spectrometry..." - DOI 404, doesn't exist
4. **@article{o2018accounting}** - "Accounting for missing values in multiplexed quantitative proteomics" - DOI 404, doesn't exist
5. **@article{varemo2013modeling}** - "Modeling the effect of cellular heterogeneity..." - DOI resolves to wrong article

## Text Replacements in jss-article.tex

### Line 74
**Before:** `\citep{wilkinson2019statistical}`
**After:** `[CITATION NEEDED]`
**Context:** Power analysis neglected in proteomics research surveys

### Line 77 (variance structure)
**Before:** `\citep{karp2010addressing, varemo2013modeling}`
**After:** `\citep{karp2010addressing} [CITATION NEEDED]`
**Context:** Dramatic differences in variance structure across peptidome

### Line 77 (distribution types)
**Before:** `\citep{hultin2005abundance, webb2015mnar}`
**After:** `[CITATION NEEDED] \citep{webb2015mnar}`
**Context:** Gamma, lognormal, inverse Gaussian distributions

### Line 77 (MS quantification)
**Before:** `\citep{cox2011principles}`
**After:** `[CITATION NEEDED]`
**Context:** Multiplicative nature of mass spectrometry quantification

### Line 79 (MCAR assumptions)
**Before:** `\citep{tyanova2016imputation, o2018accounting}`
**After:** `\citep{tyanova2016imputation} [CITATION NEEDED]`
**Context:** Missing-completely-at-random assumptions violations

### Line 86 (reproducibility challenges)
**Before:** `\citep{wilkinson2019statistical, collins2020multi}`
**After:** `[CITATION NEEDED] \citep{collins2020multi}`
**Context:** Reproducibility challenges in proteomics research

## Verification

- ✅ All 5 fabricated entries completely removed from references.bib
- ✅ All in-text citations replaced with [CITATION NEEDED] placeholders
- ✅ Remaining valid citations preserved intact
- ✅ Document structure maintained (no syntax errors detected)
- ✅ Total of 6 [CITATION NEEDED] placeholders added

## Next Steps

During revision phase, these [CITATION NEEDED] placeholders should be replaced with legitimate citations that support the claims being made.

**Status:** COMPLETE - Document is clean of fabricated citations and ready for continued content review.