# Task 13: Integrate Per-Peptide Section with JSS Figure Formatting

## PROJECT CONTEXT

**Objective**: Convert peppwR R package into publication-ready Journal of Statistical Software (JSS) article

**Current Problem**: Section 4 (Per-Peptide Mode Analysis) contains fabricated content that violates scientific integrity. Must be replaced with authentic content extracted from established package vignettes.

**Scientific Integrity Requirement**: Use only real examples from vignettes. No synthetic data or fabricated results in case studies.

## TASK OVERVIEW

Replace the current fabricated Section 4 in `jss-article.Rnw` with authentic content from `section_4_per_peptide_mode.rnw`, ensuring proper JSS publication formatting and figure generation.

## CURRENT STATE ANALYSIS

**File**: `/Users/macleand/Desktop/peppwR/paper/jss-article.Rnw`
**Problem Section**: Lines 423-551 (approximately) - Section 4: Per-Peptide Mode Analysis
**Issues with Current Content**:
- Contains fabricated examples not from established vignettes
- Inconsistent with package documentation
- Violates scientific integrity standards for academic publication

**Source File**: `/Users/macleand/Desktop/peppwR/paper/section_4_per_peptide_mode.rnw`
**Source Content**: Authentic examples extracted from `vignettes/articles/power-analysis-workflow.Rmd`

## DETAILED REQUIREMENTS

### 1. FIGURE FORMATTING REQUIREMENTS

JSS requires specific figure naming and formatting:

**Current Issues**: Figures in source file likely use generic chunk names
**Required Changes**:
- Chunk name: `<<per-peptide-power-plot, fig=TRUE, include=FALSE, width=6, height=4>>=`
  - Generates: `jss-article-per-peptide-power-plot.pdf`
- Chunk name: `<<per-peptide-sample-size-plot, fig=TRUE, include=FALSE, width=6, height=4>>=`
  - Generates: `jss-article-per-peptide-sample-size-plot.pdf`

**LaTeX Integration**:
- Add labels: `\label{fig:per-peptide-power}` and `\label{fig:per-peptide-sample-size}`
- References in text: `Figure~\ref{fig:per-peptide-power}`
- Captions: Use JSS academic style with detailed descriptions

### 2. CONTENT REQUIREMENTS

**Must Preserve from Source**:
- 100 peptides pilot data generation (`set.seed(42)`)
- N=6 per group analysis
- 73% of peptides achieve 80% power result
- Authentic vignette code and parameters
- Scientific narrative and explanations

**JSS Formatting**:
- Package references: `\pkg{peppwR}` format
- Code formatting: Proper Sweave syntax
- Section header: `\section{Per-Peptide Mode Analysis} \label{sec:per-peptide}`

## STEP-BY-STEP EXECUTION PLAN

### Step 1: Analyse Current State
1. Open `/Users/macleand/Desktop/peppwR/paper/jss-article.Rnw`
2. Locate Section 4 boundaries (search for `\section{Per-Peptide Mode Analysis}`)
3. Note exact line numbers for replacement
4. Document what fabricated content exists

### Step 2: Prepare Source Content
1. Open `/Users/macleand/Desktop/peppwR/paper/section_4_per_peptide_mode.rnw`
2. Identify figure chunks that need JSS formatting
3. Modify chunk headers to JSS requirements:
   - Add `width=6, height=4` parameters
   - Ensure proper chunk naming for PDF generation
4. Add LaTeX figure labels and references

### Step 3: Content Integration
1. Replace identified section in main file with modified source content
2. Verify section numbering consistency
3. Check package reference formatting (`\pkg{peppwR}`)
4. Ensure proper LaTeX structure

### Step 4: Compilation Testing
1. Run: `R -e "Sweave('jss-article.Rnw')"`
2. Verify figure PDFs generate with correct names:
   - `jss-article-per-peptide-power-plot.pdf`
   - `jss-article-per-peptide-sample-size-plot.pdf`
3. Check LaTeX compilation: `pdflatex jss-article.tex`
4. Verify figure references resolve correctly

### Step 5: Quality Verification
1. Confirm 100 peptides, N=6, 73% power results appear
2. Verify no fabricated content remains
3. Check figure quality and JSS formatting compliance
4. Compare against original vignette to ensure authenticity

## SUCCESS CRITERIA

- [ ] Fabricated Section 4 completely replaced
- [ ] Authentic vignette content integrated
- [ ] Two figures generate with proper JSS naming
- [ ] Figure references resolve in compiled PDF
- [ ] 73% power result preserved from vignette
- [ ] Clean compilation with no errors
- [ ] Scientific integrity maintained

## ERROR HANDLING

**If Sweave fails**: Check R code syntax and package loading
**If figures don't generate**: Verify chunk parameters and naming
**If LaTeX fails**: Check figure references and labels
**If results don't match vignette**: Re-verify source content authenticity

## DEPENDENCIES

- Requires `section_4_per_peptide_mode.rnw` to exist (completed by previous agents)
- Must be completed before Task 14 (DDA integration)
- Blocks Task 16 (final compilation)

## DELIVERABLES

1. Modified `jss-article.Rnw` with authentic Section 4
2. Generated figure PDFs with JSS naming convention
3. Verification that results match original vignette
4. Documentation of changes made