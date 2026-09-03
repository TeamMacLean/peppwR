# Task 14: Integrate DDA Case Study with JSS Figure Formatting

## PROJECT CONTEXT

**Objective**: Convert peppwR R package into publication-ready Journal of Statistical Software (JSS) article

**Current Problem**: Section 5.1 (DDA Case Study) contains fabricated content including an inappropriate "Targeted vs Discovery Proteomics" comparison that violates scientific integrity and editorial standards.

**Scientific Integrity Requirement**: Use only real examples from vignettes. No synthetic data, fabricated results, or inappropriate comparisons in case studies.

**Editorial Standards**: JSS articles require focused, academically rigorous content. Avoid tutorial-style comparisons that don't add scientific value.

## TASK OVERVIEW

Replace the current fabricated Section 5.1 in `jss-article.Rnw` with authentic content from `section_5_1_dda_case_study.rnw`, ensuring proper JSS publication formatting and removing inappropriate comparison content.

## CURRENT STATE ANALYSIS

**File**: `/Users/macleand/Desktop/peppwR/paper/jss-article.Rnw`
**Problem Section**: Section 5.1 - DDA Case Study (approximately lines 608-750)
**Issues with Current Content**:
- Contains fabricated examples not from established vignettes
- Includes inappropriate "Targeted vs Discovery Proteomics: Power Implications" comparison section
- Inconsistent missing rate claims and numbers
- Violates scientific integrity standards for academic publication

**Source File**: `/Users/macleand/Desktop/peppwR/paper/section_5_1_dda_case_study.rnw`
**Source Content**: Authentic examples extracted from `vignettes/articles/dda-time-course-power.Rmd`

## CONTENT INCLUSION/EXCLUSION CRITERIA

### ✅ CONTENT TO INCLUDE (from source vignette):

**Essential DDA Case Study Elements**:
- Real data loading: `../sample_data/dda_data.csv`
- Data filtering: Early (t=0) vs late (t=600) timepoints
- Authentic dataset statistics: 2,228 phosphopeptides, 3 replicates per condition
- Distribution fitting with small sample caveats
- Statistical test comparison (Wilcoxon, Bootstrap-t, Bayes factor)
- Three core questions: current power, sample size determination, minimum effect
- Power analysis results and interpretations

**Scientific Narrative**:
- Small sample size limitations and their impact
- Distribution fitting reliability warnings for N=6
- Authentic experimental design challenges
- Practical power analysis workflow demonstration

### ❌ CONTENT TO EXCLUDE:

**Inappropriate Comparisons**:
- Any comparison between "targeted vs discovery" approaches
- Tutorial-style mode comparisons that don't fit JSS academic tone
- Redundant explanations of analysis differences
- Claims about "superior power characteristics" between experimental approaches

**Editorial Inappropriate Content**:
- Overly tutorial-style explanations
- Content that reads like package documentation rather than academic case study
- Redundant methodology explanations already covered in Methods section
- Comparative claims without rigorous statistical backing

**Fabricated Content**:
- Any content not directly derived from established vignettes
- Synthetic discussion points or artificial conclusions
- Made-up comparative advantages or disadvantages

## DETAILED REQUIREMENTS

### 1. FIGURE FORMATTING REQUIREMENTS

**Required Changes**:
- Chunk name: `<<dda-power-plot, fig=TRUE, include=FALSE, width=6, height=4>>=`
  - Generates: `jss-article-dda-power-plot.pdf`

**LaTeX Integration**:
- Add label: `\label{fig:dda-power}`
- Reference in text: `Figure~\ref{fig:dda-power}`
- Caption: JSS academic style describing DDA power analysis results

### 2. DATA INTEGRITY REQUIREMENTS

**Must Preserve**:
- Real data path: `../sample_data/dda_data.csv`
- Authentic peptide count: 2,228 phosphopeptides
- Real experimental design: 3 replicates per condition, t=0 vs t=600
- Small sample size caveats and distribution fitting warnings
- Established narrative about experimental constraints

**Must Remove**:
- Any inconsistent missing rate claims
- Fabricated dataset statistics
- Synthetic experimental descriptions

### 3. ACADEMIC TONE REQUIREMENTS

**JSS Publication Standards**:
- Focus on demonstrating software capability with real data
- Emphasize practical power analysis workflow
- Maintain academic rigor without tutorial-style explanations
- Present results objectively without overstated claims

## STEP-BY-STEP EXECUTION PLAN

### Step 1: Analyse Current State
1. Open `/Users/macleand/Desktop/peppwR/paper/jss-article.Rnw`
2. Locate Section 5.1 boundaries (search for DDA case study)
3. Note exact line numbers for replacement
4. Identify and document inappropriate comparison content to remove

### Step 2: Prepare Source Content
1. Open `/Users/macleand/Desktop/peppwR/paper/section_5_1_dda_case_study.rnw`
2. **CRITICALLY REVIEW** content against exclusion criteria
3. Remove any tutorial-style comparisons or inappropriate content
4. Modify figure chunks to JSS requirements:
   - Add `width=6, height=4` parameters
   - Ensure proper chunk naming for PDF generation
5. Add LaTeX figure labels and references

### Step 3: Content Curation
1. **EXCLUDE** any content that compares experimental approaches
2. **EXCLUDE** tutorial-style explanations inappropriate for JSS
3. **FOCUS ON** demonstrating peppwR capabilities with real data
4. **MAINTAIN** authentic vignette results and scientific narrative

### Step 4: Content Integration
1. Replace identified section in main file with curated source content
2. Verify subsection numbering: `\subsection{DDA Phosphoproteomics Time-Course Analysis}`
3. Check package reference formatting (`\pkg{peppwR}`)
4. Ensure proper LaTeX structure and academic tone

### Step 5: Compilation Testing
1. Run: `R -e "Sweave('jss-article.Rnw')"`
2. Verify figure PDF generates with correct name: `jss-article-dda-power-plot.pdf`
3. Check LaTeX compilation: `pdflatex jss-article.tex`
4. Verify figure references resolve correctly

### Step 6: Quality Verification
1. Confirm 2,228 peptides and authentic experimental parameters appear
2. Verify no fabricated content or inappropriate comparisons remain
3. Check that content maintains JSS academic publication standards
4. Compare against original vignette to ensure authenticity maintained

## SUCCESS CRITERIA

- [ ] Fabricated Section 5.1 completely replaced
- [ ] Authentic DDA vignette content integrated
- [ ] Inappropriate "targeted vs discovery" comparison removed
- [ ] One figure generates with proper JSS naming
- [ ] Figure references resolve in compiled PDF
- [ ] 2,228 peptides and authentic parameters preserved
- [ ] Small sample caveats and scientific limitations maintained
- [ ] Academic tone appropriate for JSS publication
- [ ] No tutorial-style comparisons included
- [ ] Clean compilation with no errors
- [ ] Scientific integrity maintained

## ERROR HANDLING

**If content appears inappropriate**: Re-review against exclusion criteria
**If Sweave fails**: Check R code syntax and data paths
**If figures don't generate**: Verify chunk parameters and naming
**If LaTeX fails**: Check figure references and labels
**If tone seems too tutorial-like**: Revise to academic publication standards

## DEPENDENCIES

- Requires `section_5_1_dda_case_study.rnw` to exist (completed by previous agents)
- Must be completed after Task 13 (per-peptide integration)
- Blocks Task 16 (final compilation)

## DELIVERABLES

1. Modified `jss-article.Rnw` with authentic, academically appropriate Section 5.1
2. Generated DDA figure PDF with JSS naming convention
3. Verification that inappropriate comparison content removed
4. Documentation confirming academic publication standards met