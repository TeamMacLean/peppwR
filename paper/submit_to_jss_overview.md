# JSS Submission Plan: peppwR Paper Conversion

**Date Created:** 2026-03-25
**Current Status:** Planning Phase
**Target Journal:** Journal of Statistical Software (JSS)
**Estimated Timeline:** 2-3 weeks total effort

## Transformation Overview

### Current State
- **Format:** Bioinformatics Application Note (~2,000 words, 1 figure)
- **Structure:** R Markdown → Generic LaTeX
- **Content:** Brief software announcement format
- **Files:** `paper/paper.Rmd`, `paper/paper.pdf`, `paper/references.bib`

### JSS Target State
- **Format:** Full JSS research article (15-25 pages, 4-6 figures, 2-3 tables)
- **Structure:** JSS LaTeX class with specific markup requirements
- **Content:** Comprehensive software description with validation and benchmarking
- **Files:** JSS-compliant `.Rnw` or `.tex` with complete replication materials

## Key Transformation Requirements

### 1. Technical Infrastructure Changes

**JSS Template Setup:**
- Download `jss-article-rnw.zip` from https://www.jstatsoft.org/style
- Convert from generic `\documentclass{article}` to `\documentclass[article]{jss}`
- Implement JSS-specific markup: `\proglang{R}`, `\pkg{peppwR}`, `\code{function()}`
- Update bibliography with JSS style requirements

**Compilation Pipeline:**
- Switch from R Markdown to Sweave/knitr workflow
- Generate both LaTeX source and replication R code
- Ensure all figures compile with proper JSS formatting

### 2. Content Expansion Requirements

| Section | Current Words | JSS Target | Expansion Factor |
|---------|---------------|------------|------------------|
| Introduction | 400 | 800-1200 | 2-3x |
| Implementation | 600 | 1500-2000 | 2.5-3x |
| Examples | 450 | 1000-1500 | 2-3x |
| Discussion | 200 | 500-800 | 2.5-4x |
| **TOTAL** | **1,650** | **4,000-6,000** | **2.5-3.5x** |

### 3. New Content Requirements

**Additional Sections Needed:**
- Comprehensive literature review and positioning
- Detailed mathematical framework
- Multiple validation studies
- Benchmarking against existing tools
- Performance analysis (runtime, memory, accuracy)
- Complete algorithm descriptions
- Software architecture details

**Additional Figures/Tables:**
- Figure 2: Distribution fitting diagnostics
- Figure 3: Benchmarking results
- Figure 4: Real data application (DDA)
- Figure 5: Validation with ground truth
- Table 1: Tool comparison matrix
- Table 2: Performance benchmarks

## Detailed Section Plan

### 1. Introduction (800-1200 words)
**Current:** Brief motivation and gap identification
**JSS Needs:**
- Comprehensive review of power analysis in statistics
- Detailed coverage of proteomics-specific challenges
- Systematic comparison with existing tools (Perseus, clippda, MultiPower, G*Power)
- Clear positioning of peppwR's novel contributions
- Mathematical foundation of the problem

### 2. Methods and Implementation (1500-2000 words)
**Current:** High-level workflow description
**JSS Needs:**
- **2.1 Statistical Framework:** Mathematical foundation of per-peptide power analysis
- **2.2 Distribution Fitting:** Algorithm details, AIC selection, fallback strategies
- **2.3 Simulation Engine:** Monte Carlo methodology, convergence criteria
- **2.4 Missing Data Handling:** MNAR detection, missingness-aware simulations
- **2.5 Software Architecture:** S3 classes, method dispatch, extensibility

### 3. Examples and Applications (1000-1500 words)
**Current:** Single simulated dataset
**JSS Needs:**
- **3.1 Simulated Data Validation:** Ground truth comparison with known parameters
- **3.2 DDA Phosphoproteomics:** Real-world case study with pilot data
- **3.3 PRM Time-Course Analysis:** Alternative experimental design scenario
- **3.4 Integration Examples:** Workflow with existing proteomics tools

### 4. Benchmarking and Validation (1000 words) *[NEW SECTION]*
**JSS Requirements:**
- **4.1 Tool Comparison:** Head-to-head with existing solutions
- **4.2 Performance Analysis:** Runtime scaling, memory usage
- **4.3 Accuracy Assessment:** Simulation studies with known ground truth
- **4.4 Edge Case Testing:** Behavior with extreme parameters, small samples

### 5. Discussion and Future Work (500-800 words)
**Current:** Brief limitations and future directions
**JSS Needs:**
- Detailed analysis of strengths and limitations
- Computational complexity discussion
- Integration with broader proteomics ecosystem
- Specific roadmap for extensions (TMT, iTRAQ, etc.)

### 6. Summary (200 words) *[NEW SECTION]*
**JSS Standard:** Concise summary of contributions and availability

## Content Sources Available

**Existing Materials to Leverage:**
- `vignettes/getting-started.Rmd` - Basic usage examples
- `vignettes/find-fits.Rmd` - Distribution fitting details
- `vignettes/articles/benchmarking.Rmd` - Performance analysis foundation
- `vignettes/articles/dda-time-course-power.Rmd` - Real data example
- `vignettes/articles/prm-genotype-power.Rmd` - Alternative use case
- Package documentation and function examples

**Additional Development Needed:**
- Comprehensive tool comparison study
- Ground truth validation with simulated data
- Mathematical framework formalization
- Performance benchmarking implementation

## Success Criteria

**JSS Reviewer Expectations:**
1. **Novelty:** Clear advancement over existing power analysis tools
2. **Rigor:** Comprehensive validation with multiple datasets and scenarios
3. **Usability:** Extensive examples showing practical value to researchers
4. **Integration:** Demonstration of fit within existing proteomics workflows
5. **Reproducibility:** Complete replication materials and data availability
6. **Technical Quality:** Proper JSS formatting, markup, and style compliance

## Risk Factors

**High Risk:**
- Content expansion may reveal gaps requiring additional research/development
- JSS review process is rigorous; papers commonly undergo major revisions
- Benchmarking may require implementing comparison methods

**Medium Risk:**
- JSS LaTeX formatting can be finicky; template conversion may require iteration
- Figure quality requirements may necessitate redevelopment of existing plots

**Mitigation Strategies:**
- Start with infrastructure setup to identify technical issues early
- Leverage existing vignette content to minimize new development
- Plan for iterative submission process with JSS

## Next Steps

1. **Infrastructure Setup:** Convert to JSS template and test compilation
2. **Content Audit:** Map existing vignette content to JSS sections
3. **Task Decomposition:** Break down into specific, actionable work items
4. **Implementation:** Execute transformation in phases
5. **Review & Polish:** Final formatting and submission preparation

---

**Note:** This overview provides the strategic framework. Specific implementation tasks will be detailed in `conversion_tasks.md` for execution by specialized agents.