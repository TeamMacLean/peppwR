# JSS Conversion Tasks: Detailed Implementation Plan

**Reference Document:** `paper/submit_to_jss_overview.md`
**Target:** Transform peppwR paper from Bioinformatics application note to full JSS article
**Execution Model:** Sequential phases, self-contained agents with explicit file references

---

## Phase 1: Infrastructure Setup

### Task 1.1: JSS Template Setup and Conversion
**Agent Type:** `general-purpose` (technical setup)
**Dependencies:** None
**Estimated Time:** 2-4 hours

**Objective:** Convert current R Markdown structure to JSS-compliant LaTeX workflow

**Specific Actions:**
1. **Download JSS Templates**
   - Get `jss-article-rnw.zip` from https://www.jstatsoft.org/style
   - Extract to `paper/jss-template/` subdirectory
   - Examine `article.Rnw` template structure and requirements

2. **Create JSS Document Structure**
   - **Base Document:** Create `paper/jss-article.Rnw` using JSS template
   - **Reference Current:** Read `paper/paper.Rmd` for content structure
   - **Bibliography:** Ensure `paper/references.bib` uses JSS format with proper `\proglang{}`, `\pkg{}` markup
   - **Style Files:** Copy `jss.cls`, `jss.bst`, `jsslogo.jpg` to paper directory

3. **Initial Content Migration**
   - Map current sections from `paper.Rmd` to JSS structure
   - Convert basic content while preserving all R code chunks
   - Replace generic citations with JSS-style markup throughout

4. **Test Compilation Pipeline**
   - Verify `Sweave("jss-article.Rnw")` → `article.tex` generation
   - Test `pdfLaTeX` compilation to PDF
   - Ensure `Stangle()` produces clean replication R code
   - Resolve any LaTeX compilation errors

**Key Files to Reference:**
- `paper/paper.Rmd` (current content)
- `paper/references.bib` (bibliography)
- `paper/style_guide.txt` (JSS requirements)

**Quality Criteria:**
- Document compiles without errors via `pdfLaTeX`
- All JSS style elements properly implemented (`\proglang{R}`, `\pkg{peppwR}`, etc.)
- Bibliography renders correctly with JSS formatting
- R code chunks execute and integrate properly

**Completion Check:**
- [ ] JSS template successfully integrated
- [ ] Basic content migrated from current paper
- [ ] Compilation pipeline works end-to-end
- [ ] Generated PDF matches JSS formatting standards

---

## Phase 2: Content Development - Core Sections

### Task 2.1: Introduction Section Expansion
**Agent Type:** `general-purpose` (research/writing)
**Dependencies:** Task 1.1 complete
**Estimated Time:** 1-2 days

**Objective:** Transform current 400-word introduction into comprehensive 800-1200 word JSS-standard introduction

**Current Content Source:** `paper/paper.Rmd` lines 23-31 (Introduction section)

**Required Expansion Areas:**
1. **Statistical Power Foundation** (200-300 words)
   - Mathematical definition of statistical power
   - Importance in experimental design across fields
   - Historical context and development of power analysis

2. **Proteomics-Specific Challenges** (300-400 words)
   - **Heterogeneity:** Detailed explanation with citations
   - **Non-normality:** Distribution characteristics with examples
   - **Missing Data:** MNAR patterns and detection limits
   - **Multiple Testing:** FDR considerations in high-throughput data

3. **Existing Tool Analysis** (200-300 words)
   - **Perseus:** Comprehensive analysis platform lacking power analysis
   - **clippda:** SELDI-TOF focus, homogeneous variance assumptions
   - **MultiPower:** Multi-omics integration requirements
   - **G*Power:** Generic approach, normality assumptions
   - Create comparison table highlighting limitations

4. **peppwR Positioning** (100-200 words)
   - Novel contributions clearly stated
   - Advantages over existing approaches
   - Target user community

**Content Sources to Reference:**
- Current introduction in `paper/paper.Rmd`
- Background material in `vignettes/getting-started.Rmd`
- Literature review from existing `paper/references.bib`

**JSS Quality Standards:**
- All software tools mentioned use `\pkg{toolname}` markup
- Programming languages use `\proglang{R}` markup
- Proper citation style throughout (`\citet{}`, `\citep{}`)
- Mathematical notation follows JSS conventions
- Section headers in sentence style

**Completion Criteria:**
- [ ] Word count: 800-1200 words
- [ ] All JSS markup properly applied
- [ ] Comprehensive tool comparison included
- [ ] Mathematical foundation established
- [ ] Clear positioning of peppwR contributions
- [ ] All citations properly formatted

### Task 2.2: Methods and Implementation Section Development
**Agent Type:** `general-purpose` (technical writing)
**Dependencies:** Task 2.1 complete
**Estimated Time:** 2-3 days

**Objective:** Expand current 600-word implementation overview into comprehensive 1500-2000 word methods section

**Current Content Source:** `paper/paper.Rmd` lines 33-68 (Implementation section)

**Required Subsections:**

1. **2.1 Statistical Framework** (300-400 words)
   - Mathematical foundation of per-peptide power analysis
   - Definition of effect size in multiplicative context
   - Power calculation formulas for non-parametric tests
   - Multiple testing considerations

2. **2.2 Distribution Fitting** (300-400 words)
   - **Reference:** `vignettes/find-fits.Rmd` for detailed methodology
   - Maximum likelihood estimation approach
   - AIC-based model selection criteria
   - Candidate distributions: gamma, lognormal, normal, inverse Gaussian
   - Fallback strategies: exclude, empirical, moment-matched lognormal

3. **2.3 Simulation Engine** (300-400 words)
   - Monte Carlo methodology
   - Convergence criteria and simulation count selection
   - Effect size application (multiplicative scaling)
   - Parallelization opportunities

4. **2.4 Missing Data Handling** (200-300 words)
   - Per-peptide missingness calculation
   - MNAR detection via abundance-missingness correlation
   - Missingness-aware simulation methodology
   - Integration with detection limit concepts

5. **2.5 Software Architecture** (200-300 words)
   - S3 class system: `peppwr_fits`, `peppwr_power`
   - Method dispatch for print, plot, summary
   - Extensibility design for new distributions/tests
   - Integration with existing R ecosystem

**Content Sources to Reference:**
- `R/fits.R` - Distribution fitting implementation
- `R/power.R` - Power analysis engine
- `R/simulation.R` - Monte Carlo implementation
- `R/classes.R` - S3 class definitions
- `vignettes/find-fits.Rmd` - Distribution fitting details

**Quality Standards:**
- All function names use `\code{function_name()}` markup
- Mathematical notation properly formatted
- Algorithm descriptions with step-by-step clarity
- Code examples integrated throughout
- Technical accuracy verified against implementation

**Completion Criteria:**
- [ ] Word count: 1500-2000 words
- [ ] All 5 subsections completed
- [ ] Mathematical framework clearly established
- [ ] Implementation details match actual code
- [ ] JSS formatting standards applied throughout

### Task 2.3: Examples and Applications Section Enhancement
**Agent Type:** `general-purpose` (analysis/writing)
**Dependencies:** Task 2.2 complete
**Estimated Time:** 2-3 days

**Objective:** Transform current 450-word example into comprehensive 1000-1500 word applications section

**Current Content Source:** `paper/paper.Rmd` lines 70-88 (Example Application section)

**Required Subsections:**

1. **3.1 Simulated Data Validation** (300-400 words)
   - Ground truth comparison with known parameters
   - Validation of power estimates against theoretical expectations
   - Edge case testing (small samples, extreme effect sizes)
   - **Code Source:** Develop new validation study

2. **3.2 DDA Phosphoproteomics Case Study** (300-400 words)
   - **Reference:** `vignettes/articles/dda-time-course-power.Rmd`
   - Real pilot data from published study
   - Complete workflow: data → fitting → power analysis
   - Interpretation of per-peptide heterogeneity

3. **3.3 PRM Time-Course Analysis** (200-300 words)
   - **Reference:** `vignettes/articles/prm-genotype-power.Rmd`
   - Alternative experimental design scenario
   - Comparison with DDA approach
   - Statistical test selection rationale

4. **3.4 Integration Examples** (200-300 words)
   - Workflow with existing proteomics tools
   - Data import from common formats
   - Export of results for downstream analysis
   - Connection to pepdiff package for complete pipeline

**Content Sources to Reference:**
- `vignettes/articles/dda-time-course-power.Rmd` - Real data example
- `vignettes/articles/prm-genotype-power.Rmd` - Alternative case
- `sample_data/dda_data.csv` - Example dataset
- `sample_data/prm_data.csv` - Alternative dataset

**Quality Standards:**
- All examples fully reproducible with provided code
- Statistical interpretations accurate and meaningful
- Dataset characteristics clearly described
- Results presented with appropriate uncertainty quantification
- Practical implications for experimental design emphasized

**Completion Criteria:**
- [ ] Word count: 1000-1500 words
- [ ] All 4 subsections completed
- [ ] Real data examples included and analyzed
- [ ] Validation study demonstrates accuracy
- [ ] Integration examples show practical workflow

---

## Phase 3: New Content Development

### Task 3.1: Benchmarking and Validation Section Creation
**Agent Type:** `general-purpose` (analysis/benchmarking)
**Dependencies:** Task 2.3 complete
**Estimated Time:** 3-4 days

**Objective:** Create entirely new 1000-word section demonstrating peppwR's advantages through systematic comparison

**Content Foundation:** `vignettes/articles/benchmarking.Rmd` (performance analysis starting point)

**Required Subsections:**

1. **4.1 Tool Comparison** (300-400 words)
   - Head-to-head comparison with available alternatives
   - Create comparison matrix table (features, assumptions, outputs)
   - Quantitative accuracy assessment where possible
   - Limitations of existing tools for proteomics data

2. **4.2 Performance Analysis** (200-300 words)
   - **Reference:** `vignettes/articles/benchmarking.Rmd`
   - Runtime scaling with peptide count and simulation iterations
   - Memory usage characteristics
   - Computational complexity analysis
   - Performance recommendations for different dataset sizes

3. **4.3 Accuracy Assessment** (300-400 words)
   - Simulation studies with known ground truth parameters
   - Power estimation accuracy across parameter ranges
   - Comparison with theoretical calculations where available
   - Edge case behavior documentation

4. **4.4 Edge Case Testing** (200-300 words)
   - Behavior with extreme parameters (very high/low variance)
   - Small sample size performance
   - High missingness scenarios
   - Convergence properties of simulations

**Development Requirements:**
- Implement comparison studies if not already available
- Generate benchmark datasets for consistent comparison
- Create new analysis scripts for accuracy validation
- Develop performance profiling code

**Quality Standards:**
- All benchmarking methodologies clearly described
- Results presented with confidence intervals
- Fair comparison conditions for all tools
- Performance recommendations actionable
- Statistical significance testing of differences

**Completion Criteria:**
- [ ] Word count: 1000 words
- [ ] Comprehensive tool comparison completed
- [ ] Performance benchmarking results included
- [ ] Accuracy validation demonstrates reliability
- [ ] Edge case behavior documented

### Task 3.2: Figure and Table Development
**Agent Type:** `general-purpose` (visualization)
**Dependencies:** Tasks 2.1-3.1 complete
**Estimated Time:** 2-3 days

**Objective:** Create comprehensive figure set supporting JSS article content

**Current Assets:** `paper/figure_1.pdf` (existing multi-panel figure)

**Required Figures:**

1. **Figure 1: Workflow and Core Results** (existing, may need updates)
   - **Source:** `paper/supplementary_example.Rmd`
   - Ensure JSS-compliant formatting and resolution
   - Update captions for JSS sentence style

2. **Figure 2: Distribution Fitting Diagnostics** (new)
   - QQ plots for different distribution candidates
   - AIC comparison across peptides
   - Fit quality visualization
   - **Code Source:** Extend `vignettes/find-fits.Rmd` examples

3. **Figure 3: Benchmarking Results** (new)
   - Performance comparison bar charts/box plots
   - Runtime scaling curves
   - Accuracy comparison with confidence intervals
   - **Code Source:** Develop from Task 3.1 benchmarking analysis

4. **Figure 4: Real Data Application** (new)
   - **Source:** `vignettes/articles/dda-time-course-power.Rmd`
   - Power analysis results for real DDA dataset
   - Per-peptide power distribution
   - Sample size recommendation visualization

5. **Figure 5: Validation Results** (new)
   - Ground truth comparison plots
   - Accuracy assessment across parameter ranges
   - Edge case behavior visualization

**Required Tables:**

1. **Table 1: Tool Comparison Matrix**
   - Feature comparison across existing tools
   - Assumptions, inputs, outputs
   - Proteomics-specific capabilities

2. **Table 2: Performance Benchmarks**
   - Runtime and memory usage across dataset sizes
   - Accuracy metrics comparison
   - Statistical significance tests

**Quality Standards:**
- All figures vector format (PDF) with embedded fonts
- JSS-compliant resolution and sizing
- Captions in sentence style ending with periods
- Color schemes accessible (colorblind-friendly)
- Clear axis labels and legends
- Tables properly formatted with JSS style

**Completion Criteria:**
- [ ] All 5 figures created and JSS-compliant
- [ ] Both tables properly formatted
- [ ] Figure captions follow JSS sentence style
- [ ] All graphics integrate properly with LaTeX compilation
- [ ] Visual quality appropriate for publication

---

## Phase 4: Final Integration and Polish

### Task 4.1: Discussion and Summary Sections
**Agent Type:** `general-purpose` (writing/synthesis)
**Dependencies:** All Phase 2-3 tasks complete
**Estimated Time:** 1-2 days

**Objective:** Create comprehensive discussion (500-800 words) and summary (200 words) sections

**Current Content Source:** `paper/paper.Rmd` lines 89-105 (Discussion section)

**Discussion Section Requirements:**
1. **Synthesis of Results** (200-300 words)
   - Integration of benchmarking and validation findings
   - Practical implications for proteomics researchers
   - Cost-benefit analysis of per-peptide vs. aggregate approaches

2. **Limitations and Assumptions** (150-200 words)
   - Computational scalability constraints
   - Pilot data requirements
   - Distributional assumptions and fallbacks

3. **Integration with Proteomics Ecosystem** (150-200 words)
   - Workflow integration examples
   - Connection to downstream analysis tools
   - Position in experimental design → analysis pipeline

4. **Future Directions** (100-150 words)
   - TMT/iTRAQ extension possibilities
   - Parallel processing implementation
   - Integration with pepdiff package
   - Community feedback incorporation

**Summary Section Requirements:**
- Concise restatement of problem and solution
- Key numerical results highlighted
- Software availability and accessibility
- Primary contributions emphasized

**Quality Standards:**
- Discussion balanced between strengths and limitations
- Future work specific and actionable
- Summary standalone readable
- All claims supported by results in earlier sections

**Completion Criteria:**
- [ ] Discussion: 500-800 words
- [ ] Summary: 200 words
- [ ] Balanced analysis of strengths/limitations
- [ ] Clear future work roadmap
- [ ] Integration with broader ecosystem discussed

### Task 4.2: Final Formatting and JSS Compliance
**Agent Type:** `general-purpose` (editorial/formatting)
**Dependencies:** All previous tasks complete
**Estimated Time:** 1-2 days

**Objective:** Ensure complete JSS compliance and publication readiness

**JSS Style Checklist (from `paper/style_guide.txt`):**
- [ ] Document compiles with `pdfLaTeX`
- [ ] `\proglang{}`, `\pkg{}`, `\code{}` used throughout (including titles/references)
- [ ] References in `.bib` format with proper JSS markup
- [ ] Title in title style, sections in sentence style
- [ ] All titles in BibTeX file in title style
- [ ] Figures/tables marked with `\label{}` and referenced with `\ref{}`
- [ ] Software packages properly cited

**Additional Quality Checks:**
- Bibliography completeness and accuracy
- Figure/table numbering and referencing
- Code chunk integration and execution
- Replication materials completeness
- Abstract accuracy and completeness

**Final Tasks:**
1. **Complete Style Review**
   - Systematic check against JSS style guide
   - Fix any remaining formatting issues
   - Verify all markup applications

2. **Replication Materials**
   - Generate final `.R` file via `Stangle()`
   - Test all code execution independently
   - Ensure data availability and accessibility

3. **Submission Package Preparation**
   - Organize all required files
   - Create submission-ready archive
   - Verify compilation on clean system

**Completion Criteria:**
- [ ] All JSS style requirements met
- [ ] Document compiles without errors
- [ ] Replication materials complete and tested
- [ ] Submission package ready for JSS portal
- [ ] Final PDF matches JSS formatting standards

---

## Task Execution Notes

### General Guidelines for All Agents:
1. **File References:** Always read specified source files before beginning work
2. **Quality First:** JSS has rigorous standards; accuracy and completeness essential
3. **Iterative Refinement:** Initial drafts acceptable; plan for revision cycles
4. **Documentation:** Comment code thoroughly for replication
5. **Integration:** Ensure new content integrates seamlessly with existing sections

### Success Metrics:
- Article length: 4000-6000 words total
- 5-6 figures, 2-3 tables
- Comprehensive references (40+ citations)
- Complete replication materials
- JSS formatting compliance

### Risk Mitigation:
- Test compilation frequently during development
- Validate all numerical results independently
- Cross-reference claims with implementation
- Maintain backup copies of working versions

---

**Total Estimated Timeline:** 2-3 weeks with sequential execution
**Target Completion:** Ready for JSS submission with publication-quality materials