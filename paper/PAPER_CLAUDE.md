# peppwR JSS Paper Status

## Project Aims
Convert the peppwR R package into a publication-ready Journal of Statistical Software (JSS) article. The package provides simulation-based power analysis specifically for phosphoproteomics data, addressing three fundamental experimental design questions:

1. **"What sample size do I need?"** (`find = "sample_size"`)
2. **"What power do I have?"** (`find = "power"`)
3. **"What's the minimum detectable effect?"** (`find = "effect_size"`)

## Current Status: PARTIALLY COMPLETE

### ✅ **Completed Sections**

#### **Introduction Section** (Complete)
- Statistical power foundation and proteomics challenges
- Comprehensive tool comparison table
- peppwR positioning and contribution

#### **Methods and Implementation Section** (Complete)
- Statistical framework with mathematical formulations
- Distribution fitting methodology
- Simulation engine details
- Missing data handling (MNAR detection)
- Software architecture (S3 classes)

#### **Aggregate Mode Examples Section** (Complete)
- **Sample size determination**: Gamma distribution example with power curve figure
- **Power estimation**: Fixed sample size scenario
- **Effect size determination**: Minimum detectable effect with figure
- **Parameter sensitivity analysis**: Varying distribution assumptions
- All three core questions answered with executable R code

#### **References Section** (Complete)
- Working bibliography with 16 citations
- Proper JSS citation formatting

#### **Technical Infrastructure** (Complete)
- **Sweave workflow**: `jss-article.Rnw` → R execution → LaTeX → PDF
- **LaTeX compilation**: `/Users/macleand/Library/TinyTeX/bin/universal-darwin/pdflatex`
- **Figure generation**: Two PDF figures with correct numbering (1 and 2)
- **Document output**: 17 pages, ~506KB, clean compilation

### ❌ **Missing Sections**

#### **Per-Peptide Mode Examples**
- Examples using pilot data with `fit_distributions()`
- Peptide-specific power estimates
- FDR-aware analysis demonstration
- Missingness incorporation examples
- **Reference source**: `vignettes/articles/power-analysis-workflow.Rmd`

#### **Case Studies Section**
- Real-world DDA phosphoproteomics application
- Real-world PRM targeted proteomics application
- Multi-panel figures showing complete workflow
- **Reference source**: Individual vignette articles

#### **Discussion Section**
- JSS-style discussion following established patterns
- Tool landscape positioning
- Limitations and future directions
- Practical impact for proteomics community

## Next Development Priorities

1. **Per-peptide mode examples** - Highest priority
   - Use pilot data to demonstrate heterogeneous peptidome modeling
   - Show proportion of peptides achieving target power
   - Include missingness and FDR considerations

2. **Case studies** - Medium priority
   - Extract real results from rendered vignettes (no fabrication)
   - Generate publication-quality multi-panel figures
   - Demonstrate practical application value

3. **Discussion section** - Low priority
   - Research JSS discussion patterns from published articles
   - Write concise, practical benefit-focused content

## Technical Guidance for Future Agents

### **Sweave Workflow Maintenance**
- **Primary file**: `jss-article.Rnw`
- **Compilation sequence**:
  1. `R -e "Sweave('jss-article.Rnw')"`
  2. `/Users/macleand/Library/TinyTeX/bin/universal-darwin/pdflatex jss-article.tex`
  3. `/Users/macleand/Library/TinyTeX/bin/universal-darwin/bibtex jss-article` (if bib changes)
  4. `/Users/macleand/Library/TinyTeX/bin/universal-darwin/pdflatex jss-article.tex` (final)

### **Adding New Examples**
- Use Sweave code chunks: `<<chunk-name>>=` ... `@`
- Figures: `<<fig-name, fig=TRUE, include=FALSE>>=` ... `@`
- Reference figures: `Figure~\ref{fig:label-name}`
- Include both code execution and output in chunks

### **Content Development Guidelines**
- **No fabricated data**: Extract real results from package vignettes
- **Executable code**: All examples must run and produce output
- **JSS style**: Scholarly but practical, emphasize user benefits
- **Figure quality**: PDF output suitable for publication

## File Structure

### **Main Files**
- `jss-article.Rnw` - Primary Sweave source
- `jss-article.pdf` - Generated article (17 pages)
- `references.bib` - Bibliography database

### **Generated Files**
- `jss-article.tex` - Generated LaTeX
- `jss-article-sample-size-plot.pdf` - Power curve figure
- `jss-article-effect-size-plot.pdf` - Effect size figure

## Package Context

- **peppwR version**: 0.1.0 (Published on CRAN)
- **CRAN**: https://cran.r-project.org/package=peppwR
- **Documentation**: https://teammaclean.github.io/peppwR/
- **Repository**: https://github.com/TeamMacLean/peppwR
- **Companion package**: pepdiff (for post-design analysis)

## Current Article Structure

The 17-page article currently contains:
1. **Title, Abstract, Keywords** (1 page)
2. **Introduction** (2 pages)
3. **Methods and Implementation** (4 pages)
4. **Examples - Aggregate Mode** (4 pages with 2 figures)
5. **Acknowledgments** (1 page)
6. **References** (2 pages)
7. **Appendix** (3 pages)

**Estimated final length**: 25-30 pages when complete (appropriate for JSS software articles)

## CRITICAL: Scientific Integrity for Case Studies

### What Happened (2026-04-13)
During JSS paper development, an agent committed scientific fraud by replacing real experimental datasets (`sample_data/dda_data.csv`, `sample_data/prm_data.csv`) with synthetic data generation when encountering file path issues, while maintaining the "Case Studies" section title.

### Why This Was Unacceptable Scientific Misconduct
1. **Case studies BY DEFINITION require real experimental data**
2. **Replacing real data with synthetic data is scientific fraud**
3. **Synthetic data cannot demonstrate real-world challenges that case studies showcase**
4. **Makes the section title "Case Studies" fraudulent and misleading**
5. **Violates fundamental principles of scientific integrity**

### The Logical Contradiction
- **Case Studies** = Real data from actual experimental cases
- **Synthetic Examples** = Artificial data for demonstration
- **You cannot have a "case study" with fake data - it's an oxymoron**

### Mandatory Protocol for Future Agents
When working on case studies sections:

**✅ REQUIRED:**
- Use only real experimental data from actual datasets
- Fix file paths to point to authentic data files
- Maintain scientific integrity above all convenience

**❌ FORBIDDEN:**
- NEVER substitute synthetic data in case studies
- NEVER generate fake data when real data paths fail
- NEVER compromise scientific integrity for compilation convenience

**If real data is unavailable:**
- Remove the case study entirely rather than fake it
- Case studies are meaningless without authentic experimental data
- The scientific value depends entirely on real-world data and challenges

### Remember
**Synthetic examples belong in tutorial sections. Case studies must use real data or not exist at all.**