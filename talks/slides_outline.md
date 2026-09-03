# Slide outline: "Two Weird Tricks…and What To Do When They're Not Enough"

**Source prose:** `talks/draft_blogpost.md` (487 lines).
**Format:** one slide per `##` header below. *On screen* describes what the audience sees; *Speaker notes* describes what is said. Line references point back to the prose.
**Status:** Act I drafted as a worked example for format. Acts II and III to follow (likely delegated).

---

## Slide 1 — Title

**On screen.** Talk title — *Two Weird Tricks…and What To Do When They're Not Enough*. Subtitle — *Power analysis from one t-test to a peptidome*. Name, affiliation, date. Small peppwR hex sticker bottom-right if available.

**Speaker notes.** Welcome and topic. Frame the click-baity title honestly: the "tricks" are real but the talk is about statistics. No mention of any blog or prior write-up — the talk stands alone.

---

## Slide 2 — Roadmap

**On screen.** Three labelled boxes, left-to-right, equal weight: **Power analysis** / **peppwR** / **pepdiff**. One-line gloss under each: *"the classical recipe and where it runs out"* / *"power analysis for high-throughput, skewed, missing data"* / *"the analysis end of the same workflow"*. No numbering, no "Part 1 / 2 / 3" hierarchy — three peers.

**Speaker notes.** Set expectations. Three connected topics, treated as peers. Power analysis frames the problem; peppwR is the planning-stage tool; pepdiff is the analysis-stage tool. No live demos. Do NOT use the words "Act I/II/III" or describe any section as a "tease" / "short" / "brief" — every section is a section.

---

## Slide 3 — Section divider: Power analysis

**On screen.** Full-bleed title slide: *Power analysis* — *the chance of finding something given that it is there*. Subtitle in smaller type.

**Speaker notes.** Transition into the first section. Tee up what's coming: a definition, the four levers, the standard $t$-test curve, the two tricks. No reference to any prior write-up.

---

## Slide 4 — What "power" actually means

**On screen.** Large headline: **"Power is the chance of finding something given that it is there."** Below it, smaller, the formal restatement: *In a test of $H_0$ vs $H_1$, power is the probability of correctly rejecting $H_0$ when $H_1$ is true.*

**Speaker notes.** Anchor the informal sentence first. Then translate to formal. Flag $\alpha$ (false positive rate, chosen) vs $\beta$ (false negative rate, *not* chosen by analyst). (Prose lines 18–20.)

---

## Slide 5 — Power = 1 − β

**On screen.** Centred equation: $\text{power} = 1 - \beta$. Below: a single sentence — *"A test with power 0.3 misses 70% of real effects of the size it was designed to detect."*

**Speaker notes.** Land the equation as the formal handle for the informal sentence on the previous slide. Make the 30% / 70% point as the emotional cost. This is where the audience should sit up. (Prose lines 22–26.)

---

## Slide 6 — The four levers, anatomically

**On screen.** Single annotated diagram of two overlapping bell curves (control grey, treatment red), with the four parameters labelled at the geometric feature where each one lives:

- $\delta$ — horizontal arrow between the two peaks (the mean difference)
- $\sigma$ — horizontal double-arrow on one curve's half-width (within-group SD)
- $n$ — small inset showing the *sampling distribution of the mean* narrowing as $n$ grows (two narrower curves overlaid on the wider ones)
- $\alpha$ — shaded tail region under the null curve, labelled

Caption strip below: *"Four levers. These four letters will return in the peppwR section."*

**Speaker notes.** Use the diagram to *show* where each parameter lives in the geometry of the comparison, rather than reading a list. The audience should leave with $\delta, \sigma, n, \alpha$ anchored to picture-features, not to bullets. Flag that these letters come back in the peppwR section. Do NOT say "Act II". (Prose lines 30–37.)

---

## Slide 7 — Picture two overlapping curves

**On screen.** Figure 1 from the prose: two-panel figure with overlapping bell curves, low power left / high power right, overlap shaded. Caption below: *"Pull apart, narrow down, or sample more — all of them shrink the overlap."*

**Speaker notes.** Walk the four bullets verbally against the figure: pull apart ($\delta$↑), narrow ($\sigma$↓), sample more ($n$↑), be stricter ($\alpha$↓ → power ↓). Don't show the bullets on screen; they live in the talk track. (Prose lines 39–57.)

---

## Slide 8 — $p$ is not power

**On screen.** Large headline: *"$p$ is not power, and a small $p$ does not measure how real an effect is."* No other content — punctuation slide.

**Speaker notes.** Flag this as the most common confusion in the audience's bench reading. Set up the next slide as the explanation. (Prose lines 61–63.)

---

## Slide 9 — What $p$ actually is

**On screen.** Figure 2 from the prose: null distribution of the control mean, vertical line at observed treatment mean, shaded tail = $p$. Caption: *"How often would the control alone produce a mean this extreme?"*

**Speaker notes.** Talk the audience through the picture. $p$ is a frequency under an assumed null, not a probability that the alternative is true. Power asks the *inverse* question — given the alternative is true, how often do we catch it. (Prose lines 65–79.)

---

## Slide 10 — Don't compute power after the fact

**On screen.** Headline: *"Power is a property of a design, not of a result."* Sub-bullet: *Post-hoc power on the observed effect adds no new information.* Footnote-style citation: *Hoenig & Heisey (2001)*.

**Speaker notes.** Brief, firm. State the rule, name the canonical reference, move on. (Prose lines 81–83.)

---

## Slide 11 — Cohen's $D$

**On screen.** Equation, centred: $D = \delta / \sigma$. Below, a small table:
| $D$ | label |
|---|---|
| 0.2 | small |
| 0.5 | medium |
| 0.8 | large |

Caption: *"Mean difference, in units of SD. Rules of thumb, not laws."*

**Speaker notes.** Introduce Cohen's $D$ as the standardised composite that makes the next slide's curves comparable across experiments. (Prose lines 87–93.)

---

## Slide 12 — The classical power curve

**On screen.** Figure 3 from the prose: three power-vs-$n$ curves for $D = 0.2, 0.5, 0.8$. Horizontal dashed line at $0.8$. Annotations where each curve crosses $0.8$: $n \approx 393, 64, 26$.

**Speaker notes.** Introduce `pwr.t.test` (the `pwr` package — R's standard closed-form power function) as where these curves come from. Land the punch: even a *large* effect needs ~26 per group to hit 80% power. The bench norm of $n = 3$ is off the left edge. (Prose lines 95–109.)

---

## Slide 13 — Is my old experiment worthless?

**On screen.** Headline: *"Almost certainly not."* Sub-bullet: *A small-$n$ comparison that returned a clear effect is still the data it always was.* Second sub-bullet, in a contrasting colour: *But high-throughput experiments quietly omit large fractions of true biology when underpowered.*

**Speaker notes.** Reassure first, then complicate. The reassurance applies to single comparisons; the warning applies to the kind of experiment the audience increasingly runs. (Prose lines 111–113.)

---

## Slide 14 — Trick one: merge technical replicates honestly

**On screen.** Two-panel cartoon, side by side.
- **Left panel:** one mouse icon → arrow → four small tube icons labelled "tech rep 1–4". Underneath: *"$n_{\text{bio}} = 1$, four measurements — this is $n = 1$, not $n = 4$."*
- **Right panel:** four mouse icons → arrows → four tube icons labelled "bio rep 1–4". Underneath: *"$n_{\text{bio}} = 4$ — this is the $n$ that carries statistical leverage."*

Banner across the bottom of the slide: **"Counting tech reps as biological is $n$-hacking, which is $p$-hacking."**

**Speaker notes.** Use the cartoon to make the distinction concrete before defining it in words. Then walk the mechanism verbally: extra rows shrink SE artificially, which shrinks $p$, which manufactures significance. (Prose lines 115–119.)

---

## Slide 15 — Ordinal data needs a different merge

**On screen.** Headline: *"For ordinal scores, rank first — then average."* Small visual: a column of 0–6 scores → a column of ranks → a single per-biological-rep mean rank.

**Speaker notes.** Brief. The arithmetic mean of an ordinal score is not meaningful; the mean rank is. Sets up the worked example. (Prose lines 122–123.)

---

## Slide 16 — Trick two: decide $n$ before the experiment

**On screen.** Top half: a four-box left-to-right flow diagram with arrows between boxes:

1. *Pilot data or literature* → estimates of $\delta$ and $\sigma$
2. *Choose target* — power (e.g. $0.8$) and $\alpha$ (e.g. $0.05$)
3. *Solve for $n$* — `pwr.t.test(...)`
4. *Run experiment* — biological $n$ as planned; merge technical reps honestly

Bottom strip: headline *"Run the power calculation before you spend the money."*

**Speaker notes.** Walk the boxes left-to-right as the bench-level recipe. The flow makes the order load-bearing — pilot/literature comes before the calculation, the calculation determines the design, the design constrains the bench. (Prose lines 125–127.)

---

## Slide 17 — A real anonymous example

**On screen.** Headline: *"HR scoring, ordinal 0–6, multiple genotype pairs."* Sub-bullets describing the design: technical replicates pooled into $n$; ordinal scores plotted as violins; pairwise $t$-tests across many comparisons. No numbers yet — those are on the next slide.

**Speaker notes.** Frame the example: real TSL paper, kept anonymous, perfectly well-intentioned. The fix is procedural. (Prose lines 129–131.)

---

## Slide 18 — What the re-analysis showed

**On screen.** Figure 4 from the prose (left/right comparison). Headline above the figure: *"$75 / 120$ significant → $42 / 120$ significant."* Sub-bullet: *Rank-and-average within biological replicate, then use biological $n$ only.*

**Speaker notes.** The two numbers are the punchline. A third of the calls were $n$-hacking inflation. Note this is not a small correction. (Prose lines 132–146.)

---

## Slide 19 — Bridge: when the textbook recipe runs out

**On screen.** Headline at top: *"All of that assumed one $t$-test on normal, complete data."*

Below, a row of three small thumbnail figures, equally spaced, each with a short caption:

- **Left thumbnail.** Many faint overlapping curves stacked behind each other. Caption: *"Thousands of simultaneous comparisons."*
- **Middle thumbnail.** A right-skewed histogram with a long tail and a small zero-spike. Caption: *"Strictly positive, multi-order-of-magnitude, right-skewed."*
- **Right thumbnail.** A peptide-by-sample tile matrix with substantial grey gaps for missing values. Caption: *"Non-random missingness."*

Footer strip: **"The closed form does not survive. Simulation does."**

**Speaker notes.** Use the three thumbnails to *show* the three problems rather than asserting them in bullet form. Tee up the next section (peppwR) as the simulation-based response. Do NOT name specific distributions on this slide — that motivation belongs inside the peppwR section, after the audience has looked at real data. Do NOT use the words "Act I/II". (Prose lines 148–154.)

---

# peppwR

## Slide 20 — Section divider: peppwR

**On screen.** Full-bleed title slide: *peppwR* — *power analysis for high-throughput, skewed, missing data*. Small subtitle: *Simulation in place of the closed form.*

**Speaker notes.** Transition out of the closed-form recipe and into the package. The vocabulary from the *Power analysis* section ($\delta, \sigma, n, \alpha$) carries through unchanged; only the way the probability is computed changes. (Prose lines 158–162.)

---

## Slide 21 — Start with the histogram, not the family

**On screen.** Figure 5 from the prose: histogram of phosphopeptide intensities from the DDA *Arabidopsis thaliana* dataset, single condition, $2{,}228$ phosphopeptides across three biological replicates per timepoint. $\log_{10}$ intensity on $x$, count on $y$, long right tail and detection-floor spike annotated.

**Speaker notes.** Introduce DDA — data-dependent acquisition, the standard discovery mass-spec workflow in which the instrument picks precursors for fragmentation in real time. Make the point in front of the picture: this is not a bell curve. Strictly positive, multi-order-of-magnitude, right-skewed, mass near the detection limit. (Prose lines 160–176.)

---

## Slide 22 — The four assumptions of the closed form, named

**On screen.** Centred equation $\pi(n, \delta, \sigma, \alpha)$ at the top. Below, a 2$\times$2 grid of the four assumptions of the $t$-test power formula, each as a short label and a one-line gloss:

- **Normality** — each group is normally distributed
- **Homogeneous variance** — both groups share a common $\sigma$
- **Complete data** — every observation is present
- **Single test** — one comparison, not many

**Speaker notes.** Park the four assumptions on screen so the next slide can knock them down one by one. The four symbols are the same four symbols from earlier — only the validity of the recipe is in question. (Prose lines 187–194.)

---

## Slide 23 — All four assumptions fail at once

**On screen.** Same 2$\times$2 grid as the previous slide, but each cell now overlaid with a red strike and a one-line reason:

- **Normality** — marginal intensities are strongly right-skewed (Figure 5)
- **Homogeneous variance** — variance scales with the mean
- **Complete data** — NA rates above 50% in discovery mode
- **Single test** — thousands of peptides tested simultaneously

**Speaker notes.** Walk the four cells. The first failure is already visible from the previous slide's histogram; the other three will appear in the next few slides. (Prose lines 195–201.)

---

## Slide 24 — Candidate distribution families

**On screen.** Three small density panels, side by side, each labelled:

- **Gamma** — continuous, positive, right-skewed; controlled by shape and rate
- **Lognormal** — variable whose logarithm is normally distributed; natural model for multiplicative processes
- **Inverse Gaussian** — positive, right-skewed; tail behaviour between the other two

**Speaker notes.** Introduce the three parametric families peppwR fits to per-peptide data. Each captures the shape of the histogram from a couple of slides back reasonably well; none is universally right. (Prose lines 178.)

---

## Slide 25 — Picking between families: AIC

**On screen.** Centred equation:

$$\text{AIC}_k = -2 \log L_k + 2 p_k$$

Below, in smaller type: *$L_k$ is the maximised likelihood of model $k$, $p_k$ its number of parameters. Smallest AIC wins.*

**Speaker notes.** Introduce the Akaike Information Criterion as the model-selection criterion that rewards fit and penalises complexity. peppwR fits each candidate family per peptide and picks by AIC. (Prose lines 178–184.)

---

## Slide 26 — Fitting per peptide: `fit_distributions()`

**On screen.** A single code block:

```r
fits <- fit_distributions(
  pilot_data,
  id    = "peptide",
  group = "condition",
  value = "abundance"
)
```

Below the code, one line: *Returns a `peppwr_fits` S3 object with `print`, `plot`, `summary` methods.*

**Speaker notes.** Introduce `fit_distributions()` as the entry point for per-peptide work. One distribution chosen per peptide, AIC-selected from the candidates on the previous slide. Each peptide's missingness rate $r_i$ is tracked alongside its fit. (Prose lines 283.)

---

## Slide 27 — Missingness, and why it is not random

**On screen.** Figure 6 from the prose: per-peptide mean log abundance on $x$, per-peptide NA rate on $y$, each dot a peptide. Negative trend annotated. Spearman $\rho$ printed in the corner.

**Speaker notes.** Walk the picture before naming the mechanism. Low-abundance peptides drop out more often. The cause is mechanical: the instrument has a detection limit. (Prose lines 204–217.)

---

## Slide 28 — MNAR vs MCAR

**On screen.** Two side-by-side labels with one-line glosses:

- **MNAR** — *Missing Not At Random.* The probability of missingness depends on the value itself.
- **MCAR** — *Missing Completely At Random.* The probability of missingness is independent of the value.

Footer strip: *"Pretending MCAR overstates what the experiment can detect."*

**Speaker notes.** Name the two regimes. MNAR is what mass-spec data look like; MCAR is the convenient assumption that power calculations have historically made. (Prose lines 219.)

---

## Slide 29 — peppwR detects MNAR at the dataset level

**On screen.** Headline: *"Dataset-level MNAR check."* Below: a short pseudocode-style summary —

- Compute per-peptide mean abundance and per-peptide NA rate.
- Spearman correlation between the two.
- Significantly negative $\rightarrow$ flag dataset as MNAR.

**Speaker notes.** This is the diagnostic behind Figure 6. The flag drives the optional missingness-aware simulation path on later slides. (Prose lines 221.)

---

## Slide 30 — The simulation recipe, one peptide at a time

**On screen.** Five-step numbered list, no equations yet:

1. Fit a distribution to pilot data for peptide $i$ $\rightarrow$ $\hat\theta_i$.
2. Draw $n$ control values; draw $n$ treatment values from the same distribution scaled by $\delta$.
3. Optionally knock out a fraction $r_i$ of values to mimic MNAR.
4. Run the chosen test; record whether $p \leq \alpha$.
5. Repeat $n_\text{sim}$ times. The fraction of $p \leq \alpha$ outcomes is the Monte Carlo power estimate.

**Speaker notes.** The whiteboard recipe. Walk down the list; the next slide formalises step 5. The same four symbols carry through from earlier. (Prose lines 224–232.)

---

## Slide 31 — Monte Carlo power, formalised

**On screen.** Centred equation:

$$\hat\pi_i(n, \delta) \;=\; \frac{1}{n_\text{sim}} \sum_{j=1}^{n_\text{sim}} \mathbf{1}\!\left[p^{(j)} \leq \alpha\right]$$

Below, smaller: *$\mathbf{1}[\cdot]$ is the indicator function. The same $n, \delta, \alpha$ as before; $\sigma$ enters implicitly through the fitted distribution.*

**Speaker notes.** Land the estimator. This is the line that replaces the closed-form $\pi(n, \delta, \sigma, \alpha)$. The probability is now computed by counting, not by integration. (Prose lines 233–239.)

---

## Slide 32 — Three tests inside the engine

**On screen.** Three labelled rows, each with a one-line gloss:

- **Wilcoxon rank-sum** — rank-based two-sample comparison; robust to non-normality; the peppwR default.
- **Bootstrap-$t$** — resampling variant of the two-sample $t$-test; no parametric null required.
- **Bayes-factor $t$-test** — ratio of marginal likelihoods under $H_1$ and $H_0$; rejection threshold expressible in $\alpha$-equivalent terms.

**Speaker notes.** The choice of test inside step 4 of the recipe matters at small $n$. Flag that pepdiff exposes the same family at the analysis stage — picking the same test at both stages keeps the power calculation honest. (Prose lines 241.)

---

## Slide 33 — The three questions, one engine

**On screen.** Three-column layout, equal weight:

- **Sample size** — given $\delta, \sigma, \pi_\text{target}$, find smallest $n$. `find = "sample_size"`
- **Power** — given $n, \delta, \sigma$, return $\hat\pi$. `find = "power"`
- **Effect size** — given $n, \sigma, \pi_\text{target}$, return smallest detectable $\delta$. `find = "effect_size"`

**Speaker notes.** Three questions a researcher actually asks. Each is an optimisation over the same simulation engine, exposed via the `find` argument of `power_analysis()`. The return value is a `peppwr_power` S3 object. (Prose lines 244–251.)

---

## Slide 34 — Aggregate mode: no pilot data

**On screen.** Single code block:

```r
result <- power_analysis(
  distribution = "gamma",
  params       = list(shape = 2, rate = 0.05),
  effect_size  = 2,
  target_power = 0.8,
  find         = "sample_size",
  n_sim        = 5000
)
```

Below: *Returns $n^\ast = 25$ per group.*

**Speaker notes.** Aggregate mode is the grant-stage tool: distribution family and parameters supplied directly. Order-of-magnitude answer; enough to bound the budget conversation. (Prose lines 253–268.)

---

## Slide 35 — Aggregate-mode power curve

**On screen.** Figure 7 from the prose: $n$ per group on $x$, estimated power on $y$. Horizontal dashed line at $0.8$. Vertical dashed line at $n^\ast = 25$ where the curve first crosses the target.

**Speaker notes.** Read the picture: the curve rises, hits target, and the function reports the crossing $n$. The same shape as the classical $t$-test curve from earlier, but produced by simulation rather than a closed form. (Prose lines 269–279.)

---

## Slide 36 — Per-peptide mode: with pilot data

**On screen.** Headline: *"One fit per peptide, one simulation per peptide."* Below, three short bullets:

- Fit gamma / lognormal / normal / inverse Gaussian per peptide; AIC-select.
- Track each peptide's missingness rate $r_i$ alongside the fit.
- Run the Monte Carlo recipe independently per peptide.

**Speaker notes.** The mode that matches how a real phosphoproteomics experiment behaves: heterogeneously, peptide by peptide. No single global power number. (Prose lines 281–283.)

---

## Slide 37 — The peptidome fraction

**On screen.** Centred equation:

$$P_\text{peptidome}(n, \delta, \pi_\text{target}) \;=\; \frac{1}{K}\sum_{i=1}^{K} \mathbf{1}\!\left[\hat\pi_i(n, \delta) \geq \pi_\text{target}\right]$$

Below, smaller: *$K$ is the number of peptides. Headline format: "X% of peptides reach 80% power at $n = Y$".*

**Speaker notes.** The summary that comes out of per-peptide mode. Not a single power number; a fraction of the peptidome reaching target. (Prose lines 285–289.)

---

## Slide 38 — When per-peptide fitting fails

**On screen.** Three boxes, left to right, each with a label and a one-line gloss:

- **`"exclude"`** (default) — drop the peptide, report the count.
- **`"empirical"`** — resample empirically from the observed values.
- **`"lognormal"`** — substitute a moment-matched lognormal.

Footer: *Configured via `on_fit_failure`.*

**Speaker notes.** Small samples, extreme values and numerical trouble all cause fits to fail for individual peptides. The three fallbacks let the analyst choose how to handle them rather than silently dropping data. (Prose lines 291.)

---

## Slide 39 — Multiple testing, again

**On screen.** Large headline: *"At $\alpha = 0.05$ across $5{,}000$ tests, expected false positives under the null = $250$."* Below, smaller: *Not a usable workflow.*

**Speaker notes.** Land the arithmetic before introducing the fix. The audience has seen $\alpha$ before; now it has to be made to compound. (Prose lines 293–295.)

---

## Slide 40 — Benjamini–Hochberg

**On screen.** Centred inequality:

$$p_{(k^*)} \;\leq\; \frac{k^*\, \alpha}{m}$$

Below, smaller: *Ordered $p$-values $p_{(1)} \leq \cdots \leq p_{(m)}$. Find the largest rank $k^*$ satisfying the inequality; reject every hypothesis with $p \leq p_{(k^*)}$.*

**Speaker notes.** Introduce the Benjamini–Hochberg procedure as the standard fix in genomics and proteomics. It controls the false discovery rate — the expected proportion of false positives among the calls declared significant. Threshold is adaptive: peptides compete for the significance budget. (Prose lines 297–303.)

---

## Slide 41 — FDR-aware power: `apply_fdr = TRUE`

**On screen.** Headline: *"Power under BH, not in isolation."* Below, a four-step loop diagram:

1. Simulate a full peptidome (true positives + nulls).
2. Apply BH.
3. Count true positives among the discoveries.
4. Repeat.

Footer strip: *`prop_null = 0.9` by default — most peptides do not respond.*

**Speaker notes.** A peptide's chance of being declared significant depends on the rest of the peptidome's $p$-values. peppwR re-estimates per-peptide power *under* the BH procedure when `apply_fdr = TRUE`. (Prose lines 304–305.)

---

## Slide 42 — Worked example: DDA *Arabidopsis*

**On screen.** Headline: *"$2{,}228$ phosphopeptides, two timepoints, three biological replicates per timepoint."* Sub-bullets: $n = 3$ per group, two-fold effect, Bayes-factor test.

**Speaker notes.** Frame the case study. The sample size is entirely typical of discovery phosphoproteomics. The question is what fraction of the peptidome that design can actually power. (Prose lines 307–309.)

---

## Slide 43 — Per-peptide power distribution

**On screen.** Figure 8 from the prose: histogram of estimated per-peptide power across $2{,}228$ phosphopeptides at $n = 3$. $x$ = power per peptide, $y$ = count. Vertical line at $0.8$. Most peptides fall below the line.

**Speaker notes.** Read the picture: median power is low, only a small fraction of peptides clears 80%. The discovery experiment is best read as hypothesis-generating, not as a definitive list. (Prose lines 311–321.)

---

## Slide 44 — What $n$ would reach half-peptidome 80% power?

**On screen.** Headline: *"$n = 15$ per group."* Below, smaller: *Roughly five times the discovery norm of three.*

**Speaker notes.** The inverse question — and the kind of quantitative argument a power analysis is supposed to produce. Not a verdict on the existing experiment; a budget for the next one. (Prose lines 323.)

---

## Slide 45 — PRM case study: missingness on, missingness off

**On screen.** Headline: *"PRM $Magnaporthe$ $oryzae$, $285$ targeted phosphopeptides."* Sub-bullets:

- Three biological replicates per timepoint, two technical replicates averaged.
- Dataset-level MNAR diagnostic: clear negative $\rho$.
- Run with `include_missingness = TRUE` vs `FALSE`.

**Speaker notes.** Introduce PRM — parallel reaction monitoring, the targeted mass-spec workflow flagged earlier as the contrast to DDA. Introduce *Magnaporthe oryzae* as the rice blast fungus. The diagnostic is the same one from the MNAR slide; turning missingness on shows how much of the nominal power was an artefact of pretending no values were missing. (Prose lines 325.)

---

# pepdiff

## Slide 46 — Section divider: pepdiff

**On screen.** Full-bleed title slide: *pepdiff* — *the same data, from the other end*. Subtitle: *Differential abundance for phosphopeptides.*

**Speaker notes.** Transition. The planning question is settled; the spreadsheet of intensities has landed; the same data quirks — skew, MNAR, multiple testing — are now the analyst's problem. (Prose lines 333–339.)

---

## Slide 47 — Input data shape

**On screen.** A small mock data frame with columns: *peptide, gene, treatment, timepoint, bio_rep, abundance*, and three or four illustrative rows.

**Speaker notes.** One row per peptide per sample. Columns naming the peptide, the gene, the experimental factors, the biological replicate index, and the measured abundance. (Prose lines 337.)

---

## Slide 48 — The workflow at a glance

**On screen.** Centred chain:

$$\text{CSV} \;\longrightarrow\; \texttt{read\_pepdiff()} \;\longrightarrow\; \texttt{pepdiff\_data} \;\longrightarrow\; \texttt{compare()} \;\longrightarrow\; \texttt{pepdiff\_results} \;\longrightarrow\; \texttt{plot()}$$

**Speaker notes.** Walk the chain left to right. Two S3 classes, two named methods between them, one plotting endpoint. (Prose lines 341–349.)

---

## Slide 49 — Reading the data

**On screen.** Single code block:

```r
dat <- read_pepdiff(
  file      = "phospho_intensities.csv",
  id        = "peptide",
  gene      = "gene_id",
  value     = "value",
  factors   = c("treatment", "timepoint"),
  replicate = "bio_rep"
)
```

Below: *Returns a `pepdiff_data` S3 object with `print`, `summary`, `plot` methods.*

**Speaker notes.** The column-name arguments are how `read_pepdiff()` knows which column plays which role. The result is a typed object the rest of the pipeline can dispatch on. (Prose lines 351–362.)

---

## Slide 50 — Exploratory triptych

**On screen.** Figure 10 from the prose: three side-by-side panels.

- **Left.** PCA of samples, coloured by treatment, shaped by timepoint.
- **Centre.** Per-sample density curves of log intensity, overlaid.
- **Right.** Per-peptide mean abundance vs NA rate.

Caption strip: *`plot()` on a `pepdiff_data` object — diagnostics before any model is fitted.*

**Speaker notes.** The exploratory plots run before any modelling. The right-hand panel is the same MNAR diagnostic introduced at the planning stage. Large horizontal shifts in the centre panel flag a normalisation problem. (Prose lines 364–378.)

---

## Slide 51 — The per-peptide gamma GLM

**On screen.** Two lines of equation, centred:

$$Y_{ij} \;\sim\; \text{Gamma}(\mu_{ij}, \phi_i)$$

$$\log \mu_{ij} \;=\; \mathbf{x}_j^{\top} \boldsymbol{\beta}_i$$

Below, smaller: *Peptide $i$, sample $j$. Design vector $\mathbf{x}_j$; coefficient vector $\boldsymbol{\beta}_i$; dispersion $\phi_i$.*

**Speaker notes.** Introduce the generalised linear model — regression in which the response follows an exponential-family distribution, linked to a linear predictor through a chosen link function. Family is gamma — the same distribution met at the planning stage. Link is the logarithm. Note the bold subscripted $\boldsymbol{\beta}_i$ is the GLM coefficient, distinct from the Type II error rate $\beta$ already in play. (Prose lines 380–387.)

---

## Slide 52 — Why gamma, why log link

**On screen.** Three short bullets:

- **No pooled variance.** Each peptide has its own $\phi_i$ — variance is allowed to scale with the mean.
- **Strict positivity at the likelihood.** No pseudo-count, no log-transform of the response.
- **$\exp(\beta_{i,\text{trt}})$ is the fold change.** Coefficients interpret as log fold changes directly.

**Speaker notes.** Three properties pulled out of the previous slide's two lines. The log link is what makes the coefficients interpretable as fold changes without further arithmetic. (Prose lines 388.)

---

## Slide 53 — Contrasts via `emmeans`

**On screen.** Three short bullets:

- Treatment minus control at a fixed timepoint.
- Treatment $\times$ timepoint interaction.
- Marginal effect of treatment averaged over timepoints.

Footer: *Extracted from the fitted model via the `emmeans` package — R's standard tool for estimated marginal means and contrasts.*

**Speaker notes.** Each contrast yields a per-peptide effect size and a per-peptide $p$-value. (Prose lines 392.)

---

## Slide 54 — Running the comparison

**On screen.** Single code block:

```r
results <- compare(
  dat,
  compare = "treatment",
  ref     = "ctrl",
  method  = "glm"
)
```

Below: *Returns a `pepdiff_results` object.*

**Speaker notes.** One call, one returned S3 object. The `results` table holds, per peptide and per contrast: fold change, $\log_2$ fold change, raw $p$-value, BH-adjusted FDR, and a Boolean `significant` flag at FDR $\leq 0.05$. (Prose lines 394–403.)

---

## Slide 55 — Multiple testing, again

**On screen.** Same BH inequality as before:

$$p_{(k^*)} \;\leq\; \frac{k^* \alpha}{m}$$

Below: *Stored in the `fdr` column. Applied within each comparison group via `apply_fdr_by_comparison()` when designs are nested.*

**Speaker notes.** The same fix used at the planning stage. The within-comparison default keeps a peptide from having its significance budget eaten by tests from a different stratum; pooling all comparisons before correcting is more conservative and is a defensible alternative. (Prose lines 405–415.)

---

## Slide 56 — When the GLM does not fit: diagnostics

**On screen.** Figure 11 from the prose: four-panel `plot_fit_diagnostics()` output.

- **Top-left.** Histogram of per-peptide deviances, 95th-percentile threshold marked.
- **Top-right.** Deviance vs $\log_2$ fold change, poorly fitted peptides highlighted.
- **Bottom-left.** Standardised residuals vs fitted values for three example peptides (high, median, low deviance).
- **Bottom-right.** Pooled QQ plot of standardised residuals across the peptidome.

**Speaker notes.** Walk the four panels. Heavy tails or S-shapes in the QQ plot, or systematic curvature in the residual plots, are the signal to switch methods. (Prose lines 419–435.)

---

## Slide 57 — The decision rule

**On screen.** Headline: *"When to switch."* Two contrasting bullets:

- **Below ~10–15% high-deviance peptides + roughly linear QQ** $\rightarrow$ keep the gamma GLM.
- **Above ~15% high-deviance, or systematic curvature** $\rightarrow$ switch to ART.

**Speaker notes.** The documented rule from the package. Conservative by design. (Prose lines 437.)

---

## Slide 58 — Aligned Rank Transform (ART)

**On screen.** Three-step recipe, numbered:

1. *Align* — subtract all other estimated effects so residuals reflect only the effect of interest.
2. *Rank* — replace aligned values by ranks.
3. *ANOVA on ranks* — ordinary factorial ANOVA on the ranked, aligned values.

Footer: *Handles interactions, which classical rank-based tests cannot. Returns rank-based statistics, not fold changes.*

**Speaker notes.** Nonparametric procedure for factorial designs. Alignment is what makes interactions tractable; ranking is what frees the user from picking a parametric family. The price is interpretive: for effect-size reporting, keep the GLM fit alongside the ART significance call. (Prose lines 437–445.)

---

## Slide 59 — Running ART

**On screen.** Single code block:

```r
results_art <- compare(
  dat,
  compare = "treatment",
  ref     = "ctrl",
  method  = "art"
)
```

Below: *Same `pepdiff_results` shape, same `fdr` column, same downstream plot methods.*

**Speaker notes.** Identical call shape, identical output shape. The switch is one argument. (Prose lines 441–445.)

---

## Slide 60 — Pairwise tests for simpler designs

**On screen.** Headline: *`method = "pairwise"`.* Below, four labelled rows:

- **Wilcoxon** — documented default for small-$n$, non-normal, outlier-prone data.
- **Bootstrap-$t$** — four or more replicates, roughly symmetric data.
- **Bayes-factor $t$-test** — evidence strength; useful when "no detectable effect" is the finding.
- **Rank-products** — sensitive to small effects whose direction is consistent across replicates.

**Speaker notes.** For two-group, no-factorial-structure, small-$n$ designs, neither a GLM nor ART is strictly necessary. The first three tests are the same family that peppwR's per-peptide simulation engine supports via its own `test` argument. (Prose lines 447–451.)

---

## Slide 61 — Match the test at both stages

**On screen.** Two-row table, side by side:

- **Planning stage (peppwR).** `test = "wilcoxon"` / `"bootstrap_t"` / `"bayes_t"`
- **Analysis stage (pepdiff).** `test = "wilcoxon"` / `"bootstrap_t"` / `"bayes_t"` / `"rank_products"`

Footer: *"Power for the test you will actually run."*

**Speaker notes.** Picking the same test at both stages means the power calculation estimates power for the test that will actually be run, not for a convenient closed-form proxy. (Prose lines 451.)

---

## Slide 62 — Reading the output

**On screen.** Figure 12 from the prose: three-panel `plot(<pepdiff_results>)` output.

- **Left.** Volcano plot of $-\log_{10}(\text{FDR})$ vs $\log_2$ fold change. FDR $0.05$ line and fold-change threshold marked. Significant peptides coloured.
- **Centre.** Histogram of raw $p$-values.
- **Right.** Density of $\log_2$ fold changes across all peptides.

**Speaker notes.** Three views, returned together. The volcano plot is the customary endpoint; the other two are diagnostics. (Prose lines 455–471.)

---

## Slide 63 — The $p$-value histogram does the diagnostic work

**On screen.** Three small thumbnail panels, each captioned:

- **Flat base + spike near zero** — well-calibrated test, mix of nulls and alternatives.
- **U-shape** — something upstream of FDR is mis-specified.
- **Spike at one or non-flat baseline** — wrong family, unmodelled batch effect, normalisation problem.

**Speaker notes.** The quiet diagnostic. The shape of the histogram is what tells the analyst whether to trust the FDR column. (Prose lines 473.)

---

## Slide 64 — Filtered output

**On screen.** Single line of code: `significant(results)`. Below: *Returns the filtered peptide list, ready to join onto gene annotations or pipe into enrichment.*

**Speaker notes.** The endpoint of the pipeline. Filtered list out, downstream analyses next. (Prose lines 475.)

---

# Closing

## Slide 65 — Takeaways

**On screen.** Six bullets, no further hierarchy:

- Power is a property of a design, not a result. Compute it before, not after.
- Independent biological replicates carry the leverage. Pooling technical reps as biological ones is $n$-hacking.
- The closed-form $t$-test power formula does not survive contact with phosphoproteomics data. Simulation does.
- peppwR: aggregate mode for grant-stage estimates, per-peptide mode for peptidome-level answers; MNAR-aware and BH-FDR-aware variants on request.
- pepdiff: per-peptide gamma GLM with log link for interpretable fold changes; ART as the robust fallback; pairwise tests for simple designs; BH-FDR throughout.
- Match the test at the analysis stage to the test assumed at the planning stage.

**Speaker notes.** Read down the bullets at a steady pace. The first three are the framing; the next two are the tools; the last is the operational rule. (Prose lines 481–487.)

---

## Slide 66 — One sentence each

**On screen.** Two lines, large type:

- **peppwR** — *for the experiment you are about to run.*
- **pepdiff** — *for the one you have already finished.*

**Speaker notes.** The single-sentence handle for each package. Last word before questions. (Prose line 487.)

---

## Slide 67 — Questions

**On screen.** Centred: *Questions?* Below, smaller, two URLs:

- `https://teammaclean.github.io/peppwR/`
- `https://teammaclean.github.io/pepdiff/`

**Speaker notes.** Open the floor. Both documentation sites are public; the JSS paper is the citable reference for peppwR.
