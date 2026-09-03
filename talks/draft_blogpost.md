---
title: "Two Weird Tricks…and What To Do When They're Not Enough"
description: "Power analysis from one t-test to a peptidome: a narrative dry run for the talk."
author: Dan MacLean
date: TBD
---

A few years ago I gave a short internal talk to PhD students and post-docs with the deliberately click-baity title *"Two Weird Tricks For Improving Your Experimental Designs"*. The plan was a small bait-and-switch: tempt bench biologists into a room with the promise of a simple hack, and once inside, talk to them about statistics. It mostly worked. Someone even asked for the notes to be put online, which is how the original blog post came to exist.

This piece extends the original argument. It first retells the two tricks in compressed form, because they remain the single most useful thing I know how to tell a bench biologist about experimental design. It then picks up where the original stops: what to do when the experiment is not one comparison on roughly bell-shaped data, but several thousand comparisons on data that look nothing like a bell curve and are missing in a non-random way. And it covers what happens next, after the experiment has been designed and the measurements collected — when the data are in hand and the question becomes which peptides moved.

---

## Act I: Two weird tricks

### What "power" actually means

If only one sentence from this post stays with you, let it be this: **power is the chance of finding something given that it is there**.

The formal definition takes a little more care. In a statistical test of a null hypothesis $H_0$ against an alternative $H_1$, the test can fail in two ways. A *false positive* (Type I error) rejects $H_0$ when it is true; its rate is conventionally written $\alpha$ and fixed by the analyst, usually at $0.05$. A *false negative* (Type II error) fails to reject $H_0$ when $H_1$ is true; its rate is $\beta$, and unlike $\alpha$ it is not something the analyst chooses directly. Power is the complementary quantity,

$$
\text{power} = 1 - \beta,
$$

the probability that the test correctly rejects $H_0$ when $H_1$ is in fact true. A test with power $0.3$ will miss $70\%$ of real effects of the size it was designed to detect. That is a poor return on bench effort.

### The four parameters that determine power

For a simple two-sample comparison of means, power is a function of exactly four quantities:

1. The **effect size**, $\delta$ — the true difference in means between the two groups.
2. The **variability**, $\sigma$ — the standard deviation of the measurement within a group.
3. The **sample size**, $n$ — the number of independent replicates per group.
4. The **significance threshold**, $\alpha$ — the false-positive rate the analyst is willing to tolerate.

These four symbols ($\delta$, $\sigma$, $n$, $\alpha$) are the working vocabulary of power analysis, and they will reappear in Act II in exactly the same roles. It is worth pausing on each.

You can build intuition without any equations. Picture two overlapping bell curves on the same axis, one for the control and one for the treatment.

- Pull the curves apart (larger $\delta$ at fixed $\sigma$): the overlap shrinks, the test is more likely to distinguish them, power rises.
- Squeeze each curve narrower (smaller $\sigma$ at fixed $\delta$): the overlap shrinks for the same reason, power rises.
- Take more samples (larger $n$): the *sampling distributions of the two means* narrow even if the underlying curves do not, the means separate more cleanly, power rises.
- Lower $\alpha$ (be stricter about what counts as a hit): the rejection region shrinks, power falls.

```
+----------------------------------+
|                                  |
|       [FIGURE PLACEHOLDER]       |
|                                  |
+----------------------------------+
Figure 1. Two-panel intuition for power. Left panel: two overlapping
bell curves with small mean separation, large overlap shaded; "low
power". Right panel: same two curves pulled apart, small overlap;
"high power". Control in grey, treatment in red, shaded overlap in
pale red.
```

For continuous correlation-type analyses, the coefficient of determination $R^2$ rolls effect size, variability, and the leverage of the sample together into a single scalar — convenient, though it can hide which lever is doing the work.

### $p$ is not power

A confusion worth heading off early: $p$ and power are not the same thing, and a small $p$ is not a measurement of how big or how reliable an effect is.

The $p$-value is, loosely, a frequency under an assumed null. Imagine you have a control group and a treatment group, each measured for some continuous quantity (this is a generic two-sample comparison; pick your own example). Fit a normal distribution to the control values. Now ask: if I drew $n$ values at random from a population that looked like the control, how often would the *mean* of the draw end up at least as extreme as the observed treatment mean? That long-run frequency is what the $p$-value approximates. It does not tell you the probability that the alternative is true, and it does not tell you whether the experiment was sensitive enough to detect what it was looking for.

```
+----------------------------------+
|                                  |
|       [FIGURE PLACEHOLDER]       |
|                                  |
+----------------------------------+
Figure 2. The intuition behind a p-value. Single normal curve labelled
"null distribution of the control mean". Vertical line at the observed
treatment mean. Shaded tail beyond the line labelled "this much area =
p". Annotation: "How often would the control alone produce this?".
```

Power, by contrast, asks the inverse question. Suppose the alternative *is* true and the real difference is $\delta$. What fraction of experiments of size $n$ would correctly reject the null at threshold $\alpha$? That fraction is $1-\beta$.

A small $p$-value does not retroactively rescue a low-powered experiment, but neither does it invalidate one. If a properly designed study returns $p \leq 0.05$, the inference is what it is; computing power *after the fact* on the observed effect size and using it to argue the result is somehow more or less believable is a well-documented mistake.[^1] Power is a property of an experimental design, not of a result.

[^1]: The post-hoc-power literature is large and slightly weary. Hoenig and Heisey (2001) is the usual reference; the short version is that observed-power calculations are a deterministic transform of the $p$-value and add no new information.

### Power versus $n$, the classical curve

For the two-sample $t$-test on equal-variance normal data, the relationship between $n$, $\delta$, $\sigma$ and $\alpha$ has a closed form. In R, this is what the function `pwr.t.test`, in the `pwr` package — R's standard closed-form power function — computes. To keep effect size and variability on a common footing, the convention is to combine them into *Cohen's $D$*, a standardised composite of effect size and standard deviation,

$$
D = \frac{\delta}{\sigma},
$$

read as "the mean difference, in units of standard deviation". By convention, $D = 0.2$ is called a small effect, $D = 0.5$ a medium effect, and $D = 0.8$ a large one — rules of thumb, not laws.

Plotting power as a function of $n$ for these three values of $D$ produces the curves that get shown at the start of every statistics workshop.

```
+----------------------------------+
|                                  |
|       [FIGURE PLACEHOLDER]       |
|                                  |
+----------------------------------+
Figure 3. Power versus sample size per group, from pwr.t.test, for
Cohen's D = 0.2, 0.5, 0.8. x-axis n per group from 2 to 100, y-axis
power from 0 to 1. Horizontal dashed line at 0.8. Annotate where each
curve crosses 0.8: roughly n ~ 393, 64, 26.
```

The headline that bench biologists almost never expect is the one in the figure caption. Reaching $80\%$ power at $D = 0.8$ — a *large* effect — takes roughly twenty-six biological replicates per group. The corresponding numbers for medium and small effects are larger still. The bench norm of three biological replicates per group sits far to the left of the rising part of every one of these curves. Most experiments, on this view, are underpowered for the effects they are trying to detect.

### Is my old experiment now worthless?

Almost certainly not. A single small-$n$ comparison that returned a clear effect is not retroactively undone by a power calculation — it is still the data it always was. The places this hurts most are *high-throughput* experiments, where many tests run in parallel and the false negative rate compounds across the peptidome (or transcriptome, or whatever the "ome" is for the day). Underpowered high-throughput experiments quietly omit large fractions of true biology, and there is no way to recover what was never declared significant in the first place.

### Trick one: merge technical replicates properly

The first lever every analysis has is $n$. The first mistake every analysis makes is inflating it.

A *technical* replicate is a re-measurement of the same biological material; a *biological* replicate is an independent draw from the population the experiment claims to be about. Pooling four technical replicates from each of two mice and reporting $n = 8$ is not an experiment with eight independent observations of a treatment — it is two independent observations, measured four times each. Independence is what gives an experiment its statistical leverage; counting non-independent measurements as independent shrinks the standard deviation artificially, shrinks the standard error of the mean, and shrinks $p$. It is, in effect, $p$-hacking via the denominator.[^2]

[^2]: I have rarely if ever seen a paper that does this on purpose. Almost every case I have come across has been entirely well-intentioned. The fix is procedural, not moral.

For ordinal categorical scores — disease ratings, severity classes — the merge is slightly more subtle. The right move is to *rank* the technical replicates first and then average the ranks within each independent biological replicate. Treating an ordinal scale as if it were a continuous numeric scale and taking arithmetic means is a separate, additional mistake that often hides inside the same figures.

### Trick two: increase real biological $n$, guided by a power calculation

The second lever is also $n$, but the kind that costs real money. The point of a power calculation is to decide how much of it to spend *before* the experiment runs. With either pilot data or sensible guesses for $\delta$ and $\sigma$ from the literature, `pwr.t.test` will tell you the per-group $n$ needed to achieve a chosen power at a chosen $\alpha$. Run that calculation, plan to that number of biological replicates, and merge technical replicates correctly underneath them. Those are the two tricks.

### A worked anonymous example

To make this concrete, here is a real case from a paper I will keep anonymous. The readout was an ordinal hypersensitive-response (HR) score on a 0–6 scale, with multiple technical replicates per biological replicate, several biological replicates per genotype pair, and many genotype pairs compared in a single figure.

In the published version, the technical replicates were pooled into a per-group $n$ and the ordinal scores were plotted on a violin plot as if they were continuous. Pairwise $t$-tests across the comparisons returned $75$ significant calls out of $120$. After re-analysis — averaging *ranks* within each biological replicate to respect the ordinal scale, and using only the true biological replicates ($n = 2$) for inference — the count fell to $42$. The difference between those two numbers is the false-positive inflation produced by $n$-hacking. It is not a small correction.

```
+----------------------------------+
|                                  |
|       [FIGURE PLACEHOLDER]       |
|                                  |
+----------------------------------+
Figure 4. Anonymised HR-score re-analysis. Left: original violin plot
with ordinal scores treated as continuous, n inflated by pooling
technical replicates; small p-values everywhere. Right: count plot
with point size proportional to frequency, after rank-and-average
merge and biological n only; fewer but more defensible calls.
```

### Bridge: when the textbook recipe runs out

The first half of this post assumed something that almost never holds in practice in molecular biology beyond a single readout. It assumed one $t$-test, on roughly normal data, with all observations present. Modern high-throughput biology violates each of those assumptions at once.

A phosphoproteomics run measures the abundances of thousands of phosphopeptides per sample. The values are strictly positive, span several orders of magnitude, and are missing in patterns that are anything but random. The "experiment" is not a single comparison but several thousand simultaneous ones, each with its own scale, its own variance, and its own missingness profile. The closed-form recipe behind `pwr.t.test` no longer applies, and post-hoc inspection cannot stand in for design.

What does work is simulation. The next section is about how to do that without losing the vocabulary we just built.

---

## Act II: peppwR

### Looking at the data before naming anything

The right place to start is not with a probability distribution but with a histogram. Figure 5 shows the distribution of peptide intensities for one condition of a real data-dependent acquisition (DDA) phosphoproteomics dataset of *Arabidopsis thaliana* — DDA being the standard discovery-mode mass spectrometry workflow in which the instrument selects precursors for fragmentation in real time, as opposed to the targeted parallel reaction monitoring (PRM) workflow we will meet later. The histogram contains the abundance values for $2{,}228$ phosphopeptides across three biological replicates per timepoint.

```
+----------------------------------+
|                                  |
|       [FIGURE PLACEHOLDER]       |
|                                  |
+----------------------------------+
Figure 5. Histogram of phosphopeptide intensities from the DDA
Arabidopsis dataset, single condition. x-axis: log10 intensity (with a
linear inset). y-axis: count. Annotate the long right tail and the
spike near the detection limit.
```

This does not look like a bell curve. The distribution is strictly positive, it is heavily right-skewed, and it has visible mass close to the detection floor. Whatever the standard $t$-test power formula assumes about the data, it is not this.

A handful of parametric families capture the shape we are looking at reasonably well: the *gamma* distribution (continuous, positive, right-skewed, controlled by a shape and a rate); the *lognormal* (the distribution of a variable whose logarithm is normally distributed, a natural model for multiplicative biological processes); and the *inverse Gaussian* (positive, right-skewed, with tail behaviour intermediate between the other two). When more than one of these gives a plausible fit to a given peptide, a tie-breaker is needed. The conventional choice is the *Akaike Information Criterion* (AIC) — a model-selection criterion that rewards goodness of fit and penalises model complexity. Concretely,

$$
\text{AIC}_k = -2 \log L_k + 2 p_k,
$$

where $L_k$ is the maximised likelihood of model $k$ and $p_k$ is its number of parameters. The model with the smallest AIC wins.

### Why the closed form does not survive

The closed-form $t$-test power formula, which we can write generically as

$$
\pi(n, \delta, \sigma, \alpha)
$$

— the same four parameters as Act I, in the same roles — relies on four assumptions: that each group is normally distributed, that the two groups share a common $\sigma$, that every observation is present, and that one test is being conducted. Phosphoproteomics violates all four at once.

- *Normality* fails because the marginal distribution of intensities is strongly right-skewed, as Figure 5 shows.
- *Homogeneous variance* fails because variance in mass spectrometry data tends to scale with the mean, so peptides at different abundance levels do not share a common $\sigma$.
- *Complete data* fails because the experiment routinely returns NA at rates above 50% in discovery mode.
- *Single test* fails because thousands of peptides are tested simultaneously.

Asserting these failures is one thing; we have just seen the first one in a real histogram. The other three follow from the same dataset and will appear in turn.

### Missingness, and why it is not random

Before naming the mechanism, look at the picture. Figure 6 plots each peptide as a point: its mean (log) abundance on the $x$-axis, and its NA rate (the proportion of samples in which it was not detected) on the $y$-axis. The relationship is not flat.

```
+----------------------------------+
|                                  |
|       [FIGURE PLACEHOLDER]       |
|                                  |
+----------------------------------+
Figure 6. Mean abundance versus NA rate, per peptide, PRM dataset.
x-axis: mean log abundance. y-axis: proportion of samples with NA.
Each dot one peptide. Negative trend annotated; Spearman rho printed
in the corner.
```

Low-abundance peptides are systematically more likely to be missing. The cause is mechanical: the instrument has a detection limit, and peptides whose true abundance sits near that limit drop out probabilistically. This pattern is called *missing not at random* (MNAR) — the missingness depends on the very value that is missing — and it stands in contrast to *missing completely at random* (MCAR), in which the chance of a value being missing is independent of its size. Power calculations that pretend the data are MCAR (or that pretend there are no missing data at all) systematically overstate what the experiment can detect, because MNAR removes preferentially the values that carry the least signal-to-noise margin.

The peppwR package detects this pattern at the dataset level by computing the Spearman correlation between per-peptide mean abundance and per-peptide NA rate, and flagging a significant negative correlation as evidence of MNAR.

### Simulation as the way out

If the closed form cannot be trusted, the working alternative is to estimate power by simulating the experiment many times. The whiteboard recipe, for a single peptide $i$, is:

1. Fit a distribution to the pilot data for peptide $i$, giving estimated parameters $\hat\theta_i$.
2. Draw $n$ simulated control values from this distribution and $n$ simulated treatment values from the same distribution scaled by the effect size $\delta$.
3. Optionally, knock out a fraction $r_i$ of the simulated values according to the peptide's observed missingness rate, to mimic MNAR.
4. Run the statistical test of choice on the simulated samples. Record whether $p \leq \alpha$.
5. Repeat steps 2–4 a large number of times. The fraction of repetitions in which $p \leq \alpha$ is the Monte Carlo power estimate.

Formally, over $n_\text{sim}$ replicate simulations indexed by $j$,

$$
\hat\pi_i(n, \delta) \;=\; \frac{1}{n_\text{sim}} \sum_{j=1}^{n_\text{sim}} \mathbf{1}\!\left[p^{(j)} \leq \alpha\right],
$$

where $\mathbf{1}[\cdot]$ is the indicator function. The same four symbols ($n$, $\delta$, $\alpha$, and implicitly $\sigma$ through the fitted distribution) carry through from Act I; only the way the probability is computed has changed.

The choice of test inside step 4 is open. The peppwR package supports three: the Wilcoxon rank-sum test (the rank-based two-sample comparison, robust to non-normality and the default in peppwR); the *bootstrap-$t$* test (a resampling variant of the two-sample $t$-test that does not require a parametric null); and a *Bayes-factor $t$-test* (a Bayesian alternative that returns the ratio of marginal likelihoods under $H_1$ and $H_0$, with a rejection threshold expressible in $\alpha$-equivalent terms). For the small sample sizes characteristic of phosphoproteomics, this choice matters; we will come back to it in the case study.

### The three questions, one engine

A researcher who reaches for a power calculation is almost always asking exactly one of three questions:

- *What sample size do I need?* — given $\delta$, $\sigma$ and a target power $\pi_\text{target}$, find the smallest $n$ that meets it.
- *What power will my planned experiment have?* — given $n$, $\delta$ and $\sigma$, return $\hat\pi$.
- *What is the smallest effect I can hope to detect?* — given $n$, $\sigma$ and $\pi_\text{target}$, return the smallest $\delta$ that meets it.

Each is an optimisation over the same simulation engine. In peppwR they are exposed as the `find` argument to `power_analysis()`, with values `"sample_size"`, `"power"` and `"effect_size"` respectively. The package returns a `peppwr_power` object with print, summary and plot methods.

### Aggregate mode: no pilot data

If pilot data is not available — typical at the grant-proposal stage — the package's *aggregate mode* lets the analyst supply a distribution family and parameters directly, and returns an order-of-magnitude answer. A worked aggregate call:

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

For these parameters the simulation returns $n^\ast = 25$ per group to detect a two-fold effect at $80\%$ power. The output is enough to bound the conversation about how big the experiment needs to be, not to fine-tune it.

```
+----------------------------------+
|                                  |
|       [FIGURE PLACEHOLDER]       |
|                                  |
+----------------------------------+
Figure 7. Aggregate-mode power curve from peppwR. x-axis: n per group.
y-axis: estimated power. Horizontal dashed line at 0.8. Vertical
dashed line at n* = 25, where the curve first crosses the target.
```

### Per-peptide mode: with pilot data

With pilot data, peppwR fits one distribution per peptide via `fit_distributions()`, choosing among gamma, lognormal, normal and inverse Gaussian by AIC and tracking each peptide's missingness rate $r_i$ along the way. The simulation engine then runs the recipe above *independently* for every peptide and aggregates the results, reporting the *peptidome fraction*: the proportion of peptides reaching target power,

$$
P_\text{peptidome}(n, \delta, \pi_\text{target}) \;=\; \frac{1}{K}\sum_{i=1}^{K} \mathbf{1}\!\left[\hat\pi_i(n, \delta) \geq \pi_\text{target}\right],
$$

where $K$ is the number of peptides. The headline figure for a per-peptide run takes the form *"X% of peptides reach 80% power at $n = Y$"*. This is the mode that matches the way a real phosphoproteomics experiment behaves: heterogeneously, peptide by peptide, with no single global power number.

When fitting fails for individual peptides — small samples, extreme values, numerical trouble — peppwR offers three fallback strategies via `on_fit_failure`: drop the peptide and report the count (`"exclude"`, the default); resample empirically from its observed values (`"empirical"`); or substitute a moment-matched lognormal (`"lognormal"`).

### Multiple testing: from $\alpha$ to BH-FDR

Once thousands of tests are running side by side, the Type I error rate compounds. At $\alpha = 0.05$ across $5{,}000$ tests under the null, the expected number of false positives is $250$. That is not a usable workflow.

The standard fix in genomics and proteomics is the Benjamini–Hochberg (BH) procedure for controlling the *false discovery rate* (FDR) — the expected proportion of false positives among the calls that are declared significant. With $m$ tests producing ordered $p$-values $p_{(1)} \leq p_{(2)} \leq \cdots \leq p_{(m)}$, BH finds the largest rank $k^*$ for which

$$
p_{(k^*)} \;\leq\; \frac{k^*\, \alpha}{m},
$$

and rejects every hypothesis with $p \leq p_{(k^*)}$. The threshold for being called significant becomes adaptive: it depends on the whole vector of $p$-values, not on each one in isolation. Peptides effectively compete for the available significance budget.

This has a direct consequence for per-peptide power. A peptide's chance of being declared significant depends not only on its own simulated $p$-value but also on the rest of the peptidome's, and on what fraction of the peptidome is genuinely null. peppwR's `apply_fdr = TRUE` re-estimates per-peptide power *under* the BH procedure: a full simulated peptidome is generated each iteration, BH is applied, and the proportion of true positives among the discoveries is tracked. The proportion of peptides assumed to be under the null is set by `prop_null` (default $0.9$, reflecting the typical proteomics expectation that most peptides do not respond to a given perturbation).

### A worked example: the DDA Arabidopsis dataset

The DDA dataset already used for Figure 5 has $2{,}228$ phosphopeptides measured across two timepoints with three biological replicates per timepoint — a sample size that is entirely typical of discovery phosphoproteomics. Running per-peptide distribution fitting and then per-peptide power analysis at $n = 3$ per group, a two-fold effect, and the Bayes-factor test, produces the picture in Figure 8.

```
+----------------------------------+
|                                  |
|       [FIGURE PLACEHOLDER]       |
|                                  |
+----------------------------------+
Figure 8. Per-peptide power across 2,228 DDA phosphopeptides at n = 3
per group, two-fold effect, Bayes-factor test. x-axis: estimated
power per peptide. y-axis: count of peptides. Vertical line at 0.8.
Most peptides fall below the line.
```

At this sample size, the median per-peptide power is low and only a small fraction of peptides clear the conventional $80\%$ threshold. Asking peppwR the inverse question — what sample size *would* reach $80\%$ power for half of the peptidome at a two-fold effect — returns $n = 15$ per group. That is roughly five times the discovery norm of three, and is the kind of quantitative argument that a power analysis is supposed to produce. The figure does not say the existing experiment is worthless; it does say it is best read as a hypothesis-generating instrument, not as a definitive list of regulated peptides.

The same kind of analysis on the PRM dataset (a targeted run of $285$ phosphopeptides from *Magnaporthe oryzae*, the rice blast fungus, with three biological replicates per timepoint and two technical replicates that we average before analysis) makes the missingness story concrete. The dataset-level MNAR diagnostic returns a clear negative correlation between mean abundance and NA rate. Running the simulation with `include_missingness = TRUE` versus `FALSE` reveals how much of the nominal power is purely an artefact of pretending no values were missing.

### Bridge into Act III

A power calculation tells you whether an experiment is worth running. Once it has been run, the same statistical realities — skew, missingness, multiple testing — return at the analysis stage, and they need to be handled there too. That is what the companion package, `pepdiff`, is for.

---

## Act III: pepdiff

### The same data, now from the other end

Suppose the planning is done. The power analysis has talked the bench team into running, say, fifteen biological replicates per group rather than three; the sample preparation has gone as well as sample preparation ever does; and a spreadsheet of intensities has landed on the analyst's desk. The shape of that spreadsheet is, in the simplest case, one row per peptide per sample, with columns naming the peptide, the gene it maps to, one or more experimental factors (treatment, timepoint, genotype), a biological replicate index, and the measured abundance. The same data quirks that motivated the simulation engine of the previous section — strict positivity, right skew, variance scaling with the mean, MNAR missingness, thousands of simultaneous comparisons — are now the analyst's problem rather than the planner's.

The companion package `pepdiff` is built around the position that those quirks should be handled by the inference machinery rather than papered over before it.

### The pepdiff workflow at a glance

The whole pipeline is a short chain of S3 objects:

$$
\text{CSV} \;\longrightarrow\; \texttt{read\_pepdiff()} \;\longrightarrow\; \texttt{pepdiff\_data} \;\longrightarrow\; \texttt{compare()} \;\longrightarrow\; \texttt{pepdiff\_results} \;\longrightarrow\; \texttt{plot()}.
$$

`read_pepdiff()` takes the path to a long-format CSV plus the names of the columns that play each role (`id`, `gene`, `value`, `factors`, `replicate`) and returns a `pepdiff_data` object with `print`, `summary` and `plot` methods. The `plot` method on a `pepdiff_data` object draws three exploratory views before any model is fitted: a PCA, the per-sample intensity distributions, and a per-peptide missingness summary.

```r
dat <- read_pepdiff(
  file      = "phospho_intensities.csv",
  id        = "peptide",
  gene      = "gene_id",
  value     = "value",
  factors   = c("treatment", "timepoint"),
  replicate = "bio_rep"
)
summary(dat)
plot(dat)
```

```
+----------------------------------+
|                                  |
|       [FIGURE PLACEHOLDER]       |
|                                  |
+----------------------------------+
Figure 10. pepdiff_data exploratory triptych. Left: PCA of samples,
coloured by treatment, shaped by timepoint; replicates from the same
condition should cluster. Centre: per-sample density curves of log
intensity, overlaid; large horizontal shifts flag a normalisation
problem. Right: scatter of per-peptide mean abundance against per-
peptide NA rate; the same MNAR diagnostic introduced in Figure 6.
```

The exploratory plots are not decoration. The MNAR scatter in the right panel is the same diagnostic that peppwR uses at the planning stage.

### The per-peptide gamma GLM

For the default analysis, `compare()` fits, *for each peptide independently*, a generalised linear model (GLM) — a regression in which the response is modelled by an exponential-family distribution and connected to a linear predictor through a chosen *link function*. The family is *gamma* — the same continuous, positive, right-skewed distribution we met when fitting peptide distributions earlier — and the link is the logarithm. For peptide $i$ with intensity $Y_{ij}$ in sample $j$, design vector $\mathbf{x}_j$ (encoding treatment, timepoint, and any interactions), and coefficient vector $\boldsymbol{\beta}_i$,

$$
Y_{ij} \;\sim\; \text{Gamma}(\mu_{ij}, \phi_i), \qquad \log \mu_{ij} \;=\; \mathbf{x}_j^{\top} \boldsymbol{\beta}_i.
$$

(A notational aside: this $\boldsymbol{\beta}_i$ is distinct from the $\beta$ for Type II error rate introduced earlier — the convention forces both symbols, and we keep the GLM coefficient bold and subscripted to mark the difference.) Three things are worth pulling out of that line. First, each peptide gets its own coefficients and its own dispersion $\phi_i$, so the model does not pool variance across peptides — a deliberate choice, given that variance scales with the mean and that different peptides sit at very different abundance levels. Second, the gamma family enforces strict positivity at the level of the likelihood, so the fit does not need a pseudo-count or a log-transform to behave. Third, the log link means that coefficients are directly interpretable as log fold changes: $\exp(\beta_{i,\text{trt}})$ is the fold change of peptide $i$ between treatment and control, with no extra arithmetic.[^3]

[^3]: This is the same reason limma's voom path log-transforms before linear modelling. The gamma-log GLM gets to the same interpretive convenience by routing it through the link function rather than the response.

Contrasts of interest — treatment minus control at a fixed timepoint, the interaction between treatment and timepoint, the marginal effect of treatment averaged over timepoints — are extracted from the fitted model using the `emmeans` package, R's standard tool for estimated marginal means and contrasts. Each contrast yields a per-peptide effect size and a per-peptide $p$-value.

```r
results <- compare(
  dat,
  compare = "treatment",
  ref     = "ctrl",
  method  = "glm"
)
```

The output is a `pepdiff_results` object whose `results` table holds, per peptide and per contrast, the fold change, the log$_2$ fold change, the raw $p$-value, the BH-adjusted FDR, and a Boolean `significant` flag at the conventional FDR $\leq 0.05$ threshold.

### Multiple testing, again

The multiple-testing problem is the same one peppwR had to solve at the planning stage. With $m$ peptides tested simultaneously, raw $p$-values are not interpretable in isolation. pepdiff applies the same Benjamini–Hochberg procedure peppwR used to convert raw $p$-values into FDR-adjusted ones,

$$
p_{(k^*)} \;\leq\; \frac{k^* \alpha}{m},
$$

storing the result in the `fdr` column. When the design has nested structure — e.g., a single $\texttt{compare}$ call that produces contrasts within several timepoints — the correction is applied within each comparison group via `apply_fdr_by_comparison()`, so that a peptide does not have its significance budget eaten by tests from a different stratum.[^4]

[^4]: The choice of correction unit is a real design decision and deserves to be made explicitly. pepdiff's default is to correct within each comparison; pooling all comparisons before correcting is more conservative and is a perfectly defensible alternative.

### When the GLM does not fit

A per-peptide gamma GLM is not universally appropriate. For some peptides the family is wrong: the residuals are too heavy-tailed, the variance does not scale the way gamma expects, or extreme values pull the fit around. pepdiff exposes this directly. After running `compare(..., method = "glm")`, the function `plot_fit_diagnostics()` returns a four-panel summary: the distribution of per-peptide deviances; deviance versus estimated effect size; representative residual plots for individual peptides; and a pooled QQ plot of standardised residuals across the whole peptidome.

```
+----------------------------------+
|                                  |
|       [FIGURE PLACEHOLDER]       |
|                                  |
+----------------------------------+
Figure 11. plot_fit_diagnostics() output. Top-left: histogram of per-
peptide deviances, with the 95th-percentile threshold marked.
Top-right: deviance versus log2 fold change, peptides flagged as
poorly fitted highlighted. Bottom-left: standardised residuals
versus fitted values for three example peptides (high, median, low
deviance). Bottom-right: pooled QQ plot of standardised residuals
across all peptide models. Heavy tails or S-shapes here are the
signal to switch methods.
```

The documented decision rule is conservative: if the proportion of peptides flagged as high-deviance is below roughly $10\text{--}15\%$, and the QQ plot is roughly linear, the gamma GLM result stands. If more than around $15\%$ of peptides are flagged, or the residual plots show systematic curvature, the recommended response is to switch to the *Aligned Rank Transform* (ART) method.

ART is a nonparametric procedure for factorial designs. The recipe is: for each effect of interest, *align* the data by subtracting out all other estimated effects so that the residuals reflect only that effect; *rank* the aligned values; then run an ordinary ANOVA on the ranks. The alignment step is what allows the test to handle interactions, which classical rank-based tests cannot, and the ranking step is what frees the user from having to pick a parametric family. In `pepdiff`,

```r
results_art <- compare(dat, compare = "treatment", ref = "ctrl", method = "art")
```

returns a `pepdiff_results` object of the same shape as the GLM call, with the same `fdr` column and the same downstream plot methods. The price paid for the robustness is interpretive: ART returns ranks-based test statistics rather than fold changes, so for effect-size reporting one still wants the GLM fit (or a separate fold-change summary) alongside the ART significance call.

### Pairwise tests for simpler designs

For the simplest experimental designs — two groups, no factorial structure, small $n$ — neither a GLM nor ART is strictly necessary. `compare(..., method = "pairwise")` runs a single per-peptide two-sample test, with `test` chosen from the same family of tests that peppwR simulates over at the planning stage: the *Wilcoxon* rank-sum test (the documented default for the small-$n$, non-normal, outlier-prone regime that phosphoproteomics typically sits in); a bootstrap-$t$ for four or more replicates with roughly symmetric data; a Bayes-factor $t$-test when the question is about evidence strength rather than rejection (and in particular when "no detectable effect" is a finding one wants to express); and a *rank-products* test, which is sensitive to small effects whose direction is consistent across replicates.

peppwR's per-peptide simulation engine also exposes a `test` argument with the same three options (Wilcoxon, bootstrap-$t$, Bayes-factor $t$). Picking the same test at both stages means the power calculation is estimating power *for the test that will actually be run*, rather than for some convenient closed-form proxy. A power analysis that assumes a $t$-test but is followed by a Wilcoxon analysis is asking and answering different questions.

### Reading the output

For any of the three analysis routes (GLM, ART, pairwise) the post-analysis diagnostics are the same. `plot()` on a `pepdiff_results` object returns a volcano plot, a $p$-value histogram, and a fold-change distribution.

```
+----------------------------------+
|                                  |
|       [FIGURE PLACEHOLDER]       |
|                                  |
+----------------------------------+
Figure 12. Post-analysis triptych from plot(<pepdiff_results>). Left:
volcano plot of -log10(FDR) versus log2 fold change, with the FDR
0.05 line and a fold-change threshold marked; significant peptides
coloured. Centre: histogram of raw p-values; a flat base with a spike
near zero is the expected signature of real effects under a well-
calibrated test. Right: density of log2 fold changes across all
peptides; should be roughly centred near zero unless a global shift
is expected.
```

The $p$-value histogram is the diagnostic that quietly does the most work. A flat distribution with a spike near zero is what a well-calibrated test under a mix of nulls and alternatives is supposed to produce. A U-shape, a spike at one, or a non-flat baseline is a signal that something upstream of the FDR step is mis-specified — a wrong family, an unmodelled batch effect, a normalisation problem — and that the FDR-adjusted column should be read sceptically until the cause is found.

`significant(results)` returns the filtered peptide list, ready to be joined back onto gene annotations or piped into downstream enrichment.

---

## Takeaways

- Power is a property of an experimental design, not of a result. Compute it before, not after.
- Independent biological replicates are what carry statistical leverage; pooling technical replicates as if they were biological ones is $n$-hacking, which is $p$-hacking.
- For high-throughput phosphoproteomics, the closed-form $t$-test power formula does not survive contact with the data. Use simulation.
- peppwR provides that simulation, in aggregate mode for grant-stage estimates and per-peptide mode for realistic, peptidome-level answers, with optional MNAR-aware and BH-FDR-aware variants.
- `pepdiff` provides per-peptide gamma GLMs with a log link for interpretable fold changes, ART as the robust fallback when the GLM diagnostics fail, the same family of pairwise tests that peppwR simulates over, and BH-FDR control end-to-end.
- Match the test at the analysis stage to the test assumed at the planning stage. A Wilcoxon-powered design followed by a $t$-test analysis is asking and answering different questions.
- peppwR is for the experiment you are about to run, `pepdiff` for the one you have already finished.
