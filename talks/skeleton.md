# Talk skeleton: "Two Weird Tricks…and What To Do When They're Not Enough"

**Slot:** 45 min total (≈35 min talk + 10 min Q&A)
**Audience:** mixed, mostly early-career wet-lab. Stats / R / maths welcome to the level in the manuscript.
**Format:** Quarto reveal.js (to be drafted later).
**Goal:** Convince the audience that (a) power analysis is necessary, (b) high throughput proteomics can be approached with peppwR to estimate power/do experimental design c) pepdiff can be used for powerful analysis with well designed experiments.

---

## Act I — Power analysis, the general case (~12 min)

Condensed retelling of the "Two Weird Tricks" blog post.
Use some of the mathematical framework in the peppwr paper to allow reader to get a hook when we use something similar in the Act II. 
- **Hook.** Anonymised TSL HR-scoring example: a published, well-intentioned paper that n-hacked itself into too many significant calls. Sets up emotional stakes before any equations.
- **What is power?** Plain-English definition first, then the formal one (1 − false negative rate). Introduce $\beta$ as the false negative rate and $1-\beta$ as power so the symbol is earned before Act II uses it.
- **The four parameters.** Effect size ($\delta$), variability ($\sigma$), sample size ($n$), and significance threshold ($\alpha$). Visual intuition for each (large/small effect, large/small SD, large/small n) paired with the symbols — these four letters reappear in Act II, so they need to feel familiar by the end of Act I, not new.
- **$p$ is not power.** How $p$ is actually computed: null model + "how often does this mean show up by chance?". Use a generic two-sample case–control framing (treatment vs control, continuous measurement) — do NOT introduce a specific organism or mutant by name unless it is first defined; the source blog's *poco1* example assumed reader background we don't have.
- **Power vs $n$ curve.** The standard closed-form two-sample $t$-test result, plotted via R's `pwr.t.test` (introduce it on first mention as "R's standard closed-form power function"). Use Cohen's $D$ at 0.2 / 0.5 / 0.8, defining $D$ on first use as a standardised composite of effect size and SD (in units of SD). Headline: 80% power needs ~20 reps per group even for a large effect, far above the bench norm of three.
- **Is my existing experiment useless?** Probably not, but high-throughput experiments miss a lot when underpowered.
- **The two weird tricks.** (1) Properly merge technical replicates (rank-and-average for ordinal data). (2) Increase true independent replicates. $n$-hacking is $p$-hacking.
- **Bridge into Act II.** "All of that assumed one $t$-test on roughly normal data. What happens when the experiment is five thousand simultaneous tests on data that are skewed, span orders of magnitude, and are riddled with non-random missingness? The closed-form recipe falls over. You need simulation." Do NOT name specific distributions (gamma, lognormal, …) in the bridge — they belong inside Act II, motivated by showing the data first.

## Act II — peppwR (~12 min)

Why phosphoproteomics breaks the textbook recipe, and what peppwR does about it. Use the rigorous descriptions in the peppwR paper for the equations, but introduce each concept by showing the data first.

- **Show the data first.** Open Act II with a real phosphopeptide intensity histogram from the DDA or PRM sample data. Ask the question aloud: does this look like a bell curve? It doesn't — it is strictly positive, has a long right tail, and (often) a spike of zeros. Only *now*, having motivated the need, name the candidate families that fit this shape: gamma, lognormal, inverse Gaussian. Introduce AIC as the rule for choosing among them when it appears.
- **Why closed-form power doesn't work.** Walk through the assumptions of the standard $t$-test power formula ($\pi(n,\delta,\sigma,\alpha)$) one at a time — normality, homogeneous variance, complete data, single test — and show that the phosphoproteomics data we just looked at violates each one. The case against the textbook recipe should land *because we already saw the data*, not because we asserted it.
- **Missingness, and why it isn't random.** Before naming MNAR, show a scatter of mean abundance vs missingness rate from the sample data. The negative trend is the motivation. Then name the phenomenon (MNAR — missing not at random) and explain why pretending it is random (MCAR) inflates power estimates.
- **Simulation as the way out.** The whiteboard recipe: fit a distribution to pilot data → draw simulated control and treatment samples → run the actual test → count how often $p \leq \alpha$ → that fraction is your power. Then write it formally as $\hat\pi_i(n,\delta) = \tfrac{1}{n_{\text{sim}}} \sum_j \mathbf{1}[p^{(j)} \leq \alpha]$, using the $\alpha$, $\delta$ and $n$ symbols already introduced in Act I.
- **The three questions framing.** A researcher coming to power analysis is asking exactly one of: "What $n$ do I need?", "What's my power?", "What's the smallest effect I can detect?" peppwR exposes all three via `find = "sample_size" / "power" / "effect_size"`. State each as an optimisation over the same simulation engine.
- **Aggregate mode.** No pilot data — supply distribution + parameters, get an order-of-magnitude answer. Useful for grant proposals or sanity checks.
- **Per-peptide mode.** With pilot data, fit one distribution per peptide; report power as "X% of peptides reach 80% power at $n = Y$". This is the mode that matches how the experiment actually behaves.
- **FDR-aware mode.** Once the reader accepts that we have thousands of tests, motivate FDR control: at $\alpha = 0.05$ across 5,000 tests, 250 chance hits is unacceptable. Introduce Benjamini–Hochberg with the formula from the paper. Then explain why per-peptide power must be re-estimated under whole-peptidome BH (peptides compete for slots), and how `apply_fdr = TRUE` and `prop_null` parameterise that.
- **Worked example.** A `peppwr_power` output + a peptide-threshold curve / heatmap from the sample DDA or PRM data. Name the dataset and its dimensions when introducing it.
- **Bridge into Act III.** "Power analysis tells you whether an experiment is worth running. Once you've run it, you still have to call differential abundance from the data — and the same statistical reality (skew, missingness, multiple testing) applies. That's `pepdiff`."

## Act III — pepdiff (~11 min)

*Placeholder.* Fill in after Acts I & II are firm. Should cover what pepdiff does, how it differs from limma / MSstats, a worked example, and a closing "what to do Monday morning" takeaway.

---

## Takeaways (closing slide candidates)

- Power analysis before, never after.
- Merge technical replicates properly; n-hacking is p-hacking.
- For high-throughput phosphoproteomics, use simulation (peppwR), not closed-form.
- Design and analysis are one loop: peppwR → experiment → pepdiff.

## Open questions / decisions to revisit

- Final title.
- Whether to keep the TSL anonymised example or swap for a phosphoproteomics hook.
- How much R code to show on slides vs. point at the docs.
- pepdiff section contents (see above).
