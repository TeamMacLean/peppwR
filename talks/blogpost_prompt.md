# Agent prompt: long-form blog draft of the peppwR / pepdiff talk

This file is the explicit prompt handed to a writing agent that drafts a long-form blog post version of a 45-minute conference talk. The purpose of the draft is to let the author **read the prose and judge whether the narrative arc holds together** before any slides are produced. Treat it as a narrative dry run, not a paper and not slides.

When invoking the agent, the calling message should say something like:
> "Read `/Users/macleand/Desktop/peppwR/talks/blogpost_prompt.md` and follow it exactly. Output goes to `/Users/macleand/Desktop/peppwR/talks/draft_blogpost.md`."

---

## Mode

There are two modes. The invoking message states which one applies.

**Mode A — Full draft.** The default. Write the whole blog post end-to-end into `draft_blogpost.md`, overwriting any existing content. Follow every section below.

**Mode B — Append-only Act III replacement.** The author is happy with the existing front-matter, introduction, Act I, and Act II in `draft_blogpost.md` and considers them locked. Your job is to replace ONLY the Act III stub with a full Act III, plus update the closing "What to do on Monday morning" section if symmetry demands it. Do NOT touch any line above the existing `## Act III` heading.

In Mode B:
- Read the locked sections (front-matter through end of Act II) carefully — they are the tonal and narrative ground truth. The new Act III must continue the same voice, the same level of mathematical rigour, the same symbol conventions ($\beta, \delta, \sigma, n, \alpha$ and anything else introduced in Acts I–II), the same figure-placeholder format, the same disposition toward the reader.
- Target ≈120–150 lines for the new Act III (proportional to its ~11 min of talk time relative to Acts I+II).
- The closing "What to do on Monday morning" section can be expanded conservatively to incorporate `pepdiff` alongside `peppwR`, but do not over-rewrite it.
- All other rules in this prompt (tone ban-list, narrative continuity, formula policy, figure placeholders, what-not-to-do) still apply.

---

## Output

- Write to `/Users/macleand/Desktop/peppwR/talks/draft_blogpost.md`. Overwrite if it exists.
- Target ≈500 lines; hard ceiling 550. Markdown.
- Front-matter (title / description / author / date placeholder) at the top.
- **Do not unilaterally rewrite the working title.** The current title is *"Two Weird Tricks…and What To Do When They're Not Enough"*. Keep it. If you have a strong objection, leave it as-is and add a one-line note under the front-matter saying so — do not change the title yourself.

## What to read first

1. `/Users/macleand/Desktop/peppwR/talks/skeleton.md` — the agreed structure. Follow it bullet by bullet. The skeleton is authoritative; this prompt is supporting guidance.
2. `/Users/macleand/Desktop/peppwR/CLAUDE.md` — peppwR package overview (purpose, modes, S3 classes, distributions, tests, FDR-aware mode, missingness handling). Treat as the source of truth for what peppwR actually does. Do not invent features that aren't here or in the paper.
3. `/Users/macleand/Desktop/peppwR/paper/` — the JSS-submitted manuscript. Especially anything that describes the mathematical / simulation framework, the BH-FDR derivation, the AIC selection rule, and the case studies. Use it to pull in rigorous-but-readable equations for Act II. **Pick equations from the paper, don't invent maths.**
4. The source blog post is pasted at the bottom of this prompt. It is the source material for Act I — condense, do not expand. About one-third of the final draft is Act I.

## Structure

Follow `skeleton.md` exactly. Three acts, roughly equal length:

- **Act I** — condensed retelling of "Two Weird Tricks". Lands on the bridge into Act II.
- **Act II** — peppwR. Pull the mathematical / simulation framework from the paper at the right level (the audience is mixed but mostly early-career wet-lab; they have signed up for a bioinformatics talk, so equations are welcome if explained).
- **Act III** — `pepdiff`.

  In **Mode A**: write a SHORT stub (≈30–50 lines max) signposting what will go there later. Do NOT invent details about pepdiff's internals.

  In **Mode B**: write the full Act III (~120–150 lines). Source material is the `pepdiff` documentation site at <https://teammaclean.github.io/pepdiff/> — fetch the homepage and follow links to the Get Started / vignette, Articles, and Reference pages. Pull out:
  - what pepdiff actually does (purpose, scope, where it sits relative to peppwR);
  - the workflow (input data shape → transformations / fits → outputs);
  - the statistical / mathematical framework (model, test, multiple-testing handling) — express it in equations where it helps, using the symbol conventions already established;
  - 2–3 figure placeholders with descriptive captions, matching the convention used in Acts I–II.

  Mirror Act II's structure: motivate by showing data or a problem first, then introduce the machinery. Do NOT invent features that aren't in the docs. If the docs do not state something explicitly, do not assert it — write around it or flag the gap.

## Voice and tone (read this carefully — the previous draft missed)

The target register is *Nature Methods column* or *Royal Statistical Society Significance magazine*: scientifically careful, plainly written, occasionally dry, never matey, never performative.

**Do** write like a careful scientist explaining something to a thoughtful colleague who has not seen the data before. **Do** allow occasional dry wit, preferably in footnotes and at the author's own expense, not at the reader's or the field's. **Do** match the lightly-cheeky-but-controlled register of the source blog post — note that the source's humour lives almost entirely in headers ("Mr. Micawber's rule") and footnotes, while the body is sober.

**Do not** use any of the following — these all appeared in the previous draft and all need to go:

- "vibe", "take it as a vibe", "honestly", "lol", "tbh", "tbf", "ngl"
- "p-hacking with extra steps" or other meme/Rick-and-Morty cadence
- "please stop pretending it is" (header) or other reader-scolding constructions
- "get to feel smug", "ballpark me a number", "Watch them go quiet" — any phrasing that performs cleverness rather than explains
- Challenge-style rhetorical questions aimed at the reader ("Try plugging *that* into pwr.t.test. You can't.")
- "gorgeous and useless", "comprehensively shredded", or other overblown adjective pairs
- Click-baity headers that are not descriptive — every header should accurately name its section

If a sentence would feel out of place in a *Significance* feature, cut it.

## Narrative continuity (the rule)

**Every named concept — every statistical test, distribution family, R function, R package, dataset, biological organism, mathematical symbol — must be introduced before it is used.**

Concretely:

- The first time `pwr.t.test` appears, introduce it as "R's standard closed-form power function in the `pwr` package".
- The first time *Cohen's $D$* appears, define it as a standardised composite of effect size and SD (in units of SD), with the convention that 0.2/0.5/0.8 correspond to small/medium/large.
- The first time AIC appears, give it a one-clause gloss ("a model-selection criterion that rewards fit and penalises complexity").
- The first time a Bayes-factor test appears, give it a one-clause gloss.
- Symbols ($n$, $\sigma$, $\delta$, $\alpha$, $\beta$) must be introduced in Act I and re-used in Act II, not introduced afresh in Act II.
- Do NOT name specific distributions (gamma, lognormal, inverse Gaussian) in the Act I → Act II bridge. They belong inside Act II, motivated by showing the data first.
- Do NOT use specific organism / mutant names from the source blog (*poco1*) or anywhere else without first explaining what they are. The Act I "p is not power" section should use a generic two-sample case–control framing instead. The TSL HR-scoring example in the hook is OK because the skeleton frames it as anonymised and provides enough context.

If you find yourself naming a thing the reader has not yet met, stop and introduce it.

**Do not cite the three-act labels ("Act I", "Act II", "Act III") inside body prose.** They are scaffolding for the author, not narrative referents. When a later section wants to refer back to material from an earlier one, name the subject instead: "the planning stage", "the peppwR power analysis", "when we fitted distributions earlier", etc. The introduction may tee up the act framing; the body must not lean on it.

## No editorialising about design intent

A standing risk: when summarising two related packages or two stages of a workflow, it is tempting to invent a unifying meta-narrative ("the design intent is to mirror X", "what Y really buys you is…", "the deeper purpose of Z is…"). Don't.

- Describe what a package does, not what it was built for. *"`pepdiff` fits a gamma GLM per peptide"* is fine; *"`pepdiff` was built to mirror peppwR's statistical stance"* is editorialising and almost certainly inaccurate.
- Do not assert relationships between packages, papers, or design philosophies that the source material does not state explicitly. Two packages sharing an FDR step is not evidence they were designed as a matched pair.
- If a claim about a package's purpose or intent cannot be sourced from its docs, README, or paper, drop the claim or hedge it explicitly (*"one way to use these together is…"*).
- Throughlines between sections come from showing the same data, the same symbols, the same statistical machinery — not from authorial assertions about what the packages "really" do.

Examples of sentences that would fail this rule:

- *"What pepdiff buys, in one sentence, is the ability to run an analysis whose statistical assumptions match the ones peppwR was using…"*
- *"The design intent is to mirror peppwR."*
- *"This is the part of pepdiff where the design-analysis loop closes most tightly."*

## Formula policy

- Light formulas are welcome in Act I, not banned. Use $\beta$, $1-\beta$, $\delta$, $\sigma$, $n$, $\alpha$ alongside the visual intuition so the reader has earned the vocabulary by the time Act II uses it. One short formula per concept, well-explained.
- In Act II, use the rigorous expressions from the paper:
  - the closed-form $t$-test power $\pi(n,\delta,\sigma,\alpha)$
  - the Monte Carlo power estimator $\hat\pi_i(n,\delta) = \tfrac{1}{n_{\text{sim}}} \sum_j \mathbf{1}[p^{(j)} \leq \alpha]$
  - the AIC selection rule
  - the BH threshold $p_{(k^*)} \leq k^*\alpha/m$
  - the peptidome-fraction summary
- Inline LaTeX (`$...$`) for short, display (`$$...$$`) for longer. Don't over-equate; if a sentence does the work, don't add an equation just because you can.

## Figures

Do NOT generate figures. Wherever a figure belongs, drop a placeholder block:

```
+--------------------------------+
|                                |
|        [FIGURE PLACEHOLDER]    |
|                                |
+--------------------------------+
Figure N. <descriptive caption: what the figure should show, axes, colours, key annotations>
```

Use them liberally where they would appear in the talk — the captions are how the author judges what visuals are needed.

## What NOT to do

- Don't invent peppwR features that aren't in CLAUDE.md or the paper.
- Don't invent pepdiff details at all.
- Don't generate actual figure images.
- Don't exceed 550 lines.
- Don't add a "References" section unless it falls out naturally from footnotes.
- Don't rewrite the title.
- Don't use emojis.
- Don't write summary paragraphs at the end of every section.
- Don't address the reader as "you" in a chiding tone; explain, don't lecture.

---

## Source blog post (Act I source material)

The following is Dan MacLean's blog post "Two Weird Tricks For Improving Your Experimental Designs" (2024-05-17). Condense it for Act I. Do not import it wholesale. Do not import *poco1* without first defining it (better: replace with a generic case–control example).

---

title: "Two Weird Tricks For Improving Your Experimental Designs"
description: "The importance of assessing experimental power and treating replicates properly."
author: Dan MacLean
date: 2024-05-17

The following blog post started out life as a short talk to try and convince PhD students and Post-Docs to think again about the inherent and often deleterious weaknesses in their experimental designs. Having stared at their figures in talks and read them in papers for years I started to see patterns of behaviour that were leading them to make some conclusions that were not as robust as they believed them to be. The talk I came up with was essentially about statistics, which is never an easy sell to molecular biologists, so I gave it a click-baity title and intended to perform a bait and switch, tricking them to coming to a room for some magic answer to a thorny problem, locking the door and then lecturing them. That mostly worked, to the extent that someone even asked for the notes to be put online.

In essence the talk was really about optimising the statistical power of an experiment, with a short deviation about not abusing the p-value by taking it as some indication of reality. The example was a real TSL paper (kept anonymous) that, in good faith, n-hacked itself.

## Two tricks, one parameter

What is the power of an experiment? Definition: "the chance of finding something given that it is there". Formal definition in a test: the likelihood it correctly rejects the null hypothesis of no difference at a given p; equivalently 1 − false negative rate.

## Power is a function of four parameters

Effect size (mean difference), variability (SD), sample size (n), and the chosen p-threshold (the frequency the effect would be seen in a null model). Walk through each visually:
- Large vs small effect at fixed SD
- Large vs small SD at fixed effect
- Large vs small n at fixed effect and SD
- For correlations, R² composites all three.

## p is not power

p is not the probability of a difference. p is closer to a frequency: how often would I get the test sample's mean if I drew at random from a null model based on the control? Build the intuition with a worked example (a generic case–control framing — DO NOT use poco1 / Arabidopsis without defining them), normal distributions overlaid, asking how often a draw from the control gives the treatment mean. No need for closed-form derivation; the intuition carries it.

## Power vs sample size

Standard pwr.t.test curves at Cohen's D = 0.2, 0.5, 0.8. Headline: 80% power needs about 20 reps for a large effect. Much higher than typical bench practice. Most experiments are underpowered.

## Oh. No?

If you got p ≤ 0.05 in a case-control test, the conclusion is still valid even if power was low. Post-hoc power isn't useful for "correcting" a non-significant result — that's p-hacking. High-throughput experiments, though, miss a lot when underpowered.

## What can I do?

Power analysis a priori, with pilot data or literature values. Use it to decide n before running the real experiment.

## Increase sample size — but properly

The easiest lever is n. But independence matters: technical replicates ≠ biological replicates. Counting tech reps as biological reps inflates n, shrinks SD, shrinks p. n-hacking is p-hacking.

## TSL anonymised example

HR scoring (ordinal categorical) treated as numeric, violin-plotted, with technical replicates pooled into the n. The right plot is a count plot with size showing frequency; the right merge is rank-then-average within each independent replicate. Re-analysis: with n-hacked data, 75 of 120 pairwise comparisons significant; with merged independent reps (n=2), 42. The difference is the false positive inflation.

## Two weird tricks

1. Properly merge technical replicates (rank-and-average for ordinal data).
2. Increase independent replicates, guided by a priori power analysis.

---

End of source. Now write the blog post following `skeleton.md` and the constraints above.
