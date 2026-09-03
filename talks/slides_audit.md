# Slide audit — `talks/slides.qmd`

(Read-only audit. Severity: LEAK > FORMAT > CHECK. `[NEEDS DRAWING: …]` blocks are intentional scaffolding and are not flagged.)

## Title slide (lines 1–29)

- Clean. The HTML comment (lines 15–23) is invisible. The two hex-logo image embeds on line 25 render onto the auto-generated title slide as intended. The `.notes` block (lines 27–29) is speaker-only. No leaks, format issues, or check items.

## Slide 2 — "Roadmap" (lines 31–43)

- [LEAK] line 39: `[On-screen layout: three equal-weight boxes side by side, no numbering. See` `slides_outline.md` `slide 2.]` — pure outline-pointer directed at the author/designer. Will render on the slide as a literal bracketed paragraph telling the audience what the slide is supposed to look like.
- [CHECK] line 37: the italic "the classical recipe and where it runs out · …" line is meant as a caption to the `.r-stack` box above, but sits as a normal italic paragraph; worth confirming it does not overflow or look orphaned beneath the stack.

## Slide 3 — "Power analysis {.section-divider}" (lines 45–53)

- [LEAK] line 49: `[On-screen: full-bleed section title. See` `slides_outline.md` `slide 3.]` — outline pointer that will render on the section-divider slide as a bracketed instruction.

## Slide 4 — "What 'power' actually means" (lines 55–63)

- Clean. Blockquote, prose paragraph, and `.notes` only.

## Slide 5 — "Power = $1 - \beta$" (lines 65–75)

- Clean. Display equation and prose paragraph; `Prose lines 22–26.` reference is inside `.notes` and so is fine.

## Slide 6 — horizontal rule + "The four levers" (lines 77–87)

- The `---` on line 77 is a horizontal rule / slide break. The slide content is the `[NEEDS DRAWING: …]` block (intentional scaffolding, not flagged) plus the italic tag-line on line 83.
- [CHECK] line 83: `*Four levers. These four letters return later.*` — italic standalone line that may render as a stray caption beneath the placeholder. Confirm intent.

## Slide 7 — "Two curves, low power and high power" (lines 89–97)

- [CHECK] line 93: `*Pull apart, narrow down, or sample more — all of them shrink the overlap.*` — italic standalone tag-line that may look like a TODO/caption rather than body. Worth a glance; same pattern as previous slide.

## Slide 8 — "$p$ is not power" (lines 99–105)

- Clean. Blockquote only; `.notes` is speaker-only.

## Slide 9 — "What $p$ actually is" (lines 107–115)

- [CHECK] line 111: `*How often would the control alone produce a mean this extreme?*` — italic caption-style line following the `[NEEDS DRAWING: …]` placeholder. Probably intentional, but confirm.

## Slide 10 — "Power is a property of a design" (lines 117–125)

- [CHECK] line 119 + 121: `*Post-hoc power on the observed effect adds no new information.*` followed by a plain `Hoenig & Heisey (2001).` line. The italic line plus bare citation will render as two adjacent paragraphs; confirm this is the intended layout rather than a citation belonging in `.notes`.

## Slide 11 — "Cohen's $D$" (lines 127–141)

- [CHECK] line 137: `*Mean difference, in units of SD. Rules of thumb, not laws.*` — italic standalone caption beneath the table. Likely fine; same pattern recurs.

## Slide 12 — "The classical power curve" (lines 143–149)

- Clean. Only `[NEEDS DRAWING: …]` placeholder plus `.notes`.

## Slide 13 — "Is my old experiment worthless?" (lines 151–161)

- Clean. Bolded answer, two paragraphs, `.notes`.

## Slide 14 — "Trick one: merge technical replicates honestly" (lines 163–171)

- Clean. Placeholder is intentional; the bold line on 167 is body content for the slide.

## Slide 15 — "Ordinal scores need a different merge" (lines 173–181)

- [CHECK] line 175: `*For ordinal scores, rank first — then average.*` — italic caption-style line that precedes the `[NEEDS DRAWING: …]` placeholder. May look like a stray line above the figure scaffolding; confirm intent.

## Slide 16 — "Trick two: decide $n$ before the experiment" (lines 183–194)

- [CHECK] line 190: `*Run the power calculation before you spend the money.*` — italic tag-line below the numbered list; confirm placement vs. moving to `.notes`.

## Slide 17 — "A real anonymous example" (lines 196–206)

- [CHECK] line 198: `*HR scoring, ordinal 0–6, multiple genotype pairs.*` — italic standalone caption directly under the h2; reads as an intentional subtitle but worth confirming.

## Slide 18 — "What the re-analysis showed" (lines 208–218)

- [CHECK] line 214: `*Rank-and-average within biological replicate, then use biological $n$ only.*` — italic standalone line after a bolded headline. May render as a stray methodological caption; confirm.

## Slide 19 — "When the textbook recipe runs out" (lines 220–236)

- [FORMAT] lines 224–230: the `::: {layout-ncol=3}` fence mixes two `[NEEDS DRAWING: …]` placeholder paragraphs with a real image `![](figures/density_overlay.png)`. The placeholders are intentional, but the three-column layout will render two literal bracketed text columns next to a single real figure column. The placeholders themselves are not flagged, but the layout fence will likely look messy in mixed-content form — worth verifying that the layout is intended to display bracketed text columns rather than three figures.
- [CHECK] line 222: `*All of that assumed one $t$-test on normal, complete data.*` — italic standalone line immediately below the h2; same caption-style pattern.

## Slide 20 — horizontal rule + "peppwR {.section-divider}" (lines 238–250)

- The `---` on line 238 is a slide-break before the section divider.
- [CHECK] lines 244 + 246: two consecutive italic standalone lines (`*Power analysis for high-throughput, skewed, missing data.*` and `*Simulation in place of the closed form.*`) beneath the hex image. Renders as two italic paragraphs; confirm whether one was meant as subtitle and the other as note.

## Slide 21 — "Start with the data, not the family" (lines 252–258)

- Clean. Just an image and `.notes`.

## Slide 22 — "The four assumptions of the closed form" (lines 260–271)

- Clean. Display equation + bullets + `.notes`.

## Slide 23 — "All four assumptions fail at once" (lines 273–282)

- Clean. Bullets + `.notes`.

## Slide 24 — "Candidate distribution families" (lines 284–292)

- Clean.

## Slide 25 — "Picking between families: AIC" (lines 294–302)

- [CHECK] line 298: `*$L_k$ is the maximised likelihood of model $k$, $p_k$ its number of parameters. Smallest AIC wins.*` — italic mathematical caption directly under the display equation. Likely intentional, but italics with inline math may look mixed up against the formula; worth a glance.

## Slide 26 — "Fitting per peptide: `fit_distributions()`" (lines 304–319)

- [CHECK] line 315: `*Returns a `peppwr_fits` S3 object with `print`, `plot`, `summary` methods.*` — italic standalone line under the code block. Likely intentional caption; confirm.

## Slide 27 — "Missingness, and why it is not random" (lines 321–330)

- Clean. The HTML comment on line 325 is invisible. Just image + `.notes`.

## Slide 28 — "MNAR vs MCAR" (lines 332–341)

- [CHECK] line 337: `*Pretending MCAR overstates what the experiment can detect.*` — italic tag-line beneath the two bullets; same caption pattern.

## Slide 29 — "peppwR detects MNAR at the dataset level" (lines 343–351)

- Clean.

## Slide 30 — "The simulation recipe, one peptide at a time" (lines 353–363)

- Clean. Numbered list + `.notes`.

## Slide 31 — "Monte Carlo power, formalised" (lines 365–373)

- [CHECK] line 369: `*$\mathbf{1}[\cdot]$ is the indicator function. Same $n, \delta, \alpha$ as before; $\sigma$ enters implicitly through the fitted distribution.*` — italic standalone notation-caption under the display equation. Same pattern as slide 25.

## Slide 32 — "Three tests inside the engine" (lines 375–383)

- Clean.

## Slide 33 — "Three questions, one engine" (lines 385–393)

- Clean.

## Slide 34 — "Aggregate mode: no pilot data" (lines 395–412)

- [CHECK] line 408: `*Returns $n^\ast = 25$ per group.*` — italic standalone line after the code block. Same caption pattern.

## Slide 35 — "Aggregate-mode power curve" (lines 414–420)

- Clean. Just image + `.notes`.

## Slide 36 — "Per-peptide mode: with pilot data" (lines 422–432)

- [CHECK] line 424: `*One fit per peptide, one simulation per peptide.*` — italic subtitle directly under the h2; same caption-pattern.

## Slide 37 — "The peptidome fraction" (lines 434–442)

- [CHECK] line 438: `*$K$ is the number of peptides. Headline format: "X% of peptides reach 80% power at $n = Y$".*` — italic standalone notation-caption under the equation.

## Slide 38 — "When per-peptide fitting fails" (lines 444–454)

- [CHECK] line 450: `*Configured via `on_fit_failure`.*` — short italic standalone line; same caption pattern.

## Slide 39 — "Multiple testing, again" (lines 456–464)

- [CHECK] line 460: `*Not a usable workflow.*` — italic standalone tag-line beneath the bold headline.
- [CHECK] h2 title "Multiple testing, again" recurs verbatim on slide 49 (lines 648–656) in the pepdiff section. Likely intentional (mirroring section structure), but worth confirming the audience won't be confused.

## Slide 40 — "Benjamini–Hochberg" (lines 466–474)

- [CHECK] line 470: `*Ordered $p$-values $p_{(1)} \leq \cdots \leq p_{(m)}$. Find the largest rank $k^*$ satisfying the inequality; reject every hypothesis with $p \leq p_{(k^*)}$.*` — italic standalone notation-caption under the display equation. Same pattern.

## Slide 41 — "FDR-aware power: `apply_fdr = TRUE`" (lines 476–489)

- [CHECK] line 478: `*Power under BH, not in isolation.*` and line 485: `*`prop_null = 0.9` by default — most peptides do not respond.*` — two italic standalone caption-lines bracketing the numbered list. Confirm whether both belong on the slide or one should be a note.

## Slide 42 — "Worked example: DDA *Arabidopsis*" (lines 491–501)

- Clean. Bold headline + three bullets + `.notes`.

## Slide 43 — "Per-peptide power distribution" (lines 503–509)

- Clean. Image + `.notes`.

## Slide 44 — "What $n$ would reach half-peptidome 80% power?" (lines 511–519)

- [CHECK] line 515: `*Roughly five times the discovery norm of three.*` — italic standalone tag-line beneath the bolded answer.

## Slide 45 — "PRM case study: missingness on, off" (lines 521–531)

- Clean.

## Slide 46 — horizontal rule + "pepdiff {.section-divider}" (lines 533–545)

- The `---` on line 533 is a slide break.
- [CHECK] lines 539 + 541: two consecutive italic standalone lines (`*The same data, from the other end.*` and `*Differential abundance for phosphopeptides.*`) under the hex image — mirrors the same double-italic pattern as the peppwR divider (slide 20). Confirm intent.

## Slide 47 — "Input data shape" (lines 547–558)

- [FORMAT] line 554: the final table row is literal ellipses (`| ... | ... | ... | ... | ... | ... |`) — renders as a row of dot triplets in the table, which is a typographic placeholder. Common convention, but worth confirming the rendered output reads as "etc." rather than as malformed data.

## Slide 48 — "The workflow at a glance" (lines 560–566)

- Clean. Single inline display-math pipeline + `.notes`.

## Slide 49 — "Reading the data" (lines 568–585)

- [CHECK] line 581: `*Returns a `pepdiff_data` S3 object with `print`, `summary`, `plot` methods.*` — italic standalone caption beneath the code block. Same recurring pattern as slide 26.

## Slide 50 — "Exploratory triptych" (lines 587–595)

- [CHECK] line 591: `*`plot()` on a `pepdiff_data` object — diagnostics before any model is fitted.*` — italic standalone caption beneath the image.

## Slide 51 — "The per-peptide gamma GLM" (lines 597–607)

- [CHECK] line 603: `*Peptide $i$, sample $j$. Design vector $\mathbf{x}_j$; coefficient vector $\boldsymbol{\beta}_i$; dispersion $\phi_i$.*` — italic standalone notation-caption beneath two display equations. Same pattern as slide 31.

## Slide 52 — "Why gamma, why log link" (lines 609–617)

- Clean.

## Slide 53 — "Contrasts via `emmeans`" (lines 619–629)

- [CHECK] line 625: `*Extracted via the `emmeans` package — R's standard tool for estimated marginal means and contrasts.*` — italic standalone caption-line beneath the bullets.

## Slide 54 — "Running the comparison" (lines 631–646)

- [CHECK] line 642: `*Returns a `pepdiff_results` object.*` — short italic caption beneath the code block. Same recurring pattern.

## Slide 55 — "Multiple testing, again" (lines 648–656)

- [CHECK] line 652: `*Stored in the `fdr` column. Applied within each comparison group via `apply_fdr_by_comparison()` when designs are nested.*` — italic standalone notation-caption under the equation.
- [CHECK] h2 title repeats the earlier slide 39 title verbatim. Likely intentional mirroring; confirm acceptability for the audience.

## Slide 56 — "When the GLM does not fit: diagnostics" (lines 658–664)

- Clean.

## Slide 57 — "The decision rule" (lines 666–675)

- Clean. Bold "When to switch." line + two bullets + `.notes`.

## Slide 58 — "Aligned Rank Transform (ART)" (lines 677–687)

- [CHECK] line 683: `*Handles interactions, which classical rank-based tests cannot. Returns rank-based statistics, not fold changes.*` — italic standalone caption beneath the numbered list.

## Slide 59 — "Running ART" (lines 689–704)

- [CHECK] line 700: `*Same `pepdiff_results` shape, same `fdr` column, same downstream plot methods.*` — italic standalone caption beneath the code block. Same pattern as the other code-block slides.

## Slide 60 — "Pairwise tests for simpler designs" (lines 706–717)

- [CHECK] line 708: `method = "pairwise"` is a bare inline-code line directly under the h2 and is the only body line above the bullet list. Renders as a one-line paragraph; confirm it isn't expected to be inline with the title.
- [CHECK] this slide has four bullets, each with a substantial description; combined with the lead-in line it may approach the vertical limit of the reveal.js viewport at default sizing. Worth a sanity check on the rendered slide.

## Slide 61 — "Match the test at both stages" (lines 719–732)

- [CHECK] line 728: `*Power for the test you will actually run.*` — italic standalone tag-line beneath the table. Same caption pattern.

## Slide 62 — "Reading the output" (lines 734–740)

- Clean. Image + `.notes`.

## Slide 63 — "The $p$-value histogram does the diagnostic work" (lines 742–750)

- Clean. Three bullets + `.notes`.

## Slide 64 — "Filtered output" (lines 752–762)

- [CHECK] line 758: `*Returns the filtered peptide list, ready to join onto gene annotations or pipe into enrichment.*` — italic caption-line under the code block. Same recurring pattern.

## Slide 65 — horizontal rule + "Takeaways" (lines 764–777)

- [CHECK] line 766 onward: this slide carries six dense bullets, three of which contain semicolons separating sub-clauses. At typical reveal.js font sizes the slide may overflow vertically; worth a render-check.

## Slide 66 — "One sentence each" (lines 779–786)

- Clean.

## Slide 67 — "Questions?" (lines 788–796)

- Clean. Two URL lines + `.notes`.

## Summary

- LEAKs: 2
- FORMATs: 2
- CHECKs: 32

**Top three impactful fixes:**

1. Remove or move-to-`.notes` the two literal outline pointers on lines 39 and 49 — they are the only true leaks and will read as authorial scaffolding on two early, attention-defining slides (the Roadmap and the first section divider).
2. Decide a house style for the recurring italic standalone caption/tag-line beneath h2s, code blocks, equations, and figures, then apply it consistently. Twenty-plus slides use this pattern; right now each one is borderline-ambiguous between "intentional subtitle" and "stray TODO". A single decision (keep all / drop all / convert to a subtitle attribute) closes most of the CHECK items in one stroke.
3. Audit the `::: {layout-ncol=3}` fence on slide 19 (lines 224–230) — mixing two `[NEEDS DRAWING: …]` placeholders with one real figure in a three-column layout will render unevenly until the placeholders are replaced; meanwhile the slide looks like a debug view.

**Patterns observed:**

- Italic standalone lines used as captions/tag-lines/notation-keys appear under nearly every code block, equation, figure, and bulleted list. Pattern is consistent enough to be intentional house style, but each instance is individually borderline; treat as a single batch decision.
- `Prose lines NNN–NNN.` back-references appear at the end of almost every `.notes` block. All are correctly inside `.notes` and never leak onto a rendered slide — they are not flagged.
- Section dividers (`{.section-divider}` slides 3, 20, 46) all share the same two-italic-caption shape under a hex image / title. Either intentional or copy-paste — confirm.
- Repeated h2 title "Multiple testing, again" appears on slides 39 and 55, deliberately mirroring planning vs analysis stages. Flagged as CHECK rather than LEAK/FORMAT because the repetition reads as deliberate parallelism.
- `[NEEDS DRAWING: …]` placeholders are concentrated in the first half of the deck (slides 6, 7, 9, 12, 14, 15, 17, 18, 19); the second half is image-complete. Author is aware; not flagged.
