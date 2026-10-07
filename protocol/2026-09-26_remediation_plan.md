# Remediation plan and venue decision brief

Input: `protocol/2026-09-26_ars_peer_review_decision.md` (5-reviewer panel, unanimous Genome
Biology, zero support for Nature Communications as submitted).
Purpose: cost every required fix, separate what GB needs from what NC needs in addition, and isolate
the two decisions that are the author's to make.

---

## 0. Three facts established before planning

**F1 — The headline 0.663 survives the seed problem; the ladder's 150-nucleus point does not.**
Both the main equalised analysis and the ladder pool Mantel–Haenszel strata over cell type × seed,
and the three seeds are re-downsamplings of the same nuclei, so they are not independent strata.
The consequence differs by analysis because the number of genuinely independent strata differs.

The main analysis has 25 cell types. Widening its standard error by √3 to absorb the seed
duplication gives 0.663 → CI (0.511, 0.861); it still excludes 1. Only a full ×3 widening
(0.422, 1.042) would break it, which is the pessimistic bound, not the expected correction.

The ladder's 150-nucleus point has 2 cell types. Reconstructing the 2×2 tables from
`primary_L3_nucleus_ladder.csv` and pooling per seed instead of across seeds:

| pooling | OR | 95% CI | excludes 1 |
|---|---|---|---|
| all six strata as independent (published approach) | 0.712* | 0.527, 0.962 | yes |
| seed 0 only (2 cell types) | 0.972 | 0.611, 1.547 | no |
| seed 1 only | 0.643 | 0.380, 1.088 | no |
| seed 2 only | 0.548 | 0.303, 0.991 | yes |

\*my reconstruction of the published 0.663 (0.486–0.905); the small offset is a background-set
definition difference, the direction and magnitude agree. Between-seed SD of log-OR is 0.296 — large.
The whole point rests on **40 promoter-proximal links** across all six strata (11, 9, 5, 6, 4, 5).

So: "the two primary-cohort cell types pool to 0.663 (0.486 to 0.905)" at 150 nuclei is
significant only because three re-draws of the same nuclei were counted as independent. Two of
three seeds individually cross 1.

**F2 — The window recompute needs the four reference datasets re-downloaded, and a flag that does not exist.**
No `.bedpe` file remains on disk; only the four derived CSVs. `p5_15_promoter_OR_external.py`
requires `--bedpe`, and it *derives* the window from the data (`window = ceil(dist.max()/1e5)*1e5`,
line 104) — there is no `--window` override. This is also why the 6.591 value at ±1.7 Mb is a
hardcoded literal attributed to a flag the script never implemented. The four datasets are public
10x demonstration data (`human_brain_3k`, `pbmc_granulocyte_sorted_3k`,
`pbmc_granulocyte_sorted_10k`, `lymph_node_lymphoma_14k`), so this is a re-download plus a small
code change, not lost work.

**F3 — ALS26 needs a decision, not a wording fix.**
`Supplementary_Table_3a` (expression) covers 27 donors; `Supplementary_Table_2a` (composition)
covers 26. The extra donor is ALS26, which the manuscript's stated inclusion rule — 26 of 27
cross-chip donors "yielded at least 50 nuclei on both" — excludes. Either the expression fit should
be recomputed on 26 donors, or the inclusion rule is per cell type rather than per donor and the
Methods must say so. This is upstream of the exponent, so resolve it before W2.

---

## 1. Work packages

Effort is my own estimate of focused working time, excluding your review.
"Venue" says which submission the item is required for.

### Group A — editorial, no computation

| # | Item | Fix | Files | Effort | Venue |
|---|---|---|---|---|---|
| A1 | "309 replicate pairs" is a units error | 309 = donor × cell-type observations from 27 pairs. Restate; fix Fig. 6b legend | 3 manuscripts | 20 m | both |
| A2 | "these 26 pairs underlie every replicate analysis" false for Layer 2 | depends on F3 | 3 manuscripts + Methods | 20 m | both |
| A3 | Figure 5 never cited in main text | add a call-out in the diagnostic Results subsection; note Fig. 5a shading changes under W1 | 3 manuscripts | 15 m | both |
| A4 | Table 1 "362 samples" vs 357 elsewhere | check the resource's own count vs the deposited matrix; disambiguate both | 3 manuscripts | 15 m | both |
| A5 | Abstract: "reproduced" attached to depletion | → "the absence of promoter enrichment reproduced"; must be word-neutral in NC (150/150) | NC + full abstract | 30 m | both |
| A6 | κ "does not transfer" vs "comparable, 4.28 against 3.90" | one sentence: not established as universal, agreed in two cortex cohorts, varies more between cell types | 3 manuscripts | 15 m | both |
| A7 | Table 3 legend and Methods misdescribe the top two ladder rungs | mark n = 5,000 and 7,500 as one cell type; oligodendrocytes stop at 2,400 | 3 manuscripts + SI Fig. 1 legend | 20 m | both |
| A8 | Crossing points quoted to 3 s.f. while called order-of-magnitude | quote as ≈560 / ≈2,700, or drop the 1,340 / 4,540 extrapolations | 3 manuscripts | 15 m | both |
| A9 | Non sequitur: composition premise, expression conclusion | split into two statements | 3 manuscripts | 5 m | both |
| A10 | 4.79 denotes two unrelated quantities | disambiguate the expression coefficient from the 2,400-nucleus OR 4.791 | 3 manuscripts | 10 m | both |
| A11 | PsychAD modality never stated | state it is hashed snRNA-seq; 1,209 of 1,235 composition pairs are snRNA-seq under a "multiome" title | 3 manuscripts | 20 m | both |
| A12 | Replicate-scarcity framing overstated by the paper's own 1,235 pairs | state the narrow claim narrowly (cross-batch *multiome* replicates are rare) and present the 1,235 pairs as the answer | 3 manuscripts | 20 m | both |

**Group A total ≈ 3.5 h.** No numbers change. Nothing here is at risk.

### Group B — computation on data in hand or re-downloadable

| # | Item | Work | Risk to conclusions | Effort | Venue |
|---|---|---|---|---|---|
| B1 | Window mismatch (V1) | re-download 4 public bedpe; add `--window` to `p5_15`; recompute the four reference ORs at ±500 kb; update Introduction range, Fig. 5a shading, Discussion, Supp. Table 5, Methods line 515 | reference range will **fall** from 2.45–4.01; the gap against 0.663 narrows by an unknown amount. Ladder unaffected | 1–2 d | both |
| B2 | Exponent support-sensitivity (V2) | resolve F3; refit at min-n thresholds; **recompute the matched multinomial null on each restricted support**; deposit the per-pair motor-cortex null | if the exponent departs from its matched null under restriction, the Layer 2 headline ("tracks sampling; only the coefficient is dataset-specific") must be restated. Coefficient 4.79 → 3.46 (≥10) or 2.87 (≥25) propagates into every \|log₂FC\| threshold, Fig. 6b and the design tables | 1–2 d | both |
| B3 | Ladder seed handling (F1) | re-derive the 150-nucleus interval with seeds as repeated measures; re-check the main analysis the same way | the ladder's 150-point loses significance; the main 0.663 is expected to hold at ≈(0.51, 0.86). Headline claim becomes **absence of enrichment**, aligning it with what the external cohort already supported | 0.5–1 d | both |
| B4 | SI Fig. 3 mislabelling (V3) | replot from `OR` vs `OR_varpeaks` at fixed depth 5,564; write the missing generating script for `nabec_L3_convention_sensitivity.csv`; restate the Discussion sentence as three pipeline differences plus depth | **strengthens** the paper: the true convention effect pushes the ratio further below 1 and never across it. The 0.66–1.5 band disappears | 0.5 d | both |
| B5 | κ has no interval; κ is per cell type | report κ CI (≈4.27 [3.03, 5.58]); add per-cell-type κ as a numbered supplementary table (range 1.32–7.44, rising with abundance); give the design cost as a range (38,461 at the median type to 383,037 at oligodendrocytes, against the single 126,670 now quoted) | reframes the design table from a number to a range. More defensible, and a stronger version of your own thesis | 0.5 d | both |
| B6 | No heterogeneity statistic for "incompatible" exponents | report the donor-clustered bootstrap difference (+0.452, 95% CI 0.153–0.801, P ≈ 0.0015) | none — the claim **survives**; it currently just rests on interval inspection | 2 h | both |
| B7 | Hardcoded values falsify "every figure regenerates" | remove literals at `fig2_fig3_ggplot2.R:120`, `fig6_ggplot2.R:79`, `fig1_fig4_fig5_ggplot2.R:257`; make all figures read deposited CSVs; fix the Chinese column header in a deposited table | none | 0.5 d | both |

**Group B total ≈ 5–7 d.** B1, B2 and B3 change published numbers; B4, B5, B6 improve the paper's
position; B7 is hygiene the cover letter currently misstates.

### Group C — new analysis, required for Nature Communications only

| # | Item | Work | Effort | Why |
|---|---|---|---|---|
| C1 | External linker control | run an aggregation-based linker (Cicero or ArchR/Signac `LinkPeaks`) at 150 nuclei, **and/or** run your Pearson-on-log1p linker against one published `feature_linkage.bedpe` at matched cell number and window, reporting concordance | 1–2 wk | Your linker is per-nucleus Pearson on log1p counts with no cell aggregation — precisely what Cicero and ArchR aggregate away because sparsity defeats correlation at low cell number. All five current controls hold the linker fixed, so none excludes "this is a property of an unaggregated linker". R2 rates this **the** condition for an NC-grade Layer 3 claim |
| C2 | Package the estimator | turn floor estimation into a runnable tool (input: two replicate libraries + labels; output: κ, expression floor, design table) with tests and a worked example | 1–2 wk | The EIC's condition. It converts the paper from "a caution plus non-transferable constants" into a portable artefact, which is what Squair / scPower / scCODA actually have in common. Without it the precedent set argues *against* NC |
| C3 | Within-cohort nested decomposition | use the 12 chips / 32 wells: same well vs different well same chip vs different chip, region and protocol held fixed | 3–5 d | Replaces the cross-cohort composition decomposition — confounded by region and hashing, and conceded as unidentified — with a controlled contrast from data already in hand. Optional for GB, valuable for NC |
| C4 | Literature repair | add Cicero, ArchR, Signac/LinkPeaks, SCENT, SCARlink, Enhlink (the last two are in your own novelty audit) plus CRISPR enhancer–gene ground truth | 0.5 d | Required for both; costs nothing in GB, costs words in NC |

---

## 2. The NC word budget, quantified

The NC body is at 4,991 / 5,000 including Methods. Required additions:

| Addition | words |
|---|---|
| B2 exponent sensitivity + restricted-support null | ~80 |
| B5 κ interval + per-cell-type heterogeneity | ~70 |
| B4 convention restatement (three differences plus depth) | ~40 |
| B3 seed handling and the revised interval | ~40 |
| A7 ladder strata caveat | ~30 |
| A11 PsychAD modality | ~25 |
| A12 narrowed scarcity claim | ~15 |
| **total** | **~300** |

So NC needs roughly 300 words cut from existing body text, about 6%. R3's finding is directly
relevant: what the NC compression already removed is the design table, the per-layer cost discussion
and the reporting rules — the material that makes the deliverables adoptable. Another 300 words come
out of the same place.

---

## 3. Two decisions that are yours

**D1 — the minimum effective-nucleus threshold for the expression fit.** This sets the published
coefficient and therefore every technical-noise threshold and design table in the paper.

| threshold | N | exponent (95% CI) | coefficient | note |
|---|---|---|---|---|
| none (published) | 309 | −0.507 (−0.569, −0.443) | 4.79 | includes an observation from 1 nucleus vs 9 |
| n_eff ≥ 10 | 256 | −0.434 (−0.483, −0.389) | 3.46 | excludes −0.5; retains the null −0.482 at the interval edge |
| n_eff ≥ 25 | 186 | −0.394 (−0.448, −0.341) | 2.87 | excludes both −0.5 and −0.482 |

My recommendation: adopt **n_eff ≥ 10** as primary and report the other two as a sensitivity row.
It removes the indefensible observations (a floor from one nucleus is not a floor), it is the
threshold least likely to be read as chosen to preserve the conclusion, and the matched null is
still retained, so the Layer 2 claim survives in a weakened but honest form. Do not decide finally
until B2 has recomputed the null on each support — that is what determines whether the claim
survives at all.

**D2 — venue, which is what this document is for.** See §5.

---

## 4. Sequencing

```
day 1        F3 (ALS26)  ->  A1..A12 editorial sweep         [all numbers untouched]
day 1-2      B7 figure scripts de-hardcoded                  [prerequisite for any renumbering]
day 2-3      B3 seed handling      -> new ladder interval
             B4 SI Fig 3 replot    -> strengthens the paper
             B6 heterogeneity test -> claim survives
day 3-5      B2 exponent + restricted-support null -> D1 decision -> propagate coefficient
day 5-7      B1 window recompute -> propagate reference range -> regenerate Fig 5
day 7        B5 kappa interval + per-cell-type table
             full consistency pass, rebuild docx/pdf, GB submission ready
--- NC only ---
+2-4 wk      C1 linker control, C2 packaged estimator, C3 nested decomposition
             + 300-word cut, title scope decision
```

GB-ready in about **1.5 working weeks**. NC-ready in about **6–8 weeks** on top of that.

---

## 5. Venue recommendation

**Go to Genome Biology.** The reasoning is not that NC is unreachable; it is that the three things
NC additionally requires each work against the paper.

1. **C2 (packaged estimator) contradicts the paper's own thesis.** The EIC's point is sharp: the
   NC precedents converted a caution into a portable artefact, and this paper's conclusion is that
   the artefact cannot exist numerically — κ must be measured locally. Building the estimator is
   the right answer, and it makes the *procedure* the deliverable rather than the measurement. That
   is a different paper, and a good one, but it is 1–2 weeks of tool-building plus a reframing.
2. **C1 (linker control) may not go your way.** Your linker has no cell aggregation. If Cicero or
   ArchR recovers promoter enrichment at 150 nuclei where yours does not, the Layer 3 claim
   rescopes from "the assay cannot support this at 150 nuclei" to "unaggregated correlation cannot".
   That is still publishable and still useful, but it is a materially smaller claim, and you would
   rather discover it before submission than in review.
3. **The 300-word cut removes what makes the paper useful.** R3's judgement is that NC compression
   already stripped the adoptable material; another 6% comes from the same place.

Against that, GB takes the paper as it actually is: a measurement, a procedure, a diagnostic and a
negative result, with room for each caveat to sit at the point of claim instead of in a
supplementary note. Four of five reviewers named GB as both appropriate and likely to succeed.

**The case for still trying NC**, stated fairly: the Layer 3 finding has broad implications, the
per-cell-type κ result (B5) is arguably more original than what the paper currently leads with, and
C3 would replace the confounded decomposition with a controlled one. If you build C1, C2 and C3, the
paper becomes stronger than the current draft on its own terms, and NC becomes a reasonable bet.
That is a 6–8 week detour with one genuinely uncertain outcome (C1).

**Recommended path if you want to preserve the option:** do Groups A and B first — they are required
for both venues and none of them is wasted. Then run C1 alone, because it is the cheapest of the
three and the only one that can change the science. Decide on NC after you see its result.
