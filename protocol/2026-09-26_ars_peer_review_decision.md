# Editorial Decision Package — simulated peer review

Manuscript: "Resolution limits of single-nucleus multiome analysis: composition, expression and
regulatory inference" (sole author). Variants reviewed: `variants/MANUSCRIPT_NatureComms.md`
(primary), `MANUSCRIPT.md`, `variants/MANUSCRIPT_GenomeBiology.md`,
`variants/SUPPLEMENTARY_INFORMATION_NatureComms.md`.

Panel: ARS `academic-paper-reviewer` v1.10.0, `full` mode, 5 reviewers, two-phase sprint contract
(paper-blind pre-commitment, then paper-visible scoring). Date: 2026-09-26.

**Contract provenance disclosure.** The shipped template `shared/contracts/reviewer/full.json` was
not present on this machine. The contract used was reconstructed from `references/quality_rubrics.md`
and `references/editorial_decision_standards.md`. FC1's severity (4) and quantifier
(`at_least_1_of_5`) are therefore an orchestrator drafting choice, not a shipped baseline. This
matters: it is the sole reason the mechanical decision is Reject rather than Major Revision. See
§4.

---

## 1. Dimension scores

| Dimension (weight) | EIC | R1 Method | R2 Domain | R3 Perspective | R4 DA |
|---|---|---|---|---|---|
| D1 Originality (20%) | 82 | 76 | 68 | 78 | 70 |
| D2 Methodological Rigor (25%) | 84 | **46 block** | 70 | 68 | 60 |
| D3 Evidence Sufficiency (25%) | 76 | 52 | 68 | 76 | — |
| D4 Argument Coherence (15%) | 71 | 58 | 70 | 70 | — |
| D5 Writing Quality (15%) | 78 | 62 | 78 | 79 | — |
| D6 Literature Integration (0%) | 79 | 66 | **48 block** | 62 | — |
| D7 Significance & Impact (0%) | 84 | 64 | 66 | 72 | — |
| **Weighted** | **79** | **57.7** | **70.3** | **74** | **62.1** |

Panel mean 68.6. Spread 21.3 points, and it is not noise: the two highest scores come from the
reviewers who did not audit the statistics, the two lowest from those who did. The manuscript reads
better than it audits.

## 2. Cross-reviewer failure-condition matrix

| FC | Sev | Quantifier | EIC | R1 | R2 | R3 | R4 | Met |
|---|---|---|---|---|---|---|---|---|
| FC1 central number unverifiable/inconsistent | 4 | ≥1 of 5 | no | **yes** | declined | no | **yes** | **YES** |
| FC2 headline claim unsupported by design | 3 | ≥2 of 5 | no | no | no | no | no | no |
| FC3 negative result plausibly a pipeline artefact | 3 | ≥2 of 5 | no | **yes** | **yes** | **yes** | **yes** | **YES** |
| FC4 overstated generalisation | 2 | ≥2 of 5 | **yes** | **yes** | **yes** | **yes** | **yes** | **YES** |
| FC5 venue scope/format mismatch | 2 | ≥2 of 5 | no | n/a | n/a | n/a | n/a | no |
| FC6 reporting gaps | 1 | ≥1 of 5 | **yes** | **yes** | **yes** | **yes** | **yes** | **YES** |

R2 declined FC1 explicitly as miscalibration, routing its two items (the 6.591 literal, Table 1's
362 vs 357) to FC6 on the ground that neither is *central*. That judgement is recorded, not
overridden.

## 3. Findings independently re-verified by the orchestrator

Every item below was checked against the code and deposited files, not taken from a reviewer report.

**V1 — Window mismatch in the headline comparison (R4 C2). CONFIRMED.**
`code/analysis/p5_08_claimA_features.py:42` sets `WINDOW = 500_000`, so the 0.663 is a ±500 kb odds
ratio. All four deposited reference sets carry `window = 1000000`
(`data/derived_results/p5_promoterOR_*.csv`: 4.012, 3.575, 3.074, 2.451). The Introduction's
"ranges from 2.45 to 4.01, and in the equalised analysis here it is 0.663", Fig. 5a's reference
shading, and the Discussion therefore compare a ±500 kb value against four ±1 Mb values. Methods
line 515 states "cross-study comparisons were recomputed at a common ±1 Mb window"; no ±500 kb
reference value exists anywhere in the deposited data. The manuscript's own Supplementary Note 3,
error six, names this exact failure and reports it corrected. Direction: the paper's own
4.01 (±1 Mb) vs 6.59 (±1.7 Mb) shows OR rises with window, so the reference range is inflated
relative to ±500 kb and the reported gap is overstated by an undetermined amount.
The nucleus ladder is unaffected — it is internally window-consistent.

**V2 — Expression exponent is support-dependent and the sensitivity is unreported (R1). CONFIRMED,
with a correction to R1.**
`code/analysis/p5_09_floor_scaling.py:90` filters only `median_abs_log2FC > 0 & n_eff > 0`: there is
no minimum nucleus count. 18 of 309 observations have n_eff < 5; the smallest is HC27 Microglia at
n1 = 1, n2 = 9. Refitting with donor-clustered bootstrap:

| support | N | exponent (95% CI) | coefficient | excludes −0.5 | excludes −0.482 |
|---|---|---|---|---|---|
| all (published) | 309 | −0.507 (−0.569, −0.443) | 4.79 | no | no |
| n_eff ≥ 10 | 256 | −0.434 (−0.483, −0.389) | 3.46 | yes | **no** |
| n_eff ≥ 25 | 186 | −0.394 (−0.448, −0.341) | 2.87 | yes | yes |

The published fit reproduces exactly. **Correction to R1:** R1 reported that both restricted fits
exclude the matched null −0.482. At n_eff ≥ 10 they do not — −0.482 sits inside (−0.483, −0.389).
Only n_eff ≥ 25 excludes both.
**Caveat R1 did not state, and it matters.** The −0.482 null was computed on the same full support,
so the published observed-vs-null comparison is internally matched and defensible *on that support*.
Establishing that the observed exponent departs from its null under restriction requires recomputing
the null on each restricted support — and the per-pair motor-cortex null is not deposited (only
`seaad_L2_null.csv` and `Supplementary_Table_3c` are, both Seattle), so that check cannot currently
be run from the deposited artefacts. What is established: the exponent and coefficient are
support-dependent, no sensitivity is reported, and the coefficient 4.79 propagates into every
|log₂FC| threshold in the Discussion and into the design tables.

**V3 — Supplementary Fig. 3 mislabels a depth contrast as a convention contrast (R4 C1, R3, R1).
CONFIRMED.** In `data/derived_results/nabec_L3_convention_sensitivity.csv` the `convention` column
has two levels whose ATAC depths differ 5–22×: `unmatched` at depth_atac 27,973–120,798 (OR 0.829 to
1.631) and `p5_08 matched` at 5,564 (OR 0.601 to 1.214). The Discussion's "holding nuclei, depth and
code otherwise fixed, those two choices alone moved it between 0.66 and 1.5" is false for this
figure: depth is not fixed. The genuine convention contrast is within-row at fixed depth, `OR` vs
`OR_varpeaks`, and it moves the ratio **downward** (0.601→0.253, 1.086→0.649), never across unity.
Correcting this *strengthens* the paper's conclusion. Separately, the high-depth `unmatched` rows
reach OR 1.631, which bears on "nucleus number rather than depth was the limiting resource".

**V4 — "309 replicate pairs" is a units error, and the 26/27 statement is wrong for Layer 2.
CONFIRMED.** `Supplementary_Table_3a` has exactly 309 rows over 27 donors × 12 cell types (324
minus 15 absent combinations). The independent units are 27 donor-chip pairs, not 309. The 27th
donor is ALS26, which is present in the expression table but absent from
`Supplementary_Table_2a` (26 donors). The manuscript's "these 26 pairs underlie every replicate
analysis reported here" is therefore false for the expression analysis.

**V5 — Figure 5 is never cited in the main text. CONFIRMED.** In all three variants the only
`Fig. 5` match in Results+Discussion is "Supplementary Fig. 5". Main-text Figure 5 has zero
call-outs.

**V6 — Table 1 reports "362 samples" for NABEC/HBCC; Methods and three other places report 357.
CONFIRMED** in all three variants.

## 4. Editorial decision

**Mechanical outcome under the contract: REJECT.** FC1 (severity 4, `at_least_1_of_5`, action
`reject`) fired for two reviewers on three independently confirmed instances (V1, V2, V3). Severity
4 is highest; FC1 governs. The contract's `override_ladder` contains only `scoring_plan_dissent`,
which does not reach this case.

**Substantive panel consensus: MAJOR REVISION.** R1 and R4 — the two reviewers who fired FC1 — both
volunteered, unprompted and independently, that FC1 as drafted is too blunt for these defects, and
that Major Revision is the substantively correct outcome. Their reasoning is sound and the
orchestrator adopts it as the reportable recommendation: every FC1 instance is a labelling or
reporting failure over an analysis that reproduces, repairable with code and public data already in
hand, on a timescale of days. No reviewer fired FC2 — no one holds that the study design fails to
support its headline claim.

The two labels diverge because of FC1's drafting, which was the orchestrator's own reconstruction
(§ contract provenance). Weight the substantive consensus, not the mechanical label.

**The decision that matters for the author: do not submit to any venue in the current state.** V1
alone would be found by a competent referee, and being caught mis-comparing windows on the exact
error your own Supplementary Note says you corrected costs more credibility than the finding is
worth.

## 5. Venue verdict — 5 of 5 reviewers

| Reviewer | Verdict |
|---|---|
| EIC | Genome Biology recommended, expected to succeed. NC not as submitted; reachable only if the floor-estimation procedure is packaged as a runnable estimator and the title is either broadened beyond brain or narrowed to it. Not a Matters Arising. Nature Methods Analysis live only if the estimator is built. |
| R1 | The measurement is NC-grade; the current statistical reporting is not. Two of three blocking defects fixable without new data. |
| R2 | Genome Biology after the linker control and literature repair. NC only if the Layer 3 claim is rescoped or defended against an aggregation-based linker. |
| R3 | Genome Biology — the NC compression strips the design table, per-layer cost discussion and reporting rules that make the deliverables adoptable. |
| R4 | Not publishable as is; after reanalysis, Genome Biology over NC. |

**Unanimous: Genome Biology. Zero support for Nature Communications as submitted.**

The EIC's specific rebuttal of the author's own NC precedent argument (Squair et al., scPower,
scCODA) is the load-bearing reasoning: what those papers share is not a caution but a caution
converted into a *portable artefact*, and this manuscript's headline conclusion is that the artefact
cannot exist numerically ("κ does not transfer... must be measured in the dataset at hand"). On the
EIC's reading the precedent set argues *against* an NC Article, and the route to NC is to make the
measurement procedure itself the deliverable.

## 6. Revision roadmap, prioritised

### Binding before any submission
1. **V1.** Recompute all four reference link sets at ±500 kb and restate the range, or drop the
   numerical comparison and rest the argument on the nucleus ladder. Delete or correct the Methods
   sentence at line 515.
2. **V2.** Report the support sensitivity of the exponent and coefficient. Recompute the matched
   multinomial null on each restricted support and deposit the per-pair motor-cortex null. If the
   exponent proves support-dependent against a matched null, restate the Layer 2 claim.
3. **V3.** Replot Supplementary Fig. 3 from `OR` vs `OR_varpeaks` at fixed depth; restate the
   Discussion sentence as the three documented pipeline differences plus depth, citing the protocol
   record. Say that the corrected contrast strengthens the conclusion.
4. **V4, V5, V6.** Fix the 309 units error and the 26-vs-27 donor statement; cite Figure 5 in the
   main text or demote it; reconcile 362 vs 357.
5. **Ladder significance (R1).** The 150-nucleus primary estimate pools three seeds as
   Mantel–Haenszel strata; seed 0 alone gives 0.919 (0.570–1.482) against the published 0.663
   (0.486–0.905). Re-samplings of one dataset are not independent strata. Re-derive the interval.
6. **Abstract (R4 M1).** "reproduced in an independent cohort" attaches to depletion; the external
   pooled estimate at 150 nuclei is 1.020 (0.816–1.275) and only matched-granularity strata give
   0.656. Use "the absence of promoter enrichment reproduced".

### Strongly recommended
7. **κ heterogeneity (R4 M2).** Per cell type κ spans 1.32–7.44 in `p5_techrep_overdispersion.csv`,
   rising with abundance; the 126,670-nucleus design figure becomes 38,461 at the median type and
   383,037 at oligodendrocytes. Report per-type κ; give the design cost as a range; and note κ has
   no confidence interval anywhere (R1 obtains 4.27 [3.03, 5.58]).
8. **κ contradiction (R4 M3).** "κ does not transfer between cohorts" (Discussion) against
   "comparable, 4.28 against 3.90" (Results/Abstract). One sentence fixes it.
9. **Linker control (R2, R4 M7).** The Layer 3 controls all hold the linker fixed; the linker is
   per-nucleus Pearson on log1p counts with no aggregation, which is what Cicero and ArchR aggregate
   away. Run the pipeline against one published `feature_linkage.bedpe` at matched cell number and
   window. R2 rates this the single condition for an NC-grade Layer 3 claim.
10. **Literature (R2).** No peak-gene linking method other than cellranger-arc is cited. Absent:
    Cicero, ArchR, Signac/LinkPeaks, SCENT, SCARlink, Enhlink — the last two named in the author's
    own novelty audit — plus CRISPR enhancer-gene ground truth.
11. **Ladder strata (R4 M9).** The top two primary rungs rest on one cell type; oligodendrocytes
    stop at 2,400. Correct Table 3's legend and the Methods.
12. **PsychAD modality (R2).** PsychAD is hashed snRNA-seq, stated nowhere; 1,209 of 1,235
    composition pairs are snRNA-seq under a "multiome" title.

### Worth doing, not blocking
13. Metrology framing (R3): zero hits for repeatability / LOD / Bland-Altman / Gage R&R / MIQE
    across all variants. An opportunity, not a defect.
14. R3's cost asymmetry: the expression floor needs one replicate pair; the composition *exponent*
    needed ~10² pairs by the paper's own evidence. State a provisional κ with scope and a revision
    mechanism rather than disowning the two measured values.
15. Reproducibility gaps (FC6): hardcoded values in `fig2_fig3_ggplot2.R:120`, `fig6_ggplot2.R:79`,
    `fig1_fig4_fig5_ggplot2.R:257` falsify the cover letter's "every figure regenerates";
    `nabec_L3_convention_sensitivity.csv` has no generating script; a Chinese column header in a
    deposited table; GB abstract over length.

## 7. Credit where the panel agreed

Three reviewers independently reproduced nearly every headline number from the deposited artefacts,
several to three decimal places: overdispersions 4.28/3.90/1.36/1.22, the expression fit, all four
crossing points (562/2,733/1,337/4,544), the 6,915/126,670 design figures, all four reference ORs,
and every number in the 190-contrast reassessment. R4 recorded that its strongest briefed attack —
that 0.663 measures a defect in the author's own linker — **failed on the evidence**: the ladder's
rise holds under the alternative denominator convention in both cohorts (0.456→7.118), the depth
series excludes depth across an eightfold range, and 0.663 reproduces across two independent
implementations. R4 also found the manuscript's self-reported errors, frozen pre-specification and
recovery of the source study's finding to be load-bearing rather than decorative.

Two reviewers noted that a correct reanalysis would make the paper *stronger* than it currently
claims: the true convention effect pushes the ratio further below unity (V3), and R2's
nucleus-and-abundance-matched decomposition widens the threefold composition reduction to 4.4×.
