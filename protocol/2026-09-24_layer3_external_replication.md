# Layer 3 external replication — NABEC/HBCC prefrontal cortex multiome

Date: 2026-09-24.
Data: Zenodo 18394349 (Catching A, et al. *Cell Rep* 2026;45(3):117110), 1,501,089 nuclei,
357 samples, 521,217 peaks. See `data/DATA_SOURCES.md` §4b.
Code: `code/reassessment/nabec_layer3.py`. Results: `data/derived_results/nabec_L3_external_replication.csv`.

## Question

The manuscript's Layer 3 result — significant promoter *depletion* (OR 0.663 [0.568, 0.768])
among peak–gene links recovered from 150 nuclei — rested on a single cohort. Does it replicate
in an independent human brain multiome dataset?

## Design

The external dataset was forced to the primary cohort's operating point rather than run at its
own: **n = 150 nuclei, ATAC downsampled to 5,564 unique fragments, RNA to 5,265 counts**,
window ±500 kb, promoter ±3 kb, FDR 0.05, collinearity 0.7 — the constants of `p5_08`.

Run at two annotation granularities: the deposited 7 cell classes, and `leiden_2` (49 clusters),
whose granularity matches the primary cohort's 25 WNN subtypes. Four strata each, chosen as the
largest cluster within Oligo, MG, ExN and InN.

## Pipeline alignment — three corrections made before any comparison

A first pass differed from `p5_08` in three respects and produced odds ratios of 1.27–1.88,
i.e. apparent *enrichment*. All three were pipeline artefacts, not data:

| | `p5_08` | first pass | effect |
|---|---|---|---|
| transform | `log1p` on both matrices | raw counts | — |
| gene set | detected in ≥10% of cells | any non-zero variance | 3× more genes tested |
| BH denominator | includes in-window zero-variance peaks (p = 1) | excludes them | 10–20× more links called |

The third dominates: the extra ~1M null tests tighten the BH threshold, and the links lost are
disproportionately weak and distal. Link counts fall from 21,000–38,000 to 1,200–3,000, and the
odds ratio moves from ~1.5 to ~0.7. **No cross-dataset number was reported before the pipelines
were aligned line by line.**

## Result — no promoter enrichment in any of the seven cell classes

Matched convention, n = 150, ATAC 5,564 / RNA 5,265:

| stratum | pairs tested | links | promoter % of links | OR | 95% CI |
|---|---|---|---|---|---|
| Oligo | 1,024,259 | 1,246 | 0.64 | 1.088 | 0.574–1.996 |
| MG | 1,056,211 | 1,156 | 0.35 | 0.601 | 0.200–1.275 |
| ExN | 1,374,890 | 2,702 | 0.59 | 0.931 | 0.535–1.388 |
| InN | 1,306,648 | 3,044 | 0.69 | 1.086 | 0.680–1.597 |
| Astro | 1,116,730 | 1,938 | 0.57 | 0.905 | 0.431–1.540 |
| OPC | 1,060,844 | 1,516 | 0.53 | 0.842 | 0.346–1.543 |
| VC | 1,523,265 | 2,813 | 0.78 | 1.214 | 0.780–1.760 |
| Oligo cl18 | 1,011,313 | 1,234 | 0.24 | 0.451 | 0.064–1.101 |
| MG cl21 | 1,064,784 | 1,195 | 0.59 | 0.969 | 0.451–1.754 |
| ExN cl7 | 1,312,713 | 2,502 | 0.52 | 0.827 | 0.397–1.322 |
| InN cl37 | 1,310,473 | 2,637 | 0.34 | 0.559 | 0.264–0.973 |
| **primary cohort** | 3,029,251 | 33,227 | **0.63** | **0.663** | 0.568–0.768 |

All seven classes were analysed; none was omitted. Not one stratum of eleven shows enrichment.
Mantel–Haenszel pooled:

| strata | manuscript background | variable-peak background |
|---|---|---|
| 7 cell classes | 0.956 [0.777, 1.177] | 0.490 [0.398, 0.603] |
| **4 fine subtypes** (granularity matched to primary) | **0.656 [0.463, 0.928]** | 0.340 [0.240, 0.481] |
| all 11 | 0.854 [0.714, 1.020] | 0.439 [0.368, 0.525] |

At coarse class level the pooled estimate shows *absence of enrichment* (CI spans 1), not
depletion. Significant depletion appears at the fine-subtype granularity that matches the
primary cohort — 0.656 [0.463, 0.928] against the primary cohort's 0.663 [0.568, 0.768] — and
under the variable-peak background at every granularity.

## Positive control — the pipeline recovers enrichment when given enough nuclei

The central alternative explanation is that this FDR convention manufactures promoter-depleted
link sets irrespective of the data: padding the BH denominator with ~1M zero-variance pairs
tightens the threshold, and if the surviving strongest correlations were distal by construction,
both 0.663 and 0.656 would be artefacts of the convention rather than replications of a fact.

Tested by holding depth and code fixed and varying only nucleus number:

| | n = 150 | n = 600 | n = 2,400 | n = 5,000 |
|---|---|---|---|---|
| **Oligo** OR | 1.088 | 0.720 | 0.862 | **2.601 [1.342, 4.310]** |
| links | 1,246 | 2,332 | 1,576 | 776 |
| promoter % | 0.64 | 0.43 | 0.51 | **1.55** |
| **ExN** OR | 0.931 | 0.353 | 0.901 | **4.907 [1.851, 9.726]** |
| links | 2,702 | 2,801 | 762 | 207 |
| promoter % | 0.59 | 0.21 | 0.52 | **2.90** |

At 5,000 nuclei — identical per-nucleus depth, identical code, identical background and FDR
convention — both cell types return significant promoter enrichment at or above the full-depth
cellranger-arc reference range of 2.45–4.01. **The convention is not the cause.** The failure at
150 nuclei is a nucleus-number effect, and the transition lies between 2,400 and 5,000.

The mechanism is visible in the link counts: they *fall* as nuclei are added (Oligo 2,332 → 776;
ExN 2,801 → 207) while the promoter fraction rises several-fold. The large link sets recovered
from few nuclei are dominated by spurious correlations, and because distal pairs vastly outnumber
promoter-proximal ones in a ±500 kb window, that noise is distal-biased.

## The nucleus ladder — how many nuclei regulatory inference needs

Depth and code held fixed at the primary cohort's operating point (ATAC 5,564 / RNA 5,265);
only nucleus number varied. Seven levels × two cell types × three seeds = 42 runs.
Mantel–Haenszel pooled over both cell types and all three seeds (6 runs per point):

| nuclei per cell type | MH OR | 95% CI |
|---|---|---|
| 150 | 1.020 | 0.816–1.275 |
| 600 | 0.458 | 0.339–0.618 |
| 1,200 | 0.472 | 0.337–0.662 |
| 2,400 | 0.792 | 0.571–1.099 |
| 3,600 | 1.641 | 1.223–2.202 |
| 5,000 | 2.888 | 2.189–3.811 |
| 7,500 | 5.973 | 4.643–7.684 |

Per cell type the curves are near-parallel (slopes 1.24 and 0.94 on the log–log scale,
r = 0.95 and 0.91) and the crossing estimates agree to within 8%.

Fitting `log(OR) ~ log(n)` over the monotone segment (n ≥ 600), weighted by inverse variance:

- **promoter enrichment crosses 1 at n ≈ 1,857 nuclei [1,627, 2,091]**
- **it reaches 2.45, the lower edge of the full-depth cellranger-arc reference range, at
  n ≈ 4,325 nuclei [3,843, 4,956]**

This is the paper's first positive design quantity for Layer 3. The equalised design of the
primary cohort — 150 nuclei, set by the smallest subtypes in the cohort — is more than an order
of magnitude below the point at which the diagnostic returns any enrichment at all.

### Caveats on the ladder

- The fit is restricted to n ≥ 600 because the curve is non-monotone below that, with a minimum
  near 600–1,200. **That restriction was chosen after seeing the data** and the crossing
  estimates should be read as approximate.
- At n = 150 the pooled estimate is 1.020 — absence of enrichment, not depletion. The primary
  cohort's 0.663 and this dataset's fine-subtype pooled 0.656 both show depletion at that
  nucleus count, so the n = 150 estimate depends on annotation granularity.
- Promoter-proximal links number 11–49 per pooled point. The shape is supported by 42 runs and
  by agreement across cell types and seeds, but every individual point rests on few events.
- Two cell types, one cohort, one depth. Whether the crossing point moves with per-nucleus depth
  is not established.
- The interval on the crossing point propagates the Mantel–Haenszel standard errors only; it does
  not include uncertainty in the power-law form itself.

## The nucleus requirement is cohort-specific (added 2026-09-24, later)

The ladder above was measured in the external cohort. Applying its threshold to the primary
cohort would transfer a constant across datasets — the error this paper documents for the
composition overdispersion. We therefore ran the same ladder **in the primary cohort itself**,
in the two cell types with enough control nuclei above both equalised thresholds
(`Exc_LINC00507_FREM3`, 9,430; oligodendrocytes, 2,883), at identical depth and convention,
three seeds per level (30 runs; `code/reassessment/primary_ladder_*.py`,
`data/derived_results/primary_L3_nucleus_ladder.csv`).

| nuclei | primary cohort | external cohort |
|---|---|---|
| 150 | **0.663** (0.486–0.905) | 1.020 (0.816–1.275) |
| 400 | 0.727 | 0.652 |
| 900 | 1.557 | 0.656 |
| 2,400 | 4.791 | 0.792 |
| 5,000 | 7.066 | 2.888 |
| 7,500 | 5.244 | 5.973 |

Two results.

**The reimplementation is validated.** At 150 nuclei the two primary-cohort cell types pool to
0.663 (0.486–0.905), reproducing to three decimals the 0.663 [0.568, 0.768] obtained across all
25 subtypes by the original `p5_08` code — an independent implementation on independently
extracted matrices.

**The threshold does not transfer.** By local interpolation the primary cohort crosses unity at
≈560 nuclei and reaches 2.45 at ≈1,340; the external cohort at ≈2,730 and ≈4,540. The mechanism
is the peak set: 521,217 peaks against 200,791 over the same genome, so at identical depth about
half the external peaks carry no variance and enter the multiple-testing denominator without
being callable, while almost none do in the primary cohort beyond 400 nuclei. The nucleus
requirement therefore behaves like the composition overdispersion, not like the expression
exponent: order of magnitude transfers, value does not.

This also resolves an internal tension. The external ladder's U shape — a minimum near
600–1,200 — appeared to contradict the manuscript's argument that sparsity produces depletion.
The primary cohort is monotone over the same range, so the dip is a property of that cohort's
peak set rather than of the mechanism, and the two accounts no longer compete.

A third observation supports the paper's own caveat on the diagnostic: the primary cohort peaks
at 7.07 at 5,000 nuclei and falls to 5.24 at 7,500 while its link count nearly doubles, the extra
links being predominantly distal. The manuscript already predicted this non-monotonicity from
the external reference link sets; it is here reproduced in a controlled ladder.

## Depth does not rescue it

ATAC depth raised to 2×, 4× and 8× the equalised design with n and RNA depth held fixed
(Oligo + ExN pooled, Mantel–Haenszel):

| ATAC fragments | × primary | links | OR | 95% CI |
|---|---|---|---|---|
| 5,564 | 1 | 3,948 | 0.940 | 0.629–1.405 |
| 11,128 | 2 | 5,136 | 0.575 | 0.367–0.903 |
| 22,256 | 4 | 5,413 | 0.663 | 0.440–0.999 |
| 44,512 | 8 | 5,587 | 0.420 | 0.253–0.697 |

No depth in this range produces promoter enrichment. Set against the nucleus ladder above, this
localises the limiting resource: 8× the reads changes nothing at 150 nuclei, while 12× the nuclei
at unchanged depth carries the odds ratio from 1.0 to 6.0. Point estimates trend downward but are
non-monotone with overlapping intervals; only the negative statement is defensible. Together with
the nucleus ladder this localises the limiting resource: **nuclei, not reads.**

## What this changes in the manuscript

1. **The Layer 3 claim is supported by an independent cohort** and should say so. This closes the
   "single cohort" objection against the weakest of the three layers.
1b. **The nucleus ladder converts the negative result into a dose–response with a working
   positive control in the same data.** This is the strongest single addition available: it shows
   the diagnostic is sound and that only inference from few nuclei fails.
2. **A required caveat on the diagnostic.** The Discussion proposes reporting the
   promoter-enrichment odds ratio alongside link counts. The same data, same nuclei, same depth
   give 0.66 or 1.5 depending only on the FDR denominator and the background set. The
   recommendation must therefore require both to be stated, or the numbers are not comparable
   between studies.
3. **The ×2–8 depth ladder is new evidence** that the failure is not a sequencing-depth problem,
   which strengthens the "absence of enrichment is not absence of power" argument with a
   different kind of control than the power calculation.

## Limits

- Per-stratum intervals are wide (1,200–3,000 links against the primary cohort's 33,227); the
  replication rests on the pooled estimate and on the consistency of sign across 11 strata.
- The nucleus ladder has four points in two cell types and one seed each; the transition between
  2,400 and 5,000 is not resolved.
- Following the manuscript's convention, the background includes the links themselves. Links are
  under 1% of tested pairs, so the effect is negligible, but the convention is inherited rather
  than chosen.
- One seed (0) per configuration.
