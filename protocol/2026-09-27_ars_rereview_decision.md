# Re-review decision package (ARS `re-review`, 2026-09-27)

Round 1: `protocol/2026-09-26_ars_peer_review_decision.md` (5 reviewers, Minor/Major split, unanimous GB).
Revision log: `protocol/2026-09-27_revision_log.md`.
Panel this round: EIC (per protocol) + methodology verifier (deliberate addition — FC1 fired on
numerical grounds in round 1 and the revision changed those numbers; an EIC-only re-review could not
verify them). Contract unchanged, including its disclosed FC1 drafting caveat.

## 1. Where the two reviewers agree

Both recomputed independently and both report **the numerical work is sound**. Between them they
re-derived ~50 quantities from raw or deposited data — the four ±500 kb ORs to four decimals, the
seed-recombined ladder, κ and its interval, the B2 support table, B6, the word counts — and **found no
number wrong**. The methodology reviewer went to the raw 10x bundles and independently confirmed the
inconsistent-window diagnosis, naming the mechanism: Cell Ranger ARC's link reach is ≤1 Mb (max bedpe
distance 999,997) while the script tests |peak midpoint − GENCODE TSS|, which reaches 2.18 Mb, so the
old default bounded the denominator and not the numerator.

Both also agree on the failure mode: **the reporting layer lags the analysis layer.** Three odds-ratio
conventions, two κ estimators and two seed conventions now coexist across the manuscripts, the
Supplementary Information, the figure legends and the data documentation.

## 2. Where they disagree, and the arbitration

**FC1.** EIC: fires. Methodology: does not fire. Both are correct on their own axis — the EIC judges
internal consistency ("contradicted by the manuscript's own supplementary material", which is FC1's
literal wording), the methodology reviewer judges numerical correctness. Under the contract's
`at_least_1_of_5` quantifier FC1 fires, and mechanically that is `reject`. **Disregard the mechanical
label**, for the same reason recorded in round 1: FC1's severity-4 drafting is the orchestrator's own
reconstruction, not a shipped baseline, and no numerical claim is wrong.

**FC4** fires for both reviewers, meeting `at_least_2_of_5` (severity 2 → major revision).
**FC6** fires for both.

**Substantive decision: MAJOR REVISION**, one step worse than round 1's Minor Revision consensus —
driven entirely by defects introduced or left by the revision itself, not by anything in the original
science. Score 79 → 75.6 (EIC); the only regressed dimension is argument coherence (71 → 64).

## 3. The D1 resolution neither option offered

The author was choosing between keeping the published full support and adopting an n_eff ≥ 10
threshold. The methodology reviewer identified a third option that is not a new convention at all:
**apply the inclusion rule the manuscript already states.** ALS26 contributes 9 of the 309 expression
observations despite carrying 37 nuclei on its second chip, failing the stated "at least 50 nuclei on
both" criterion.

| motor-cortex rule | N | donors | exponent | coefficient | cross-cohort *P* |
|---|---|---|---|---|---|
| as published | 309 | 27 | −0.507 | 4.79 | 0.98 consistent |
| **the paper's own ≥50-nuclei rule (drop ALS26)** | **300** | **26** | **−0.497** | **4.58** | **0.79 consistent** |
| n_eff ≥ 10 | 256 | 27 | −0.434 | 3.46 | 0.012–0.033 **inconsistent** |

Applying the paper's own rule removes the rule violation, removes the 26-vs-27 contradiction, and
**preserves the Layer 2 claim**. The n_eff ≥ 10 threshold is a post-hoc convention that breaks it.
Recommendation: adopt the stated rule as primary and report n_eff ≥ 10 and ≥ 25 as a sensitivity row.

## 4. Binding items for the next pass

### Estimator reporting (methodology reviewer)
1. **Table 3's body and its legend are different estimators.** Body = Mantel–Haenszel without
   Haldane–Anscombe, background not excluding detected links. The legend's 0.700 uses HA + exclusion.
   So the 0.663 → 0.700 move is **the HA correction, not the seed handling**; on the body's own
   convention seed recombination gives **0.649 (0.284–1.484)**. Rewrite the legend to contrast
   like with like, and state one convention for the whole paper.
2. **Methods misdescribes the estimator behind every odds ratio in the paper** — it specifies HA plus
   2,000 parametric bootstrap replicates, whereas the headline, all of Table 3, the external pooling
   (0.656 / 0.956) and the depth series are computed without HA, and the headline's 25-subtype cluster
   bootstrap is not described anywhere.
3. **Depth series still pools seeds as strata.** Corrected: 2× 0.570 (0.349–0.930),
   4× **0.589 (0.327–1.061), crossing one**, 8× 0.503 (0.265–0.956). The claim "significantly below one
   throughout" fails and must be restated.
4. **Regression introduced by B7**: the figures now read the pair-level κ (4.2678 → prints 4.27) while
   the text and Table 2 use the RMS over 12 cell types (4.28). Two κ estimators in one paper. Also the
   two intervals quoted in one sentence resample different units (12 cell types vs 80 donors).
5. The ±1.7 Mb bar in Fig. 5b reintroduces the numerator/denominator mismatch just fixed; recompute it
   under the corrected rule or drop it.
6. 0.663 / 0.568 / 0.768 remain hardcoded in three scripts with no deposited 2×2.

### Propagation (EIC; independently confirmed)
7. **Supplementary Figure 3's legend still asserts the withdrawn 0.66–1.5 band and still says "at the
   same depth"** — in all five Supplementary Information files. This is the exact claim the revision
   identified as false. The revision log's "the band is gone" is **false as written**; the fix reached
   the plotting code and the Discussion only.
8. Figure 5b's legend gives two windows; the figure and the text give three. All three variants.
9. `data/DATA_SOURCES.md` still carries the ±1 Mb table and states the quoted range "was always"
   2.45–4.01. The deposit now holds two disagreeing ±1 Mb reference-OR sets.
10. "overdispersion does not transfer" survives in `MANUSCRIPT.md:384` and
    `MANUSCRIPT_GenomeBiology.md:413`, contradicting the corrected wording elsewhere.
11. The 26-vs-27 arithmetic contradiction sits inside single documents (resolved by item D1 above).
12. **The Genome Biology variant is the dirtiest of the three** and carries the most residual
    contradictions — which matters, because GB is the recommended venue.

## 5. Venue: unchanged

**Genome Biology**, with a higher expectation of success than round 1, conditional on the bookkeeping
pass above being applied *to the GB variant in particular*.

**Nature Communications**: unchanged. C2 (the packaged estimator) was explicitly declined, so the gate
stands exactly where it stood.

**Nature Methods (Analysis)**: the EIC judges this the stronger Nature-tier route, on a sharper
argument than was available before — the format rewards **a portable diagnostic rather than a portable
constant**, which dissolves the tension between this paper's central finding (the constants do not
transfer) and the NC precedent set (which rewards portable artefacts). It is gated by the same
estimator, plus a second linker (a harder requirement at NM than at NC), plus the literature repair,
which is triage-disqualifying and is the cheapest item on the entire roadmap.

## 6. Credit

The EIC's summary, which the methodology reviewer's recomputation supports: the revision recomputed
from 1.31 GB of re-downloaded source, found a second-order bug inside its own fix (the window bounded
the background but not the detected set), and accepted a **weaker** headline gap as a result; disclosed
the SEA-AD null attribution against interest; corrected the convention analysis in the direction that
helps the paper without overclaiming; and added three unprompted improvements. The 323-word cut was
verified by bidirectional sentence diff against the pre-revision PDF as substantially restatement, with
two imprecisions in the claim that no caveat was removed — the |log₂FC| series (recoverable from the
retained power law) and one clause that survives in the cited Supplementary Note 1. Nothing was
smuggled in. The measurement is sound; the bookkeeping is not.
