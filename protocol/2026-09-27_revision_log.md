# Revision log — response to the 2026-09-26 panel

Input roadmap: `protocol/2026-09-26_ars_peer_review_decision.md` §6.
Plan: `protocol/2026-09-26_remediation_plan.md`.
All changes applied to `MANUSCRIPT.md`, `variants/MANUSCRIPT_NatureComms.md`,
`variants/MANUSCRIPT_GenomeBiology.md` and the two SUBMISSION_* folders unless noted.

**No Response-to-Reviewers letter exists.** This log is the traceability input. Every claim
below states where to verify it.

## Binding items (roadmap §6.1)

| # | Roadmap item | What was done | Verify at |
|---|---|---|---|
| 1 | V1 window mismatch | Downloaded the four 10x ARC 2.0.0 analysis bundles (1.31 GB); added `--window` to `p5_15_promoter_OR_external.py` **and** made it restrict the detected link set as well as the background; recomputed all four reference sets at ±500 kb (matched to `p5_08`'s `WINDOW = 500_000`) and at ±1 Mb. Reference range **2.45–4.01 → 2.26–3.56**. Methods sentence corrected. | `code/analysis/p5_15_promoter_OR_external.py`, `code/analysis/p5_19_recompute_refs_at_window.sh`, `data/derived_results/window_500000/`, `window_1000000/`, `supplementary_tables/Supplementary_Table_5_promoterOR_*.csv` |
| 2 | V2 exponent support-dependence | Reported the support sensitivity in Results; established that the published null (−0.482) and the 1.30 excess are **SEA-AD's**, computed pair-for-pair, and said so; reported the observed-vs-null difference on identical support. Threshold **not** changed — see "Deferred" below. | `p5_17_B2_exponent_support_sensitivity.csv`; Results, Layer 2 section |
| 3 | V3 SuppFig 3 mislabelling | Replotted from `OR` vs `OR_varpeaks` at fixed depth 5,564. Discussion restated: the convention moves the ratio **0.60–1.21 → 0.25–0.65**, further below one, never across it. The "0.66 and 1.5" band is gone. | `code/figures/supp_figs_ggplot2.R` SF3 block; Discussion "Three properties" paragraph |
| 4 | V4/V5/V6 | "309 replicate pairs" → "309 donor by cell-type observations" (6 places); deleted the false "these 26 pairs underlie every replicate analysis"; Figure 5 now cited (5a/5b in NC, 5a/5b/5c in the long versions); Table 1 "362 samples" → "357 samples". | grep the three manuscripts |
| 5 | Ladder seed handling | Seeds are re-draws, not independent strata. Recomputed by pooling within seed and combining across seeds: 150 nuclei **0.700 (0.319–1.536)**, 2,400 nuclei **4.89 (3.61–6.63)**. Stated in the Table 3 legend. The 150-nucleus point loses significance; the positive control holds. | `p5_17_B3_ladder_seed_handling.csv`; Table 3 legend |
| 6 | Abstract "reproduced" | → "and no enrichment appeared in an independent cohort either". | NC + full abstracts |

## Strongly recommended (roadmap §6.2)

| # | Item | What was done | Verify at |
|---|---|---|---|
| 7 | κ heterogeneity | κ now carries a bootstrap interval (4.28, 2.56–5.72); reported that κ runs 1.53–8.50 across cell types and rises with abundance (r = 0.85). Design cost range deposited: 6,915 / 44,435 / 126,593 / 499,575. | `p5_17_B5_kappa_by_celltype.csv`; Results, composition paragraph |
| 8 | κ contradiction | "κ does not transfer" → "agreed to within 10% in the two cortex cohorts measured here, which is not enough to treat it as universal". | Discussion |
| 9 | Linker control (C1) | **Not done.** Requires an aggregation-based linker run. | — |
| 10 | Literature repair | **Not done.** | — |
| 11 | Ladder strata | Table 3 legend now states that the 5,000 and 7,500 rows rest on `Exc_LINC00507_FREM3` alone (oligodendrocytes provide 2,883 control nuclei in total). | Table 3 legend |
| 12 | PsychAD modality | Stated: hashed single-nucleus RNA-seq, bearing on the composition and expression layers only. | Results, replicate-hierarchy section |

## Also done, not on the roadmap

- B6: the "incompatible" composition exponents now carry a donor-clustered bootstrap test, **0.432 (0.140–0.791), P = 0.003**.
- B7: the three hardcoded figure literals removed. All four cohort κ values are recomputed from the deposited per-pair tables by `p5_17` and read from CSV by `fig2_fig3` and `fig6`; Fig 6b's coefficient label is derived from the fit. Chinese column header `过度离散倍数` → `overdispersion_ratio`.
- Fig 5b's previously unsourced 6.591 was recomputed as **6.593** at ±1.7 Mb and is now reproducible; the panel shows three windows (3.56 / 4.05 / 6.59).
- A9–A12: non sequitur split; the 4.79 collision disambiguated; replicate-scarcity claim narrowed to cross-batch multiome replicates.
- Seven earlier editorial items (A1–A8) as listed in the plan.

## Length

NC body cut from 5,322 to **4,999 / 5,000** words including Methods; abstract 149/150.
Cuts were restatement only. Removed: the Introduction's results-preview paragraph (122 w, which
also still asserted the "significant promoter depletion" that item 5 withdrew); the Discussion's
one-sentence paragraph restating its own opening; the `|log₂FC|` endpoint values (readable from
Fig. 6b, with the 100-nucleus operating point retained); the TSS-enrichment aside; the crossing-point
property sentence duplicated from the Discussion; Supplementary Note 1's paragraph compressed to a
pointer. **No evidence, caveat or limitation was removed.**

## Deferred, with reasons

- **D1 threshold.** The author chose n_eff ≥ 10. On testing, that threshold **inverts the Layer 2
  claim**: exponents go from consistent (difference −0.002, P = 0.98) to different (+0.060,
  P = 0.033), while coefficients converge from 31% apart (4.79 vs 3.65) to 0.8% apart
  (3.464 vs 3.438). Applying it silently would have rewritten the paper's Layer 2 narrative from
  "the exponent transfers" to "the coefficient transfers". Held for an explicit decision. The
  practical thresholds barely move (0.46 → 0.47 at 100 nuclei) and the excess over null is stable
  (1.30 → 1.29).
- **C1 external linker control, C2 packaged estimator, C3 within-cohort nested decomposition.**
  Not attempted; these were scoped as Nature Communications / Nature Methods requirements.
- **ALS26.** Its second chip carries 37 nuclei in total, so it fails the stated ≥50-nuclei rule, yet
  contributes 9 of the 309 expression observations. Coupled to D1 and deferred with it.
