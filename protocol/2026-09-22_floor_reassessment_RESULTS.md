# Results: re-assessing composition claims against the technical floor

Date: 2026-09-22.
Pre-specification: `2026-09-22_floor_reassessment_prespecification.md`
(SHA-256 `b5efb71d878b11c1c531f7408d6d9430c1d545cb0554649fdca41dfd03c75e94`, frozen
2026-09-22T04:26:34Z, before any statistic was computed).
Code: `code_methods/floor_reassessment.py`. Output: `results_methods/floor_reassessment/`.

## Datasets analysed

190 contrasts (dataset × annotation level × cell type × group pair), from five datasets.
No dataset was dropped silently; `dropped.json` is empty after one recorded fix (below).

| ID | Dataset | Tissue | Contrast | Contrasts |
|---|---|---|---|---|
| D1 | ALS motor cortex multiome | brain | ALS vs HC; ALS_FTD vs HC | 28 |
| D2 | SEA-AD MTG | brain | High vs Not AD (ADNC) | 27 |
| D3 | GSE290359 BA44/BA46 | brain | Disease yes vs no | 23 |
| D4 | GSE330130 MTC | brain | ALS and variants vs control | 44 |
| D5 | GSE158055 COVID PBMC | blood | severe/critical vs control | 68 |

Recorded deviations from the frozen plan:
- **GSE330130 LSC** dropped: the deposited `obs` has no diagnosis field (prespec §9).
- **GSE330130 MTC** was initially dropped by a coding error — the control label is
  `"Non-neurological control"` and the matcher used exact rather than substring matching.
  Fixed before any outcome was inspected; this changed an implementation detail, not a
  criterion.

---

## 1. The pre-specified primary outcome, and why it was the wrong question

**Result as specified.** Among the 53 contrasts reaching conventional significance
(donor-level Welch *t*, P < 0.05), the fraction falling below the floor was:

| Floor factor κ | Below floor |
|---|---|
| **1 (multinomial bound) — pre-specified primary** | **0 / 53 = 0.0%** |
| 3.90 (SEA-AD, measured) | 13 / 53 = 24.5% |
| 4.28 (ALS, measured) | 14 / 53 = 26.4% |

⚠️ **The primary outcome as I specified it is close to tautological, and I should have seen
this before freezing.** A donor-level *t*-test uses the *observed* between-donor standard
deviation, which already contains the technical variance. Anything that survives such a test
has therefore already cleared an empirical noise bar that is, in these data, a median of 4.5
times the multinomial floor. Reporting 0% is honest but uninformative: it says that
donor-level testing is self-protecting against sampling noise, which was never in doubt.

The κ = 3.9–4.3 rows are more interesting — about a quarter of nominally significant findings
have a between-donor spread *smaller* than the technical floor measured in brain predicts.
But since the manuscript's own result is that overdispersion does not transfer, these can only
be read as "what the floor would be if this dataset behaved like the two brain cohorts", not
as a measurement.

## 2. The informative analysis (post-hoc, not pre-specified)

The floor bites where an analysis does **not** use donor-level variance. Pooling cells across
donors and testing the resulting 2×2 table treats cells as independent — precisely the
assumption the manuscript's overdispersion result refutes. This is a common practice in the
published literature.

For every contrast we added a naive pooled test (χ² on cell-type vs rest × group A vs B,
cells pooled within group).

| | Contrasts | Share |
|---|---|---|
| Significant by **cell-pooled** χ² | 166 / 190 | 87.4% |
| Significant by **donor-level** *t* | 53 / 190 | 27.9% |

Among the 166 that the pooled test calls significant:

| Fails | n / 166 | Share |
|---|---|---|
| donor-level *t*-test (P ≥ 0.05) | 113 | **68.1%** |
| below κ = 1, the multinomial bound | 14 | 8.4% |
| below κ = 3.90 | 65 | 39.2% |
| below κ = 4.28 | 69 | **41.6%** |

Restricted to the four brain datasets, where the measured κ is at least comparable in tissue:
105 pooled-significant contrasts, of which **87 (82.9%)** fail donor-level testing and
**48 (45.7%)** fall below the κ = 4.28 floor.

Per dataset (pooled-significant → fails donor-level):

| Dataset | Pooled sig. | Fails donor-level | Below κ=4.28 |
|---|---|---|---|
| D1 ALS motor cortex | 24 | 22 (91.7%) | 12 |
| D2 SEA-AD MTG | 24 | 16 (66.7%) | 10 |
| D3 GSE290359 | 20 | 20 (100%) | 12 |
| D4 GSE330130 MTC | 37 | 29 (78.4%) | 14 |
| D5 GSE158055 COVID (blood) | 61 | 26 (42.6%) | 21 |

### The clearest category

Fourteen contrasts are called significant by the pooled test yet have a donor-level mean
difference **smaller than pure multinomial sampling noise** — a bar no real dataset attains,
so these cannot be real regardless of any overdispersion assumption. A worked example:

> **SEA-AD, L5 IT, High ADNC vs Not AD.** Pooled proportions 9.11% against 8.80%.
> With 1.18 million nuclei, the pooled χ² gives P = 2.1 × 10⁻⁴. The donor-level mean
> difference is 0.00072 — **0.40 times** the multinomial floor of 0.00179 — and the
> donor-level *t*-test gives P = 0.91 (42 against 9 donors).

The pooled test is not detecting a difference in composition; it is detecting that 1.18
million is a large number.

## 3. What this does and does not say

- It is **not** a re-assessment of the original papers' conclusions. Those rest on multiple
  cell types and usually on additional evidence, and several of these datasets' authors did
  use donor-level statistics. The statement is about how much *per-cell-type* evidence
  survives when cells are not treated as independent.
- The 41.6% figure imports a brain-measured κ into datasets where it was not measured. It is
  a sensitivity statement, not a measurement, and must be labelled as such.
- D5 is blood. An AIDA analysis of 562 same-suspension replicate donors gives a calibrated
  overdispersion of 1.31× for blood composition, which is a *lower* bound because it omits
  dissociation. The κ = 3.9–4.3 rows are almost certainly too strict for D5.

## 4. Assessment: does this move the manuscript to Nature Communications?

**Partly, and less than I projected.**

What works: the post-hoc analysis is a clean, quantitative demonstration that the floor
changes how a common analytical practice should be read, across five datasets and two
tissues, with a worked example that any reader can check. "Two-thirds of composition
differences that a cell-pooled test calls significant do not survive donor-level testing, and
a sixth fall below the floor we measure" is a concrete, actionable claim of the kind NC
rewards.

What does not: the strongest number (41.6%) depends on importing a κ the manuscript itself
says does not transfer, so it can only be a sensitivity. The unconditional number — below the
multinomial bound — is 8.4%, which is real and defensible but not headline-grabbing. And the
pre-specified primary outcome returned nothing, which has to be reported.

**Revised venue assessment:**

- **Genome Biology / Genome Medicine — still first choice, now stronger.** This section adds
  a direct "so what" to the composition result without needing any imported constant: the
  8.4% unconditional figure plus the 68.1% donor-level failure rate stand on their own.
- **Nature Communications — still borderline, moved a little closer but not over.** The
  honest framing weakens the headline. To actually clear NC I think the analysis would need
  to move from deposited data to **published claims as stated in papers** — i.e. showing that
  specific reported findings fall below the floor — which is a different and more laborious
  exercise, and one that carries real risk of misrepresenting other people's work.
- **Nature Methods — unchanged, not attainable.** Still a measurement/audit paper.
- **Bioinformatics / NAR GaB — floor unchanged.**

**My recommendation:** fold section 2 into the manuscript as a short Results subsection after
the composition floor, with the worked example, and keep the framing unconditional
(multinomial bound + donor-level failure rate). Do not chase NC by hardening the imported-κ
number.
