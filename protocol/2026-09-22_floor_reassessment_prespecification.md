# Pre-specification: re-assessing published composition claims against the technical floor

Date: 2026-09-22. Frozen before any statistic was computed.
Purpose: supply the increment identified in `2026-09-22_novelty_audit_resolution_limits.md`
as the single most effective way to move the manuscript toward *Nature Communications* —
a demonstration that the measured floor changes how existing composition claims should be read.

---

## 0. The binding constraint

The manuscript's own composition result is that **overdispersion does not transfer between
datasets** (4.28-fold against 3.90-fold in magnitude, and opposite scaling exponents).
It is therefore **not permitted** to apply either cohort's overdispersion factor to a third
dataset as if it were a constant. Doing so would contradict the paper's central claim.

What *is* transferable is the **multinomial expectation**, which is a mathematical lower
bound on sampling noise that no real dataset attains. All primary results are stated against
that bound. The measured brain factors enter only as a labelled sensitivity range.

## 1. Question

Among cell-type composition differences that would be reported as findings by a conventional
donor-level test, what fraction cannot be distinguished from technical variation?

## 2. Datasets (fixed; no additions after freeze)

| # | Dataset | Tissue | Grouping variable | Floor measured in this manuscript? |
|---|---|---|---|---|
| D1 | ALS motor cortex multiome (primary cohort) | brain | `Case`: ALS / ALS_FTD / HC | yes (4.28×) |
| D2 | SEA-AD MTG snRNA-seq | brain | AD neuropathology group | yes (3.90×) |
| D3 | GSE290359 (BA44/BA46) | brain | `Disease?` | no |
| D4 | GSE330130 (MTC, LSC) | brain | `Diagnosis` | no |
| D5 | GSE158055 COVID PBMC | blood | CoVID-19 severity | no |

A dataset is dropped only if, on inspection, it lacks a usable donor identifier, cell-type
label or grouping variable. Any drop is recorded with its reason.

## 3. Unit of analysis

One **contrast** = (dataset × annotation level × cell type × pair of groups).
Annotation levels: the coarse and the fine level native to each dataset, analysed separately
and reported separately. No level is chosen after seeing results.

## 4. Inclusion rules (fixed)

- Donor contributes ≥ 50 cells in total (matching the manuscript's replicate threshold).
- Each group has ≥ 5 donors.
- The cell type is present (pooled proportion > 0) in the dataset.
- Cell types are **not** filtered on abundance; rare types are retained and reported.

## 5. Statistics

For a cell type with pooled proportion *p*, donor *d* contributing *n_d* cells, and groups
A and B with *D_A*, *D_B* donors:

**Observed effect** Δ = | mean(f_d | A) − mean(f_d | B) |, where f_d is donor *d*'s proportion.

**Technical variance of Δ** under floor factor κ:

  Var_κ(Δ) = κ² · p(1−p) · [ Σ_{d∈A}(1/n_d)/D_A² + Σ_{d∈B}(1/n_d)/D_B² ]

**Floor** at factor κ: F_κ = 1.96 · √(Var_κ(Δ)).
A contrast is *above the floor* if Δ > F_κ.

κ values:
- **κ = 1** — primary. The multinomial bound; unattainable in practice, so this is the most
  permissive possible test and any claim failing it fails unconditionally.
- **κ = 3.90 and κ = 4.28** — sensitivity, labelled as measured in brain in this manuscript.
  For D3–D5 these are explicitly *imported* values, not measurements, and are reported as
  "what the floor would be if this dataset behaved like the two brain cohorts".

**Conventional test** (what determines whether a difference would be reported): two-sided
Welch t-test on donor-level proportions; nominal P < 0.05. Benjamini–Hochberg across all
cell types within a dataset × level × contrast is reported alongside but is **not** the
primary filter, because published claims are frequently made on nominal significance.

## 6. Primary outcome

Among contrasts with conventional P < 0.05, the proportion with Δ ≤ F_κ, reported for each
κ, overall and per dataset.

## 7. Secondary outcomes

- The same proportion after BH correction.
- Δ / F₁ (the effect expressed in multiples of the multinomial floor), distribution.
- Whether contrasts failing the floor are concentrated in rare cell types.
- The minimum detectable Δ at each dataset's median cell number, as a design table.

## 8. Stopping rule and honesty commitments

- Every contrast passing inclusion is reported. No cell type, dataset or level is dropped
  after seeing its result.
- If the primary outcome is small (few claims fall below the floor), that is the result and
  it is reported as such — it would mean published composition claims are generally robust
  to sampling noise, which is a useful and publishable negative.
- The analysis is **not** a re-assessment of the original papers' conclusions, which rest on
  more than one cell type and often on additional evidence. It is a statement about how much
  of the *per-cell-type* evidence survives a sampling-noise floor. Wording must reflect this.
- ⚠️ D5 is blood; the manuscript's floor was measured in brain. The κ = 3.9–4.3 sensitivity
  is reported for D5 only to show the range, and the text must state that blood's own
  overdispersion is unmeasured here. Context available but not to be presented as measured
  for this dataset: an AIDA analysis of 562 same-suspension replicate donors gave a calibrated
  overdispersion of 1.31×, which is a lower bound for blood because it omits dissociation.

## 9. What would make this fail

- If donor identifiers cannot be reconstructed for a dataset, it is dropped (recorded).
- If a grouping variable is ordinal with no natural dichotomy, the extreme groups are
  contrasted and that choice is recorded here in advance.
- No result is reported from a dataset whose inclusion rules were relaxed to obtain it.
