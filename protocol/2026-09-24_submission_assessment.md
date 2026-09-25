# Submission assessment — novelty, venue, exposure

Date: 2026-09-24, after the Layer 3 external replication, the nucleus ladder, the design
tables and the prior-art attribution were written in. Supersedes the venue section of
`2026-09-22_floor_reassessment_RESULTS.md` §4.

## 1. What is new

**The core claim, and the one nothing in the literature contests.** Every method cited in
the composition literature — scDC [12], scCODA [13], propeller [14], sccomp [15] — establishes
that between-subject variance exists and must be modelled, and Squair et al. [16] quantified
the same error in its expression form. None of them measures **how large the technical
component is**, because none had replicate libraries from the same donor. This paper measures
it, for three inference layers, in cohorts where the replicate provenance is documented.

Four results follow that could not have been obtained otherwise:

1. **Transferability is layer-specific.** The expression floor's scaling reproduces across
   cohorts (−0.507 and −0.505), so one calibration point fixes the curve. The composition
   floor's magnitude reproduces (4.28, 3.90) but its scaling does not, so whether more nuclei
   help must be established within each dataset.
2. **The excess decomposes by replicate provenance.** Against 1,129 PsychAD donors whose two
   libraries are split aliquots of one suspension, removing tissue sampling and dissociation
   cuts the composition excess about threefold and leaves the expression excess unchanged. The
   two layers fail for different reasons and take different remedies.
3. **Regulatory inference at an equalised 150-nucleus design returns promoter depletion**,
   reproduced in an independent cohort of 1,501,089 nuclei forced to the same operating point
   (0.656 [0.463, 0.928] against 0.663 [0.568, 0.768]) across eleven strata.
4. **The failure is a nucleus-count effect with a measured threshold.** Holding depth and code
   fixed, enrichment crosses one at ≈560 nuclei per cell type in the primary cohort and ≈2,730
   in the external one, while an eightfold increase in reads at 150 nuclei never produces
   enrichment. ⚠️ Revised 2026-09-24 (later): the earlier figures ≈1,860 / ≈4,330 came from the
   external cohort alone. Running the same ladder in the primary cohort showed the threshold is
   cohort-specific by four- to fivefold, tracking peak-set granularity, so only its order of
   magnitude transfers — which makes it parallel to the composition constant rather than to the
   expression exponent.

**Two secondary contributions.** The promoter-enrichment odds ratio as a zero-cost diagnostic,
with its convention dependence characterised (the same nuclei give 0.66 or 1.5 depending on the
background set and the multiple-testing denominator). And three design tables converting each
floor into a requirement in nuclei and donors.

## 2. What is not new, and is credited as such

- Cells within a donor are not independent observations: scDC 2019, scCODA 2021,
  propeller 2022, sccomp 2023; the expression-level version, Squair et al. 2021.
- Prospective parametric power for single-cell designs: scPower [2].
- Precision of quantification from scRNA-seq: Dai et al. 2025 [3], the closest prior art.

The manuscript states this explicitly at the head of the relevant Results section and uses no
unqualified novelty language anywhere (checked).

## 3. Borderline

The 190-contrast re-assessment. Its unconditional figures — 8% of pooled-significant contrasts
below the multinomial bound, 68% failing a donor-level test — are defensible measurements. The
larger 41.6% figure imports an overdispersion the paper itself shows does not transfer, and is
labelled a sensitivity. The pre-specified primary outcome returned nothing informative for a
reason visible in advance, and is reported as such.

## 4. Venue

**Nature Communications first; Genome Biology as the fallback.**

The 2026-09-22 assessment put NC as "borderline, not over" on two grounds, both now addressed:
Layer 3 rested on one cohort, and the strongest number depended on an imported constant. Layer 3
now has an independent replication, a working positive control and a dose–response with numbers;
the PsychAD decomposition supplies a positive result that borrows nothing. The paper also now
delivers something rather than only warning: three design tables and a threshold.

Against it: the work is a measurement and audit, negative-result-shaped, with no new method,
no software and no orthogonal validation. That is a real ceiling.

**Not Nature Methods** — no method is proposed and no tool delivered; its Analysis format
benchmarks methods, and this measures data. **Not Nature/Science/Cell/PNAS** — this changes how
existing results should be read; it is not a discovery or a technology.

Genome Biology remains the highest-probability home, the right readership, and `declarations.md`
is already in BMC format, so the fallback costs nothing but the NC review cycle.

## 5. Exposure, ranked by how much damage each could do

**A. The primary cohort's replicates span tissue sub-sampling.** Each cross-chip pair separates
a fresh tissue aliquot and an independent dissociation, so 4.28 is an upper bound on purely
technical variation and conflates within-region biological heterogeneity with technique. The
PsychAD decomposition answers this and should lead, rather than the bare 4.28.

**B. The Layer 3 external test is underpowered per stratum, and the ladder is non-monotone.**
⚠️ Partly closed 2026-09-24: the ladder was repeated in the primary cohort, where it is monotone
over the same range and where the 150-nucleus value reproduces the manuscript's 0.663 by an
independent implementation. Crossing points are now quoted from local interpolation in each
cohort rather than from a fit whose window was chosen post hoc.
1,200–3,000 links against the primary cohort's 33,227; promoter-proximal links number 4–49 per
run. The ladder fit is restricted to n ≥ 600 because the curve is non-monotone below that, and
that restriction was chosen after inspecting the data. The non-monotonicity itself is
unexplained. **This is the most attackable technical point in the paper.**

**C. At 150 nuclei the sign depends on annotation granularity.** Coarse classes give 0.956 —
absence of enrichment, not depletion — and only the granularity-matched strata give 0.656. The
defensible headline at that operating point is "no enrichment", and the text now says so.

**D. The full-depth reference range is not a like-for-like comparator.** 2.45–4.01 comes from
cellranger-arc link sets on full datasets: different algorithm, different multiple-testing
handling, vastly more nuclei. It carries weight in Fig. 5a and Fig. 6a. The nucleus ladder
largely defuses this — the same pipeline reaches that range at 5,000 nuclei — but the objection
will still be raised.

**E. κ non-transferability undercuts the re-assessment's strongest number**, leaving the
"so what" resting on the smaller unconditional figures. Already labelled in the text.

**F. A failed pre-specified primary outcome.** Reported honestly, which we would argue is a
credit; some reviewers will read it as a weakness regardless.

**G. Single seed** for the depth series and for the original per-class matched runs; the
nucleus ladder has three.

**H. Length.** 12,083 words. Roughly half is needed for an NC main text.

**I. Risk that this paper's own constant is misapplied.** The composition table's κ columns are
presented as a range with an instruction to measure locally, but a reader may lift 4.28 as a
universal constant — the precise error the paper warns against.

**J. Outstanding administrative items**, none of which we can supply: authors, affiliations,
ORCIDs, corresponding author, funding, ethics determination, competing interests, CRediT, and
whether to disclose generative-AI use. Plus a Zenodo DOI for the code.

## 6. What would most improve the odds, in order

1. Cut to ~5,000 words, methods and sensitivity analyses to supplementary.
2. Reframe Layer 1 to lead with the PsychAD decomposition rather than the raw overdispersion,
   which neutralises exposure A.
3. Add seeds to the depth series and, if cheap, resolve the ladder's non-monotone region below
   n = 600, which is exposure B.
