# The statistical unit of inference sets resolution limits in single-nucleus multiome analysis

Analysis code for *The statistical unit of inference sets resolution limits in
single-nucleus multiome analysis*.

Single-nucleus multiome studies routinely present three kinds of claim side by
side: that a cell type changes in proportion, that genes change within a cell
type, and that accessible regions are linked to genes within that type. Each
rests on a different statistical unit, and therefore has its own resolution
limit. This repository holds the code that measures all three, using
donor-matched process replicates from four cohorts.

## Scope of this repository

**Scripts only.** No data, no figures, no tables. Every input is either public
and documented below, or deposited with the manuscript as supplementary data.

That means nothing here runs end to end out of the box: each script needs
inputs you have to obtain first. What the repository gives you is the exact
procedure, including the conventions that turn out to matter and the mistakes
that were made and corrected along the way. Two of the three constants in the
paper do not transfer between datasets, so the procedure is the part worth
publishing.

## Getting the inputs

| Input | Where |
|---|---|
| GENCODE v32 annotation | `analysis/p5_05_download_refs.sh` |
| 10x Cell Ranger ARC reference link sets | `analysis/p5_18_download_arc_refs.sh` |
| Primary cohort (Ulm multiome) | from its depositors; see the manuscript's data availability statement |
| Seattle atlas (SEA-AD) | public; portal link in the manuscript |
| NABEC/HBCC prefrontal cortex | Catching et al., Cell Rep 2026 |
| PsychAD split-aliquot replicates | public |
| Derived tables the figures are built from | Supplementary Tables 1-21 of the manuscript |

The two download scripts are the only ones that fetch anything; everything
else reads from disk.

## Layout

```
analysis/            primary pipeline, run in numerical order
                     (p5_00 ... p5_21; there is no p5_11)
reassessment/        later analyses: PsychAD split aliquots, the Layer 3
                     external test, the two nucleus ladders, the 190-contrast
                     re-assessment, the three-arm linker comparison
figures/             figure generation (R / ggplot2), shared theme in
                     fig_common.R
technical_floor/     the packaged floor estimator (see its own README)
```

### Analysis pipeline

| Step | Produces |
|---|---|
| `p5_01_pseudobulk_stream.py` | pseudobulk by donor x chip x cell type |
| `p5_04_technical_replicates.py` | composition and expression discrepancies per replicate pair |
| `p5_05_download_refs.sh`, `p5_06_build_tss.py` | GENCODE v32 and the TSS table |
| `p5_07_extract_cells.py`, `p5_08_claimA_features.py` | equalised 150-nucleus substrate; peak-gene links |
| `p5_09`, `p5_10` | scaling fits for the expression and composition layers |
| `p5_12`-`p5_14` | Seattle atlas replication; matched null and stratified fits |
| `p5_15_promoter_OR_external.py` | promoter-enrichment odds ratios for external link sets |
| `p5_17_reviewer_response_stats.py` | recomputations from the deposited tables |
| `p5_20_calibration_donor_requirement.py` | how many replicate donors pin down the expression floor |
| `p5_21_kappa_abundance_model.py` | tests whether a constant-CV model explains the kappa-abundance slope |

`p5_16_figures.py` is a superseded matplotlib script kept only for provenance.
The published figures come from `figures/`.

### Re-assessment and external tests

| Script | What it does |
|---|---|
| `psychad_floors.py`, `psychad_L2.py` | composition and expression floors in the PsychAD split-aliquot donors |
| `floor_reassessment.py` | 190 deposited composition contrasts against the floor |
| `nabec_layer3.py` | Layer 3 in the external cohort; `--log1p --min-det-frac --fdr-all-pairs` align it with `p5_08` |
| `nabec_L3_report.py`, `nabec_L3_ladder_fit.py` | Mantel-Haenszel pooling and crossing-point estimates |
| `primary_ladder_extract.py`, `primary_ladder_run.py` | the same nucleus ladder inside the primary cohort |
| `primary_ladder_aggregated.py` | that ladder under three linking procedures: single-nucleus correlation, aggregation with ArchR-like retrieval (`--archr-defaults`), aggregation with FDR only |
| `c1_run_grid.py`, `c1_summary.py` | drive that comparison across the ladder and pool it |
| `c1_depth_grid.py`, `c1_depth_summary.py` | the nuclei x depth surface at constant product |
| `c1_trans_null_validation.py`, `c1_trans_null_summary.py` | stability of the trans-pairing null over nine reassignments and a leave-one-chromosome-out jackknife |

### Figures

```bash
Rscript figures/fig2_fig3_ggplot2.R      <derived_results> <out>
Rscript figures/fig6_ggplot2.R           <derived_results> <out>
Rscript figures/fig7_ggplot2.R           <supplementary_tables> <out>
Rscript figures/supp_figs_ggplot2.R      <derived_results> <out>
# Figures 1, 4 and 5 additionally need the primary cohort's donor metadata:
Rscript figures/fig1_fig4_fig5_ggplot2.R <derived_results> <out> <metadata_dir>
```

Supply `<derived_results>` and `<supplementary_tables>` from the manuscript's
supplementary data. Requires R >= 4.5 with `ggplot2`, `dplyr`, `readr`,
`patchwork`, `scales`, `ragg` and `tidyr`. Every label is Arial and no panel
carries an in-figure title, per the journal's figure requirements.

Two rendering faults are worth knowing about, since both fail silently and both
are handled in `fig_common.R`:

* Rscript inherits the caller's locale, and under the C locale R transliterates
  any character outside Latin-1 to `..` on its way to the graphics device. A
  kappa becomes two periods in both the PNG and the PDF, with no warning.
* A theme's `base_family` does not reach text drawn by a geom, which keeps the
  device default and leaves a second font embedded in the PDF.

## Measure the limits in your own data

The floor-estimation procedure is packaged as `technical_floor`. Two of the
three constants reported in the paper do not transfer between datasets, which
is precisely why the tool exists: what travels is the procedure, not a table of
numbers.

```bash
python3 -m technical_floor \
  --composition your_composition_pairs.csv \
  --expression  your_expression_floors.csv \
  --out results/
```

It takes two tidy CSVs of replicate pairs and returns kappa per cell type, the
expression scaling fit, and the design tables derived from both.
`technical_floor/README.md` documents the input formats. Requires only `numpy`
and `pandas`.

The paper's own per-pair tables are deposited as Supplementary Tables 2, 3 and
10, and running the tool on them reproduces every published value.

Read `composition_kappa_by_celltype.csv` before the pooled kappa. Across cell
types kappa varied 5.6-fold within the primary cohort against 1.1-fold between
cohorts, and rose with abundance (r = 0.87); because cost scales as kappa
squared, that is a 31-fold range in required nuclei within one dataset. The
pooled value describes no individual cell type.

## Three definitions that are easy to confuse

**Composition floor.** For a replicate pair and cell type with proportions
f1, f2 and pooled proportion p: `z = |f1 - f2| / sqrt(p(1-p))`. The multinomial
expectation is `sqrt(2/n_eff)` with `n_eff = 2/(1/n1 + 1/n2)`. Overdispersion is
the **root mean square** across pairs of `z / sqrt(2/n_eff)`, not a median.

**Expression floor.** Per pair and cell type, the median over genes of
`|log2(CPM1+1) - log2(CPM2+1)|`, over genes with CPM > 0 in **either** member
(union, >= 200 genes). The matched null resamples both members multinomially
from their pooled profile at the observed depths.

**Promoter-enrichment odds ratio.** The odds of promoter-proximal
(`|peak midpoint - TSS| <= 3 kb`) among detected links against the same odds
among all eligible (peak, TSS) pairs in the window, Haldane-Anscombe corrected.
This is **not** the TSS enrichment score used in ATAC quality control: that
score is computed on reads and asks whether a library captured open chromatin
at promoters, while this one is computed on inferred peak-gene pairs. A library
can have an excellent TSS enrichment score and still produce a link set with an
odds ratio at or below one, which is what the paper reports.

The odds ratio depends on two conventions that must be reported with it: which
pairs constitute the background, and whether pairs that cannot be called (peaks
with no variance at the depth analysed) enter the denominator of the
multiple-testing correction. Holding nuclei, depth and code otherwise fixed,
those two choices alone move the ratio between 0.66 and 1.5.

## Citation

Guo, K. The statistical unit of inference sets resolution limits in
single-nucleus multiome analysis. *Manuscript in preparation* (2026).

## License

MIT, see `LICENSE`. The primary datasets the analyses read remain under the
licences of their respective depositors.
