# Data sources

No primary data are redistributed in this package. Everything below is publicly
available (or, where marked ⚠️, must be supplied by the authors before submission).
`derived_results/` holds the per-pair and per-contrast tables our own analyses produced;
those are the files the Supplementary Tables are drawn from.

---

## 1. Primary cohort — ALS motor cortex single-nucleus multiome

Used for: all three floors; the 26 cross-chip technical replicate pairs.

| | |
|---|---|
| Content | 79 donors, 180,016 nuclei passing the original authors' QC, 12 chips, 32 wells |
| Donor cohorts | Boston, Edinburgh, Hannover, Netherlands Brain Bank |
| Groups | ALS, ALS-FTD, healthy control (`Case` field) |
| Annotation | `WNN_L1` … `WNN_L4` |
| Files used | `Multiome_Dataset_Metadata.txt`, `Multiome_Dataset_RNA_counts_raw.mtx` (35,367 × 180,016; 665,652,320 stored entries), `Multiome_Dataset_ATAC_counts_raw.mtx` (200,791 × 180,016), plus barcodes/features |
| Local path | `ResearchD/Data/Other_Datasets/Multiome_Dataset/files/` (10 files, 37 GB) |
| **Source** | **Zenodo https://zenodo.org/records/18370646 — doi:10.5281/zenodo.18370646** |
| Deposited | 2026-01-29 (v1), 10 files, 40.2 GB, CC-BY-4.0 |
| Depositors | Veselin Grozdanov (contact); Karin Danzer (project lead), Ulm |
| Title | *Multi-omic single-nucleus ATAC-seq + RNA-seq of nuclei from the post-mortem human primary motor cortex from patients with ALS, ALS-FTD, C9ORF72-ALS-FTD and neurologically unaffected controls* |

**Primary publication:** Ruf WP, Kühlwein JK, Meier L, et al. Multi-modal dissection of
cell-type specific TDP-43 pathology in the motor cortex. *Nat Commun* 2026;17:2406.
doi:10.1038/s41467-026-69944-6. Code: github.com/DanzerLab/Ruf_et_al_2026.

**Provenance of the cross-chip pairs — resolved 2026-09-22** from the paper's Supplementary
Data 3 (`Well_Layout`), verbatim:

> "Two to four samples' nuclei were pooled during preparation, processed together on a well
> of the chip and demultiplexed in silico after sequencing based on sex ... and
> single-nucleotide polymorphisms (SNPs). **29/79 samples were processed on different 10X
> Genomics chip and wells on different experimental runs.**"

Parsing that table: 79 samples across 12 chips, 109 sample-by-well rows; 29 samples appear in
≥2 wells and **27 of those span ≥2 chips**, matching the manuscript exactly. Example: ALS21–24
were pooled as one group twice, on Chip8/Well1 and again on Chip12/Well1.

So each pair separates a fresh tissue aliquot, an independent nuclei isolation (mechanical
homogenisation, OptiPrep gradient), pooling, loading, library preparation and sequencing run.
This confirms the manuscript's "technical by construction" claim and, if anything,
understates it. ⚠️ It also means the floor includes within-region tissue sub-sampling, so it
is an upper bound on purely technical variation — worth stating in the text.

## 2. External replication cohort — SEA-AD, middle temporal gyrus

Used for: replication of the composition and expression floors (80 donors with ≥2
`library_prep`); the Layer 3 linkage-saturation panel (28 multiome libraries).

| | |
|---|---|
| Content | 84 donors, 1,178,694 nuclei, 205 library preparations |
| Annotation | `Class` (3), `Subclass` (24), `Supertype` (148) |
| Portal | https://sea-ad.org |
| Browser | CZ CELLxGENE collection `1ca90a2d-2943-483d-b678-b809bf464c30` |
| Citation | Gabitto MI, Travaglini KJ, et al. *Nat Neurosci* 2024. doi:10.1038/s41593-024-01774-5 |
| Files used | `SEAAD_MTG_RNAseq_final-nuclei.2026-06-22.h5ad` (raw counts in `layers/UMIs`), `SEAAD_MTG_RNAseq_final-nuclei_metadata.2026-06-22.csv` |
| Local path | `ResearchD/Data/Other_Datasets/SEAAD_MTG/files/` |

## 3. Reference peak–gene link sets (Layer 3 diagnostic)

Four public 10x Genomics multiome demonstration datasets, processed with Cell Ranger ARC
2.0.0, distributed under CC BY 4.0 from https://www.10xgenomics.com/datasets

Reported at ±500 kb, the window of the primary analysis. The background and
the detected link set are restricted to the same window; the earlier ±1 Mb
figures bounded only the background, which is why they differ. Both windows
are deposited, under `derived_results/window_500000/` and `window_1000000/`;
regenerate with `code/analysis/p5_18_download_arc_refs.sh` then
`p5_19_recompute_refs_at_window.sh`. | Dataset | Peaks | links used (±500 kb)
| OR at ±500 kb | OR at ±1 Mb | |---|---|---|---|---| | `human_brain_3k` |
133,986 | 105,471 | 3.558 [3.465, 3.660] | 4.054 [3.947, 4.168] | |
`pbmc_granulocyte_sorted_3k` | 98,290 | 61,520 | 2.723 [2.622, 2.833] | 3.087
[2.974, 3.210] | | `pbmc_granulocyte_sorted_10k` | 143,836 | 274,818 | 2.261
[2.212, 2.312] | 2.463 [2.410, 2.518] | | `lymph_node_lymphoma_14k` | 109,714
| 82,732 | 3.127 [3.031, 3.231] | 3.596 [3.486, 3.715] |

Local path: `ResearchD/Data/Other_Datasets/TenX_Multiome_Examples/files/`

⚠️ **Two acquisition pitfalls, both encountered here:**

1. The standalone `human_brain_3k_feature_linkage.bedpe` offered for download can fail
   silently and be saved as a 111-byte S3 `AccessDenied` XML page. Take the file from
   `analysis/feature_linkage/feature_linkage.bedpe` **inside the analysis tarball** instead
   (847,032 lines).
2. The complete peak set must come from the **Feature ID column of the per-cluster
   differential accessibility output**, not from the transcription-factor motif mapping file,
   which omits peaks without motif matches and would bias the background.

Prepared inputs are regenerated by the shell block recorded in
`protocol/2026-09-22_methods_written_from_original_code.md`.

### Fifth link set — dropped 2026-09-24

A fifth external link set (OR 5.40 [5.30, 5.51]) appeared in an earlier draft, attributed to a
published cell-type-resolved Alzheimer's disease multiome study. No citation, no input file and
no code capable of recomputing it exists anywhere in the project, and `derived_results/` holds
only the four sets above. A number no reader could reproduce has no place in the paper, so it
was removed rather than reconstructed; the manuscript now reports four link
sets, and the quoted range is the range of those four, 2.26–3.56 at the ±500
kb window of the primary analysis.

## 4. Datasets used only in the floor re-assessment

| ID | Accession / source | Tissue | Contrast | Citation |
|---|---|---|---|---|
| D3 | GEO **GSE290359** (BA44/BA46) | brain | disease vs control | — |
| D4 | GEO **GSE330130** (mid-temporal cortex) | brain | ALS vs non-neurological control | — |
| D5 | GEO **GSE158055** | blood (PBMC) | severe/critical vs control | Ren X, et al. *Cell* 2021 |

GEO records: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<ACCESSION>`

The lumbar spinal cord component of GSE330130 was examined and excluded: its deposited `obs`
carries no diagnosis field.

## 4b. External test of Layer 3 — NABEC/HBCC prefrontal cortex multiome

Surveyed and excluded as a source of technical replicates (§6), this resource is
nonetheless the only public human brain dataset with **paired ATAC and gene expression in
the same nucleus** at a scale that permits the Layer 3 question to be asked at cell numbers
far above the equalised design of the primary cohort.

| | |
|---|---|
| Content | 1,501,089 nuclei, 357 samples, 7 cell classes; 38,606 genes, 521,217 peaks |
| Pairing | RNA and ATAC barcodes match 1:1 for all 1,501,089 nuclei once the duplicated sample suffix is stripped from the ATAC index (`…-1_HBCC-1193_HBCC-1193` → `…-1_HBCC-1193`) |
| Nuclei per class | Oligo 610,235; ExN 373,139; InN 197,917; Astro 128,599; MG 87,384; OPC 86,659; VC 17,156 |
| Depth | ATAC 5,615–18,588 unique fragments (median by class); RNA ≈4,300–4,900 counts |
| Source | Zenodo https://zenodo.org/records/18394349 (processed multiome matrix, open) |
| Code | github.com/NIH-CARD/scMAVERICS; doi:10.5281/zenodo.18135365 |
| Raw sequencing | dbGaP `phs004202.v1.p1` (NABEC) and NDA collection #3151 (HBCC), both controlled — not needed here |
| Local files | `Data/final_rna_data.h5ad`, `Data/final_atac_data.h5ad` |

**Donor metadata** (`Data/TableS1_CohortDemographics (1).xlsx`, 362 rows) links to the
matrices once HBCC identifiers are prefixed: TableS1 writes them as bare numbers while the
deposited `SampleID` writes `HBCC-<number>`. After that normalisation all 357 analysed
samples match, and the derived demographics reproduce the published values exactly —
NABEC 202 donors, mean age 48.73 (reported 48.8), mean PMI 12.34 (reported 12.34); HBCC 155
donors, mean age 45.16 (reported 45.2), mean PMI 33.79 (reported 33.79). The 5 unmatched
TableS1 rows are the samples dropped at quality control, which is what separates the 362
collected from the 357 analysed.

## 5. Gene annotation

GENCODE release 32 — the annotation matching the 10x GRCh38-2020-A reference used by the
primary cohort.

```
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz
```

The derived TSS table (`derived_results/p5_tss_gencode_v32.csv`, 35,334 genes) is included.
Citation: Frankish A, et al. *Nucleic Acids Res* 2019;47(D1):D766–D773.

## 6. Resources examined in the replicate survey and excluded

These were surveyed, not analysed. None contained within-donor cross-batch replicates.

| Resource | Why excluded | Citation |
|---|---|---|
| Prefrontal cortex multiome atlas (NABEC/HBCC) | one library per donor — no replicate structure, so it cannot measure a floor. ⚠️ Excluded **as a replicate source only**; it is used as the external test of Layer 3 and as the substrate for the nucleus ladder — see §4b | Catching A, et al. *Cell Rep* 2026;45(3):117110 |
| SEA-AD multiome subset | 28 libraries from 28 donors | Gabitto MI, et al. (as above) |
| ROSMAP multi-region | one library per donor × region | Mathys H, et al. *Nature* 2024;632:858–868 |
| C9orf72 ALS/FTD multiome (GSE212630) | 14 donors, below pre-specified minimum | Wang H-LV, et al. *PNAS* 2025;122(9):e2419818122 |

Resolved 2026-09-22 against the primary paper and its Table S1: the collected set is
**362** (204 NABEC + 158 HBCC) and the analysed set is **357** (202 + 155). The manuscript
now uses 357, matching the cited report.

---

## Derived results included here

`derived_results/` contains the intermediate tables our analyses produced. They are the
direct upstream of the Supplementary Tables and allow every figure to be regenerated without
re-processing any primary data.

| File | Content |
|---|---|
| `p5_fig1_donor_table.csv` | donor, chip, well assignment |
| `p5_L1_scaling_pairs.csv`, `seaad_L1_pairs.csv` | per-pair composition discrepancies |
| `p5_floor_scaling_pairs.csv`, `seaad_L2_obs.csv`, `seaad_L2_null.csv` | per-pair expression floors with matched null |
| `p5_floor_scaling_boot.csv`, `p5_two_layer_boot.csv` | bootstrap replicates of the scaling exponents |
| `seaad_L2_by_celltype.csv` | within-cell-type exponents |
| `window_500000/p5_promoterOR_*.csv` | promoter-enrichment odds ratios, four
external link sets, at the ±500 kb window of the primary analysis (also
`window_1000000/`; the superseded inconsistent-window originals are in
`superseded_inconsistent_window/`) |
| `seaad_L3_libraries.csv` | feature linkages per external multiome library |
| `p5_claimA_features.csv` | per-subtype independent link counts (equalised design) |
| `p5_tss_gencode_v32.csv` | TSS table derived from GENCODE v32 |
| `floor_reassessment_contrasts.csv` | all 190 composition contrasts against the floor |
| `nabec_L3_external_replication.csv` | Layer 3 external replication, 7 classes + 4 fine clusters + depth series |
| `nabec_L3_nucleus_ladder.csv` | Layer 3 nucleus ladder, external cohort, 10 levels × 2 cell types × 3 seeds |
| `primary_L3_nucleus_ladder.csv` | Layer 3 nucleus ladder, primary cohort, 6 levels × 2 cell types × 3 seeds |
