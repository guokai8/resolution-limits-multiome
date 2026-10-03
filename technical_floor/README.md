# technical_floor

Measure the resolution limits of your own single-nucleus data from replicate libraries,
and convert them into a study design.

The paper this accompanies reports that two of the three constants do **not** transfer
between datasets. That is the reason this package exists: if the constants had transferred,
a table of numbers would have been enough. They do not, so what travels is the procedure.

## Install

No installation. The package has two dependencies, `numpy` and `pandas`.

```bash
export PYTHONPATH=code
python3 -m technical_floor --help
```

## Input

Two tidy CSVs, either of which may be given alone.

**Composition** — one row per (donor, cell type) replicate pair, in either form:

| form | columns | when to use |
|---|---|---|
| raw counts | `donor, celltype, n1, n2` | starting from your own pipeline |
| precomputed | `donor, celltype, z, n_eff` | reproducing a published analysis |

Prefer raw counts. The precomputed form exists because the denominator of a cell-type
fraction depends on which cell types entered the analysis, so reconstructing fractions
from counts does not always match what was published.

**Expression** — one row per (donor, cell type) pair:

`donor, celltype, n_eff, floor` (plus optional `n1, n2`)

`floor` is the median over union genes of |log2(CPM+1) difference| between the two
libraries. Build it from pseudobulk counts with `floor_from_counts`. **Supply `n1` and
`n2` if you have them**: without them the tool cannot apply the minimum-nuclei pair
criterion to this layer, and will say so. That gap is not hypothetical — it is the defect
this package was written after finding in the accompanying paper's own analysis.

## Run

```bash
python3 -m technical_floor \
    --composition pairs.csv \
    --expression floors.csv \
    --out results/
```

Useful flags: `--min-nuclei` (pair criterion, default 50), `--min-n-eff` (support floor for
the scaling fit, default 0), `--abundance` and `--target-pp` (what the design table is
asked to resolve).

## Output

- `composition_kappa_by_celltype.csv` — κ per cell type. **Read this before the pooled value.**
  κ varied 5.6-fold across cell types in the paper's own cohort and rose with abundance, so
  the pooled figure describes no individual cell type.
- `composition_design_table.csv` — resolvable difference by donors × nuclei, at the
  multinomial bound, the pooled κ, the median cell type and the worst.
- `composition_requirement.csv` — the same read backwards: nuclei per group for a target.
- `expression_design_table.csv` — |log2 fold change| thresholds, flagged where extrapolated.

## Worked example

The paper's own deposited tables are in `data/example_inputs/`. Running on them reproduces
every published number, which is also the acceptance test:

```
κ = 4.268 (95% CI 2.98–5.66), 300 pairs / 26 donors
  per cell type 1.53 to 8.50, r = +0.87 against log abundance
  nuclei per group for 1 percentage point:
      multinomial bound              6,915
      pooled κ                     125,949
      median cell type              44,435
      most overdispersed           499,575
expression floor = 4.58 × n^-0.497, 300 observations / 26 donors
  |log2FC| threshold 0.93 at n=25 … 0.15 at n=1,000
```

## Tests

```bash
python3 -m unittest tests.test_technical_floor -v
```

24 tests. One of them, `test_the_pair_criterion_changes_the_answer`, asserts that skipping
the minimum-nuclei criterion returns the superseded −0.507 rather than the published
−0.497 — it exists so that the failure mode cannot come back silently.

## What this does not do

Layer 3, regulatory inference, needs the raw peak × cell matrices and a linking algorithm,
so it is out of scope here. The paper's promoter-enrichment diagnostic runs separately from
a link table and a TSS annotation.
