# Environment

## Python

Tested on Python 3.11. Install with `uv` (preferred) or pip:

```bash
uv venv && uv pip install -r requirements.txt
# or
python3 -m venv .venv && . .venv/bin/activate && pip install -r requirements.txt
```

`anndata` and `h5py` are needed only by the scripts that stream `.h5ad` input
(`reassessment/nabec_layer3.py`, `reassessment/psychad_*.py`). The scripts that
stream MatrixMarket files (`analysis/p5_01`, `analysis/p5_07`,
`reassessment/primary_ladder_extract.py`) need only numpy, scipy and pandas.

## R

Tested on R 4.5.

```r
install.packages(c("ggplot2", "dplyr", "readr", "patchwork",
                   "scales", "ragg", "tidyr"))
```

`ragg` is used for PNG output at 450 dpi. Figures are written at 180 mm width
for a two-column layout.

## Memory and runtime

Most steps run in minutes on a laptop. Two are heavier:

- `analysis/p5_01_pseudobulk_stream.py` and
  `reassessment/primary_ladder_extract.py` stream a 30 GB MatrixMarket file in a
  single pass. Peak memory stays under ~4 GB because nothing is held dense;
  the ATAC pass takes about 10 minutes on an SSD.
- `reassessment/nabec_layer3.py` at 5,000+ nuclei computes peak–gene
  correlations by sparse matrix–vector products rather than densifying the
  n × 521,217 accessibility matrix. Densifying it would need roughly 10 GB.

## Reproducibility

Every script that samples takes a `--seed`. The paper reports three seeds for
each configuration of the nucleus ladders and of the equalised design. Set
`PYTHONHASHSEED` if you need bit-identical dictionary ordering across runs.
