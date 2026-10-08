#!/usr/bin/env python3
"""Primary-cohort nucleus ladder, step 1: extract the nuclei.

Selects control nuclei of the target subtypes that clear the equalised depth
thresholds, then pulls their columns out of the ATAC and RNA count matrices.

The matrices are tens of gigabytes of MatrixMarket text, so each is streamed
once, line by line, keeping only the columns wanted. Loading them into memory
first is not an option on any machine this was run on.
"""
import argparse, gzip, json, logging, time
from pathlib import Path
import numpy as np, pandas as pd
from scipy import sparse

logger = logging.getLogger("extract")
ATAC_MIN, RNA_MIN = 5564, 5265

def scan(path, keep_cols, n_feat, tag):
    """Stream one MatrixMarket file, keeping only the wanted columns.

    Args:
        keep_cols: 1-based MTX column number -> new 0-based index.
        n_feat: row count of the output, i.e. the feature count.
        tag: label for the progress log.

    Returns:
        A CSC matrix of the kept columns.
    """
    # A dense lookup array rather than a dict: this is hit once per non-zero
    # entry, hundreds of millions of times, and the dict lookup dominates.
    # -1 marks a column that is not wanted. 180017 is the cell count of the
    # source matrices, so every valid 1-based column number indexes into it.
    remap = np.full(180017, -1, dtype=np.int64)
    for c, j in keep_cols.items(): remap[c] = j
    rows, cols, vals = [], [], []
    t0 = time.time(); n = 0
    with open(path) as fh:
        fh.readline(); fh.readline()   # MatrixMarket banner and dimension line
        for line in fh:
            n += 1
            if n % 200_000_000 == 0:
                logger.info("  %s %d 行 (%.0f 分)", tag, n, (time.time()-t0)/60)
            # Parse by hand rather than with split(): at this line count the
            # allocation of a list per line is the single largest cost.
            i = line.find(' '); j = line.find(' ', i+1)
            c = int(line[i+1:j])
            k = remap[c]
            if k < 0: continue
            v = int(line[j+1:])
            if v == 0: continue
            rows.append(int(line[:i])-1); cols.append(k); vals.append(v)
    logger.info("  %s 完成：%d 行扫描，%d 个非零保留 (%.0f 分)", tag, n, len(vals), (time.time()-t0)/60)
    return sparse.csc_matrix((vals, (rows, cols)), shape=(n_feat, len(keep_cols)), dtype=np.int32)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="Data/Other_Datasets/Multiome_Dataset/files")
    ap.add_argument("--out", default="results_methods/primary_ladder")
    ap.add_argument("--types", default="Exc_LINC00507_FREM3,Oligodendrocytes")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S")
    D = Path(a.data); O = Path(a.out); O.mkdir(parents=True, exist_ok=True)

    m = pd.read_csv(D/"Multiome_Dataset_Metadata.txt", sep="\t", low_memory=False,
                    usecols=["CellId","Case","WNN_L4","atac_peak_region_fragments","nCount_RNA"])
    m["row"] = np.arange(1, len(m)+1)          # MTX 列号 = 元数据行序（1-based）
    types = a.types.split(",")
    # Control nuclei only, of the requested subtypes, that already clear both
    # equalised depth thresholds. Filtering before the scan is what keeps the
    # kept-column set small enough for the streaming pass to be worth doing.
    sel = m[(m.Case == "HC") & (m.WNN_L4.isin(types)) &
            (m.atac_peak_region_fragments >= ATAC_MIN) & (m.nCount_RNA >= RNA_MIN)]
    logger.info("选中 %d 个核: %s", len(sel), sel.WNN_L4.value_counts().to_dict())
    keep = {int(r): j for j, r in enumerate(sel.row.to_numpy())}
    sel[["CellId","WNN_L4","atac_peak_region_fragments","nCount_RNA"]].to_csv(O/"cells.csv", index=False)

    A = scan(D/"Multiome_Dataset_ATAC_counts_raw.mtx", keep, 200791, "ATAC")
    sparse.save_npz(O/"atac.npz", A.tocsr()); logger.info("ATAC 已存 %s", A.shape)
    R = scan(D/"Multiome_Dataset_RNA_counts_raw.mtx", keep, 35367, "RNA")
    sparse.save_npz(O/"rna.npz", R.tocsr()); logger.info("RNA 已存 %s", R.shape)
    logger.info("EXTRACT COMPLETE")

if __name__ == "__main__":
    main()
