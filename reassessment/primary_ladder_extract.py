#!/usr/bin/env python3
"""原队列核数梯队 · 第 1 步：抽取目标亚型中通过等化阈值的对照核，单遍流式扫描 MTX。"""
import argparse, gzip, json, logging, time
from pathlib import Path
import numpy as np, pandas as pd
from scipy import sparse

logger = logging.getLogger("extract")
ATAC_MIN, RNA_MIN = 5564, 5265

def scan(path, keep_cols, n_feat, tag):
    """keep_cols: 1-based 列号 -> 新下标。返回 CSC。"""
    remap = np.full(180017, -1, dtype=np.int64)
    for c, j in keep_cols.items(): remap[c] = j
    rows, cols, vals = [], [], []
    t0 = time.time(); n = 0
    with open(path) as fh:
        fh.readline(); fh.readline()
        for line in fh:
            n += 1
            if n % 200_000_000 == 0:
                logger.info("  %s %d 行 (%.0f 分)", tag, n, (time.time()-t0)/60)
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
