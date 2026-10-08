#!/usr/bin/env python3
"""
Cross-chip technical replicate calibration
=========================================
The Ulm design contains a resource its own authors did not use, and that this
project very nearly missed:

    12 chips / 32 wells
    52 donors appear on one chip; 27 donors appear on TWO

One person, two chips: every difference between the two measurements is
technical. That is a free floor on technical noise, and any disease effect
smaller than it cannot be believed.

Two layers are calibrated here.

  L1, composition -- metadata only, no matrix computation.
      The difference in cell-type proportion for one donor across two chips.
      Gives a noise floor for every proportion-based conclusion.

  L2, expression -- needs pseudobulk grouped by (donor x chip x cell type).
      Expression correlation and variability for one donor across two chips.
      Gives a noise floor for every differential expression or accessibility
      conclusion. Run p5_01_pseudobulk_stream.py --group-extra Well first.

Why this matters: the source paper itself notes that DEG counts correlate
strongly with total reads per cell type, so comparing the magnitude of change
between cell types is meaningless without power balancing. A measured technical
floor is the yardstick that makes that balancing possible.

用法：
  python3 p5_04_technical_replicates.py --data DIR --out results/
  python3 p5_04_technical_replicates.py --data DIR --out results/ \
      --pseudobulk results/pseudobulk_RNA_WNN_L25_byWell.npz     # 加做 L2
"""

import argparse
import itertools
import os
import sys

import numpy as np
import pandas as pd


def load_meta(data_dir, level):
    m = pd.read_csv(os.path.join(data_dir, "Multiome_Dataset_Metadata.txt"),
                    sep="\t", usecols=["CellId", "ID", "Well", "Case", level],
                    dtype=str)
    m["Chip"] = m["Well"].str.extract(r"^(Chip\d+)")
    return m


def composition_level(m, level, min_cells=50):
    """L1: cross-chip proportion differences for one donor.

    This is the technical noise floor for any composition analysis.
    """
    per = (m.groupby(["ID", "Chip", level]).size()
             .rename("n").reset_index())
    tot = per.groupby(["ID", "Chip"])["n"].sum().rename("tot").reset_index()
    per = per.merge(tot, on=["ID", "Chip"])
    per = per[per["tot"] >= min_cells]
    per["frac"] = per["n"] / per["tot"]

    # Donors appearing on two or more chips -- the replicate pairs
    nchip = per.groupby("ID")["Chip"].nunique()
    reps = nchip[nchip >= 2].index.tolist()
    print(f"  上了 >=2 个芯片且每片 >= {min_cells} 核的供体: {len(reps)}")
    if not reps:
        sys.exit("[FATAL] 没有可用的技术重复。")

    rows = []
    for did in reps:
        sub = per[per.ID == did]
        chips = sorted(sub.Chip.unique())
        for c1, c2 in itertools.combinations(chips, 2):
            a = sub[sub.Chip == c1].set_index(level)["frac"]
            b = sub[sub.Chip == c2].set_index(level)["frac"]
            ct = a.index.union(b.index)
            a = a.reindex(ct).fillna(0)
            b = b.reindex(ct).fillna(0)
            for t in ct:
                rows.append(dict(ID=did, chip1=c1, chip2=c2, celltype=t,
                                 f1=a[t], f2=b[t],
                                 absdiff=abs(a[t] - b[t]),
                                 logratio=np.log2((a[t] + 1e-3) / (b[t] + 1e-3))))
    d = pd.DataFrame(rows)

    summ = (d.groupby("celltype")
              .agg(n_pairs=("absdiff", "size"),
                   mean_frac=("f1", "mean"),
                   median_absdiff=("absdiff", "median"),
                   p90_absdiff=("absdiff", lambda x: x.quantile(0.90)),
                   median_abs_log2FC=("logratio", lambda x: x.abs().median()),
                   p90_abs_log2FC=("logratio", lambda x: x.abs().quantile(0.90)))
              .reset_index()
              .sort_values("mean_frac", ascending=False))
    return d, summ


def expression_level(pb_path, m, level):
    """L2: cross-chip expression agreement for one donor.

    This is the technical noise floor for any differential expression result.
    """
    z = np.load(pb_path, allow_pickle=True)
    counts, groups = z["counts"], [str(g) for g in z["groups"]]
    parts = [g.split("||") for g in groups]
    if len(parts[0]) < 3:
        sys.exit("[FATAL] 伪批量不是按 (供体×Well×细胞类型) 分组的。"
                 "请用 --group-extra Well 重跑 p5_01。")
    gi = pd.DataFrame(parts, columns=["donor", "well", "celltype"])
    gi["chip"] = gi["well"].str.extract(r"^(Chip\d+)")
    gi["col"] = np.arange(len(gi))

    # CPM normalise then log, so the comparison is on the scale DE tests use
    lib = counts.sum(0).astype(float)
    ok = lib > 0
    X = np.zeros_like(counts, dtype=np.float32)
    X[:, ok] = np.log1p(counts[:, ok] / lib[ok] * 1e6)

    rows = []
    for (donor, ct), sub in gi.groupby(["donor", "celltype"]):
        chips = sub.chip.dropna().unique()
        if len(chips) < 2:
            continue
        for c1, c2 in itertools.combinations(sorted(chips), 2):
            i = sub[sub.chip == c1].col.values
            j = sub[sub.chip == c2].col.values
            if len(i) == 0 or len(j) == 0:
                continue
            x, y = X[:, i[0]], X[:, j[0]]
            keep = (x > 0) | (y > 0)
            if keep.sum() < 200:
                continue
            r = np.corrcoef(x[keep], y[keep])[0, 1]
            rows.append(dict(donor=donor, celltype=ct, chip1=c1, chip2=c2,
                             pearson_logCPM=r, n_genes=int(keep.sum()),
                             median_abs_log2FC=float(
                                 np.median(np.abs(x[keep] - y[keep])) / np.log(2))))
    d = pd.DataFrame(rows)
    if d.empty:
        return d, d
    summ = (d.groupby("celltype")
              .agg(n_pairs=("pearson_logCPM", "size"),
                   median_r=("pearson_logCPM", "median"),
                   min_r=("pearson_logCPM", "min"),
                   noise_floor_log2FC=("median_abs_log2FC", "median"))
              .reset_index().sort_values("median_r", ascending=False))
    return d, summ


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", default="results")
    ap.add_argument("--level", default="WNN_L2.5")
    ap.add_argument("--min-cells", type=int, default=50)
    ap.add_argument("--pseudobulk", default=None,
                    help="按 (供体×Well×细胞类型) 分组的伪批量 npz，加做 L2")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    print("[1/3] 读元数据 …")
    m = load_meta(args.data, args.level)
    d0 = m.drop_duplicates(["ID", "Chip"])
    print(f"  {m.ID.nunique()} 供体 / {m.Chip.nunique()} Chip / {m.Well.nunique()} Well")
    print("  每供体的芯片数分布:", d0.groupby("ID").size().value_counts().to_dict())

    print("[2/3] L1 组成层标定 …")
    d1, s1 = composition_level(m, args.level, args.min_cells)
    d1.to_csv(os.path.join(args.out, "p5_techrep_composition_pairs.csv"), index=False)
    s1.to_csv(os.path.join(args.out, "p5_techrep_composition_summary.csv"), index=False)

    print("\n" + "=" * 86)
    print("L1 · 组成层技术噪声下限（同供体不同芯片）")
    print("=" * 86)
    print(f"{'细胞类型':<20}{'平均占比':>10}{'|Δ占比|中位':>13}{'|Δ占比|P90':>12}"
          f"{'|log2FC|中位':>13}{'|log2FC|P90':>12}")
    for _, r in s1.iterrows():
        print(f"{r.celltype:<20}{r.mean_frac:>10.3f}{r.median_absdiff:>13.4f}"
              f"{r.p90_absdiff:>12.4f}{r.median_abs_log2FC:>13.3f}"
              f"{r.p90_abs_log2FC:>12.3f}")
    print("\n解读：任何小于 |log2FC| P90 的组成学效应，都在技术噪声之内，不可报告为发现。")

    if args.pseudobulk:
        print("\n[3/3] L2 表达层标定 …")
        d2, s2 = expression_level(args.pseudobulk, m, args.level)
        if d2.empty:
            print("  无可用的跨芯片表达对。")
        else:
            d2.to_csv(os.path.join(args.out, "p5_techrep_expression_pairs.csv"), index=False)
            s2.to_csv(os.path.join(args.out, "p5_techrep_expression_summary.csv"), index=False)
            print(f"{'细胞类型':<20}{'配对数':>8}{'r中位':>9}{'r最小':>9}{'噪声下限log2FC':>16}")
            for _, r in s2.iterrows():
                print(f"{r.celltype:<20}{int(r.n_pairs):>8}{r.median_r:>9.3f}"
                      f"{r.min_r:>9.3f}{r.noise_floor_log2FC:>16.3f}")
    else:
        print("\n[3/3] 跳过 L2（未提供 --pseudobulk）。先跑：")
        print("  python3 code/p5_01_pseudobulk_stream.py --data DIR --out results \\")
        print("      --modality RNA --level WNN_L2.5 --group-extra Well")

    print(f"\n[done] 输出写入 {args.out}/")


if __name__ == "__main__":
    main()
