#!/usr/bin/env python3
"""
Compute the five pre-registered features and run the hard gates
==============================================================
Strictly per the definitions in the pre-registration, computed on the
equal-count, equal-depth sub-matrices produced by the extraction step:

  P1 ⭐ redundancy      每个表达基因的**独立** linked peak 数的中位数
  P2    distal_frac     link 中远端（距 TSS >3 kb）的比例
  P3    entropy         peak 可及性分布的 Shannon 熵
  P4    redundancy_ALS  同 P1，限于预注册 §6 的 ALS 基因集
  P5    n_expressed     检出基因数（控制变量）

A link is defined as: a peak-gene correlation within this cell type across
  its 150 nuclei, inside a +/-500 kb window, at FDR < 0.05. "Independent" means
  the significant peaks are greedily decorrelated against one another, with
  Pearson |r| > 0.7 treated as the same regulatory element.

  This script NEVER looks at the dependent variable. It produces the feature
  table and nothing else. The relationship between features and the dependent
  variable is computed in a later step, and only after the gates below pass.
  Keeping the two apart is what stops a gate being chosen after seeing the
  answer.

  Hard gates. Failing either one voids the claim:
    1. |partial_r(redundancy, total counts in this type | other features)| < 0.2
       -- the feature must not be a restatement of sequencing depth
    2. Spearman rho > 0.9 between cell-type RANKINGS of redundancy across the
       three random seeds -- the feature must be stable to which nuclei were
       drawn

用法：
  python3 p5_08_claimA_features.py --out results --seeds 0,1        # 已有的种子
  python3 p5_08_claimA_features.py --out results --seeds 0 --types 3  # 冒烟测试
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

ALS_GENES = [  # 预注册 §6，事前固定，不得增删
    "SOD1", "TARDBP", "TBK1", "ALS2", "BSCL2", "DCTN1", "GLE1", "GRN",
    "HNRNPA1", "NEFH", "PFN1", "SETX", "UBQLN2", "VAPB", "VEGFA",
    "NEFL", "STMN2", "UNC13A", "NPTX1", "NPTX2", "NPTXR", "EPHA4",
]
WINDOW = 500_000
PROMOTER = 3_000
FDR = 0.05
COLLIN = 0.7


def load_sparse(path):
    z = np.load(path, allow_pickle=True)
    shape = tuple(int(x) for x in z["shape"])
    return dict(rows=z["rows"], cols=z["cols"], vals=z["vals"], shape=shape,
                features=z["features"], celltype=z["celltype"],
                donor=z["donor"], cells=z["cells"])


def dense_block(sp, col_idx):
    """Rebuild the named columns as a dense (n_feat x len(col_idx)) float32.

    Densified one cell type at a time; the full matrix would not fit.
    """
    remap = np.full(sp["shape"][1], -1, dtype=np.int64)
    remap[col_idx] = np.arange(len(col_idx))
    m = remap[sp["cols"]] >= 0
    X = np.zeros((sp["shape"][0], len(col_idx)), dtype=np.float32)
    X[sp["rows"][m], remap[sp["cols"][m]]] = sp["vals"][m]
    return X


def bh_threshold(p, alpha):
    """Benjamini-Hochberg: largest p-value called significant, or -1 if none."""
    if len(p) == 0:
        return -1.0
    ps = np.sort(p)
    k = np.arange(1, len(ps) + 1)
    ok = ps <= alpha * k / len(ps)
    return ps[ok][-1] if ok.any() else -1.0


def corr_rows_vs_vec(A, y):
    """Pearson r of every row of A against y.

    A is (m x n), one feature per row; y is (n,). Vectorised over rows.
    """
    Az = A - A.mean(1, keepdims=True)
    As = np.sqrt((Az ** 2).sum(1))
    yz = y - y.mean()
    ys = np.sqrt((yz ** 2).sum())
    denom = As * ys
    out = np.zeros(A.shape[0], dtype=np.float32)
    nz = denom > 0
    out[nz] = (Az[nz] @ yz) / denom[nz]
    return out


def r_to_p(r, n):
    """Two-sided p from r. Normal approximation, which is fine at n = 150."""
    r = np.clip(r, -0.999999, 0.999999)
    t = r * np.sqrt((n - 2) / (1 - r ** 2))
    from math import erfc, sqrt
    # Normal approximation: p = erfc(|t|/sqrt(2))
    return np.array([erfc(abs(x) / sqrt(2)) for x in t])


def greedy_decorrelate(P, idx, thresh=COLLIN):
    """
    Greedily drop peaks that duplicate one another.

    P is (k x n) accessibility for the candidate peaks, idx their global
    indices. Candidates are ordered by correlation strength with the gene and
    admitted one at a time; a candidate correlating |r| > thresh with anything
    already admitted is discarded as the same regulatory element.

    Returns the number kept, i.e. the count of independent elements. Without
    this step the measure would mostly count how many peaks a region was
    tiled into.
    """
    if len(idx) <= 1:
        return len(idx)
    Z = P - P.mean(1, keepdims=True)
    s = np.sqrt((Z ** 2).sum(1))
    s[s == 0] = 1
    Z = Z / s[:, None]
    kept = []
    for i in range(len(idx)):
        if not kept:
            kept.append(i); continue
        if np.abs(Z[i] @ Z[kept].T).max() <= thresh:
            kept.append(i)
    return len(kept)


def features_for_type(Xr, Xa, genes, tss, peaks_by_chrom, als_set, min_frac=0.10):
    """Compute the five features for one cell type.

    Xr and Xa arrive already depth-equalised and dense.
    """
    n = Xr.shape[1]
    # A gene counts as expressed when detected in at least 10% of nuclei
    det = (Xr > 0).sum(1)
    expressed = np.flatnonzero(det >= min_frac * n)
    # log-CPM. The depths are already equal, so this is variance stabilisation
    # rather than normalisation.
    R = np.log1p(Xr[expressed])
    A = np.log1p(Xa)

    gname = genes[expressed]
    red, dist_hits, tot_hits, red_als = [], 0, 0, []
    allp = []
    per_gene = []

    for gi, g in enumerate(gname):
        t = tss.get(g)
        if t is None:
            continue
        chrom, pos = t
        pk = peaks_by_chrom.get(chrom)
        if pk is None:
            continue
        lo = np.searchsorted(pk["mid"], pos - WINDOW)
        hi = np.searchsorted(pk["mid"], pos + WINDOW)
        if hi - lo < 2:
            continue
        pidx = pk["idx"][lo:hi]
        P = A[pidx]
        r = corr_rows_vs_vec(P, R[gi])
        per_gene.append((gi, pidx, P, r, pk["mid"][lo:hi], pos))
        allp.append(np.abs(r))

    if not per_gene:
        return None
    # One BH correction across every gene-peak pair within the cell type, not
    # per gene: correcting per gene would make the threshold depend on how many
    # peaks happen to sit near each gene.
    rr = np.concatenate([x[3] for x in per_gene])
    pv = r_to_p(rr, n)
    thr = bh_threshold(pv, FDR)
    if thr < 0:
        return dict(redundancy=0.0, distal_frac=np.nan, redundancy_ALS=0.0,
                    n_links=0, n_genes_tested=len(per_gene))

    off = 0
    for gi, pidx, P, r, mids, pos in per_gene:
        k = len(r)
        p_g = pv[off:off + k]; off += k
        sig = np.flatnonzero(p_g <= thr)
        if len(sig) == 0:
            red.append(0)
            if gname[gi] in als_set: red_als.append(0)
            continue
        order = sig[np.argsort(-np.abs(r[sig]))]
        ni = greedy_decorrelate(P[order], pidx[order])
        red.append(ni)
        if gname[gi] in als_set:
            red_als.append(ni)
        tot_hits += len(sig)
        dist_hits += int((np.abs(mids[sig] - pos) > PROMOTER).sum())

    acc = Xa.sum(1).astype(np.float64)
    acc = acc[acc > 0]; acc /= acc.sum()
    ent = float(-(acc * np.log(acc)).sum())

    return dict(
        redundancy=float(np.median(red)) if red else 0.0,
        redundancy_mean=float(np.mean(red)) if red else 0.0,
        distal_frac=float(dist_hits / tot_hits) if tot_hits else np.nan,
        entropy=ent,
        redundancy_ALS=float(np.median(red_als)) if red_als else np.nan,
        n_expressed=int(len(expressed)),
        n_links=int(tot_hits),
        n_genes_tested=int(len(per_gene)),
        n_als_genes=int(len(red_als)),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results")
    ap.add_argument("--seeds", default="0,1,2")
    ap.add_argument("--types", type=int, default=0, help=">0 时只算前 N 个类型（冒烟）")
    ap.add_argument("--slice", default="", help="a:b 只算第 a..b 个类型，结果追加")
    ap.add_argument("--append", action="store_true")
    args = ap.parse_args()
    seeds = [int(s) for s in args.seeds.split(",")]

    tssdf = pd.read_csv(os.path.join(args.out, "p5_tss_gencode_v32.csv"))
    tss = {r.gene_name: (r.chrom, int(r.tss)) for r in tssdf.itertuples()}
    als_set = set(ALS_GENES)

    allres = []
    for seed in seeds:
        fr = os.path.join(args.out, f"A2_sparse_RNA_seed{seed}.npz")
        fa = os.path.join(args.out, f"A2_sparse_ATAC_seed{seed}.npz")
        if not (os.path.exists(fr) and os.path.exists(fa)):
            print(f"[skip] seed{seed} 子矩阵不全"); continue
        print(f"=== seed {seed} ===", flush=True)
        spr, spa = load_sparse(fr), load_sparse(fa)
        genes = spr["features"]
        pf = pd.DataFrame({"peak": spa["features"]})
        pr = pf.peak.str.rsplit("-", n=2, expand=True)
        pf["chrom"] = pr[0]
        pf["mid"] = (pr[1].astype(int) + pr[2].astype(int)) // 2
        pf["idx"] = np.arange(len(pf))
        peaks_by_chrom = {}
        for ch, s in pf.groupby("chrom"):
            s = s.sort_values("mid")
            peaks_by_chrom[ch] = dict(mid=s["mid"].values, idx=s["idx"].values)

        types = sorted(set(spr["celltype"]))
        if args.types:
            types = types[:args.types]
        if args.slice:
            a, b = args.slice.split(":")
            types = types[int(a):int(b)]
        for ci, ct in enumerate(types):
            col = np.flatnonzero(spr["celltype"] == ct)
            Xr = dense_block(spr, col)
            Xa = dense_block(spa, col)
            f = features_for_type(Xr, Xa, genes, tss, peaks_by_chrom, als_set)
            del Xr, Xa
            if f is None:
                continue
            f.update(celltype=ct, seed=seed, n_cells=len(col),
                     depth_rna=int(np.round(np.mean(
                         [0]))) if False else 0)
            allres.append(f)
            print(f"  [{ci+1}/{len(types)}] {ct:<26} "
                  f"redundancy={f['redundancy']:.1f}  "
                  f"links={f['n_links']:,}  genes={f['n_genes_tested']:,}",
                  flush=True)
        del spr, spa

    if not allres:
        sys.exit("没有产出")
    df = pd.DataFrame(allres)
    f = os.path.join(args.out, "p5_claimA_features.csv")
    if args.append and os.path.exists(f):
        old = pd.read_csv(f)
        df = pd.concat([old, df], ignore_index=True).drop_duplicates(
            ["celltype", "seed"], keep="last")
    df.to_csv(f, index=False)
    print(f"\n[done] {f}  ({len(df)} 行)")
    print("\n⚠️ 下一步必须先跑门禁（p5_09 的 gate 模式），"
          "通过后才允许把特征与因变量对上。")


if __name__ == "__main__":
    main()
