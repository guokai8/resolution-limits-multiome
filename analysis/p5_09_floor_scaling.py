#!/usr/bin/env python3
"""
Scaling of the noise floor: can the floor be predicted from nucleus count?
=========================================================================
This is the step that turns "we measured a floor in this dataset" into "the
floor at a given nucleus count is predictable" -- the difference between one
data point and a ruler someone else can use.

The question: how does the technical noise floor vary with the number of nuclei
in a cell type?

    floor(n) = a * n^b

  b ~ -0.5      pure multinomial sampling noise: collecting more nuclei fixes it
  b ~  0        pure batch effect: more nuclei does NOT help, the design must change
  -0.5 < b < 0  a mixture, with an irreducible batch component

If the third case holds, that matters more than the value of the floor itself:
it means piling on nuclei cannot make within-cell-type expression conclusions
reliable, however many you collect.

Method: the unit of observation is (donor x chip pair x cell type). Effective
nucleus count n_eff is the HARMONIC mean of the two chips' counts, because the
variance of a paired difference goes as 1/n1 + 1/n2. Fit log(floor) on
log(n_eff), with the slope interval from a bootstrap clustered BY DONOR, since
the several cell types of one donor are not independent.

Usage:
  python3 p5_09_floor_scaling.py --data DIR --out results
"""

import argparse
import itertools
import os

import numpy as np
import pandas as pd


def load_pairs(out_dir):
    f = os.path.join(out_dir, "p5_techrep_expression_pairs.csv")
    if not os.path.exists(f):
        raise SystemExit(f"[FATAL] 缺 {f}，先跑 p5_04 的 L2 部分")
    return pd.read_csv(f)


def attach_cell_counts(pairs, data_dir, level="WNN_L2.5"):
    m = pd.read_csv(os.path.join(data_dir, "Multiome_Dataset_Metadata.txt"),
                    sep="\t", usecols=["ID", "Well", level], dtype=str)
    m["Chip"] = m["Well"].str.extract(r"^(Chip\d+)")
    n = (m.groupby(["ID", "Chip", level]).size()
           .rename("n").reset_index()
           .rename(columns={"ID": "donor", level: "celltype"}))
    p = pairs.merge(n.rename(columns={"Chip": "chip1", "n": "n1"}),
                    on=["donor", "chip1", "celltype"], how="left")
    p = p.merge(n.rename(columns={"Chip": "chip2", "n": "n2"}),
                on=["donor", "chip2", "celltype"], how="left")
    p = p.dropna(subset=["n1", "n2"])
    # Effective nucleus count is the harmonic mean, because the variance of the
    # paired difference goes as 1/n1 + 1/n2 -- the smaller library dominates.
    p["n_eff"] = 2.0 / (1.0 / p.n1 + 1.0 / p.n2)
    return p


def fit_loglog(x, y, groups, n_boot=2000, seed=0):
    """OLS in log-log space, bootstrapped over `groups` (donors).

    Returns slope and intercept with confidence intervals. Resampling donors
    rather than observations is what keeps the interval honest.
    """
    lx, ly = np.log(x), np.log(y)
    A = np.vstack([lx, np.ones_like(lx)]).T
    beta = np.linalg.lstsq(A, ly, rcond=None)[0]
    rng = np.random.default_rng(seed)
    uniq = np.unique(groups)
    bs = []
    for _ in range(n_boot):
        pick = rng.choice(uniq, len(uniq), replace=True)
        idx = np.concatenate([np.flatnonzero(groups == g) for g in pick])
        if len(idx) < 5:
            continue
        Ab = np.vstack([lx[idx], np.ones_like(idx, dtype=float)]).T
        try:
            bs.append(np.linalg.lstsq(Ab, ly[idx], rcond=None)[0])
        except np.linalg.LinAlgError:
            continue
    bs = np.array(bs)
    return beta, np.percentile(bs, [2.5, 97.5], axis=0), bs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", default="results")
    ap.add_argument("--n-boot", type=int, default=2000)
    args = ap.parse_args()

    p = attach_cell_counts(load_pairs(args.out), args.data)
    p = p[(p.median_abs_log2FC > 0) & (p.n_eff > 0)]
    print(f"观测数 {len(p)}（(供体 × 芯片对 × 细胞类型)），"
          f"{p.donor.nunique()} 个供体，{p.celltype.nunique()} 个细胞类型")
    print(f"n_eff 范围 {p.n_eff.min():.0f}–{p.n_eff.max():.0f}（中位 {p.n_eff.median():.0f}）")
    print(f"floor 范围 {p.median_abs_log2FC.min():.3f}–{p.median_abs_log2FC.max():.3f}")

    beta, ci, bs = fit_loglog(p.n_eff.values, p.median_abs_log2FC.values,
                              p.donor.values, args.n_boot)
    b, loga = beta
    print("\n" + "=" * 66)
    print("标度律  floor(n) = a · n^b")
    print("=" * 66)
    print(f"  斜率 b      = {b:+.3f}   95% CI [{ci[0,0]:+.3f}, {ci[1,0]:+.3f}]")
    print(f"  系数 a      = {np.exp(loga):.2f}")
    print(f"  纯抽样噪声预测 b = −0.500")
    print(f"  纯批次效应预测 b =  0.000")
    p_half = float((bs[:, 0] <= -0.5).mean())
    p_zero = float((bs[:, 0] >= 0.0).mean())
    print(f"\n  b ≤ −0.5 的 bootstrap 比例 = {p_half:.3f}"
          f"  → {'不能' if p_half>0.025 else '可以'}排除纯抽样")
    print(f"  b ≥  0.0 的 bootstrap 比例 = {p_zero:.3f}"
          f"  → {'不能' if p_zero>0.025 else '可以'}排除纯批次")

    if -0.5 < ci[0, 0] and ci[1, 0] < 0:
        frac = 1 - abs(b) / 0.5
        print(f"\n  ⭐ b 显著介于 −0.5 与 0 之间 → 存在**不可约的批次分量**")
        print(f"     约 {100*frac:.0f}% 的噪声无法通过增加细胞数消除")

    print("\n  外推：达到指定下限所需的每类型细胞数")
    print(f"  {'目标 floor':>12}{'所需 n_eff':>14}")
    for tgt in (0.5, 0.3, 0.2, 0.1):
        need = np.exp((np.log(tgt) - loga) / b)
        print(f"  {tgt:>12.2f}{need:>14,.0f}")

    # Observed against predicted, per cell type, as a fit diagnostic
    p["pred"] = np.exp(loga) * p.n_eff ** b
    s = p.groupby("celltype").agg(n_pairs=("n_eff", "size"),
                                  n_eff=("n_eff", "median"),
                                  obs=("median_abs_log2FC", "median"),
                                  pred=("pred", "median")).reset_index()
    s["ratio"] = s.obs / s.pred
    s = s.sort_values("n_eff", ascending=False)
    print("\n  逐细胞类型：观测下限 vs 标度律预测")
    print(f"  {'细胞类型':<20}{'配对':>5}{'n_eff中位':>10}{'实测':>8}{'预测':>8}{'实测/预测':>10}")
    for _, r in s.iterrows():
        flag = "  ← 显著偏高" if r.ratio > 1.5 else ""
        print(f"  {r.celltype:<20}{int(r.n_pairs):>5}{r.n_eff:>10.0f}"
              f"{r.obs:>8.3f}{r.pred:>8.3f}{r.ratio:>10.2f}{flag}")

    p.to_csv(os.path.join(args.out, "p5_floor_scaling_pairs.csv"), index=False)
    s.to_csv(os.path.join(args.out, "p5_floor_scaling_bytype.csv"), index=False)
    pd.DataFrame(dict(b=bs[:, 0], loga=bs[:, 1])).to_csv(
        os.path.join(args.out, "p5_floor_scaling_boot.csv"), index=False)
    print(f"\n[done] {args.out}/p5_floor_scaling_*.csv")


if __name__ == "__main__":
    main()
