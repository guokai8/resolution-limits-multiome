#!/usr/bin/env python3
"""
P5 · 两层噪声结构的严格对照：同一个统计量，同一个回归，两个指数
================================================================
⚠️ 这个脚本是在**检验一个已经写进提纲的论断**，不是补一张表。

提纲里写了"L1 组成层是批次驱动的，加细胞无效"。
但那是从"过度离散 2.36×"推出来的，**而过度离散本身不能区分两种情况**：

  (a) **乘性**放大：噪声 = 2.36 × 抽样噪声
      → 斜率仍为 −0.5，**加细胞仍然有效**（只是需要 2.36² ≈ 5.6 倍）
  (b) **加性**批次：噪声 = 抽样噪声 + 与 n 无关的常数项
      → 斜率被压平（>−0.5），**加细胞无效**

**两者的实践含义完全相反。** 不检验就断言，是逻辑漏洞。

做法：把两层化成**同一个可比的量**。

  L2 表达层： floor_expr(n) = median |Δ log2CPM|
  L1 组成层： 标准化的组成差异
              z = |f1 − f2| / sqrt( p(1−p) )
              纯多项抽样下 SD(f1−f2) = sqrt(p(1−p)(1/n1+1/n2))
              ⇒ z ∝ n_eff^(−0.5)，**与 L2 同样的 −0.5 预测**

于是两层都拟合 log(量) ~ log(n_eff)，斜率直接可比。

**一个佐证性的先验**：`FINDING_P5_04` 观察到过度离散**随丰度上升**
（少突 7.44×，稀有抑制性类型仅 1.3×）。在乘性模型下这不该发生；
在加性模型下，丰度高的类型抽样噪声小、固定批次项占比更大 → 比值更高。
所以先验倾向 (b)，但**必须实测**。

用法：
  python3 p5_10_two_layer_scaling.py --data DIR --out results
"""

import argparse
import itertools
import os

import numpy as np
import pandas as pd


def fit_loglog(x, y, groups, n_boot=2000, seed=0):
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


def composition_pairs(data_dir, level="WNN_L2.5", min_cells=50):
    m = pd.read_csv(os.path.join(data_dir, "Multiome_Dataset_Metadata.txt"),
                    sep="\t", usecols=["ID", "Well", level], dtype=str)
    m["Chip"] = m["Well"].str.extract(r"^(Chip\d+)")
    n = m.groupby(["ID", "Chip", level]).size().rename("n").reset_index()
    tot = n.groupby(["ID", "Chip"])["n"].sum().rename("tot").reset_index()
    n = n.merge(tot, on=["ID", "Chip"])
    n = n[n.tot >= min_cells]
    n["f"] = n.n / n.tot
    rows = []
    for did, sub in n.groupby("ID"):
        chips = sorted(sub.Chip.unique())
        if len(chips) < 2:
            continue
        for c1, c2 in itertools.combinations(chips, 2):
            a = sub[sub.Chip == c1].set_index(level)
            b = sub[sub.Chip == c2].set_index(level)
            for ct in a.index.union(b.index):
                f1 = a.f.get(ct, 0.0); f2 = b.f.get(ct, 0.0)
                n1 = a.tot.get(ct, np.nan); n2 = b.tot.get(ct, np.nan)
                if np.isnan(n1) or np.isnan(n2):
                    continue
                p = (f1 * n1 + f2 * n2) / (n1 + n2)
                if not (0 < p < 1):
                    continue
                # n_eff 用**总核数**的调和平均（组成的抽样单位是全部核）
                neff = 2.0 / (1.0 / n1 + 1.0 / n2)
                z = abs(f1 - f2) / np.sqrt(p * (1 - p))
                if z <= 0:
                    continue
                rows.append(dict(donor=did, celltype=ct, chip1=c1, chip2=c2,
                                 p=p, n_eff=neff, z=z))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", default="results")
    ap.add_argument("--n-boot", type=int, default=1500)
    args = ap.parse_args()

    # ---------------- L1 组成层 ----------------
    c = composition_pairs(args.data)
    print(f"L1 组成层：{len(c)} 个观测，{c.donor.nunique()} 供体，"
          f"{c.celltype.nunique()} 类型")
    print(f"  n_eff（总核数）范围 {c.n_eff.min():.0f}–{c.n_eff.max():.0f}")
    b1, ci1, bs1 = fit_loglog(c.n_eff.values, c.z.values, c.donor.values,
                              args.n_boot)

    # ---------------- L2 表达层（已有） ----------------
    e = pd.read_csv(os.path.join(args.out, "p5_floor_scaling_pairs.csv"))
    b2, ci2, bs2 = fit_loglog(e.n_eff.values, e.median_abs_log2FC.values,
                              e.donor.values, args.n_boot)

    print("\n" + "=" * 70)
    print("两层对照：同一统计量、同一回归、两个指数")
    print("=" * 70)
    print(f"{'层':<16}{'斜率 b':>10}{'95% CI':>22}{'纯抽样预测':>12}")
    print(f"{'L1 组成':<16}{b1[0]:>10.3f}{f'[{ci1[0,0]:.3f}, {ci1[1,0]:.3f}]':>22}{-0.5:>12.3f}")
    print(f"{'L2 表达':<16}{b2[0]:>10.3f}{f'[{ci2[0,0]:.3f}, {ci2[1,0]:.3f}]':>22}{-0.5:>12.3f}")

    # 两个斜率是否不同？按供体配对 bootstrap 求差
    k = min(len(bs1), len(bs2))
    d = bs1[:k, 0] - bs2[:k, 0]
    print(f"\n斜率之差 (L1 − L2) = {b1[0]-b2[0]:+.3f}   "
          f"95% CI [{np.percentile(d,2.5):+.3f}, {np.percentile(d,97.5):+.3f}]")
    sep = np.percentile(d, 2.5) > 0
    print(f"  两层斜率是否显著不同: {'✅ 是' if sep else '❌ 否'}")

    p1_half = float((bs1[:, 0] <= -0.5).mean())
    print(f"\nL1: b ≤ −0.5 的 bootstrap 比例 = {p1_half:.3f}")
    if p1_half < 0.025:
        print("  ⭐ L1 斜率**显著平于** −0.5 → **加性批次分量存在，加细胞无法消除**")
        print(f"     不可约分量占比估计：1 − |b|/0.5 = {1-abs(b1[0])/0.5:.1%}")
    elif p1_half > 0.975:
        print("  L1 斜率显著陡于 −0.5 —— 异常，需追查")
    else:
        print("  ⚠️ L1 与纯抽样**不可区分** → 提纲里"
              "'加细胞无效'的论断**必须撤回**，改为'乘性放大 2.36×'")

    c.to_csv(os.path.join(args.out, "p5_L1_scaling_pairs.csv"), index=False)
    pd.DataFrame(dict(b_L1=bs1[:k, 0], b_L2=bs2[:k, 0])).to_csv(
        os.path.join(args.out, "p5_two_layer_boot.csv"), index=False)
    print(f"\n[done] {args.out}/p5_L1_scaling_pairs.csv, p5_two_layer_boot.csv")


if __name__ == "__main__":
    main()
