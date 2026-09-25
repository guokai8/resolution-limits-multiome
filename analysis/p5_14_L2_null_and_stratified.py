#!/usr/bin/env python3
"""
P5 审查 · L2 复现结果的两项必做检验（动笔前）
================================================
背景：SEA-AD 的 L2 给出 b=−0.505 [−0.524,−0.483]，与 Ulm 的 −0.507 几乎相同。
在把它写成"跨队列可复现的标度律"之前，必须先排除两个使其失效的解释。

检验 1 ⭐ **零模型**：b=−0.5 是不是这个统计量的数学必然？
  做法：把每对文库的伪批量**用同一个共同谱**做多项式重抽样（保持各自实际深度），
        再用完全相同的流程算 floor 与斜率。
  · 若模拟 b ≈ −0.5 且残差量级与实测相当
      → L2 的"复现"只是确认零模型成立，**不是发现**。
        必须改写为"表达层的噪声与纯抽样不可区分"，并明说这是理论预期。
  · 若模拟 b ≈ −0.5 但**截距 a 显著低于实测**
      → 实测比纯抽样更嘈杂，超出部分是真实的技术噪声，a 有信息量。
        → 这才支持"标定 a + 用 b=−0.5 外推"这个协议。

检验 2 ⭐ **类型内分层**：斜率是不是被"细胞类型丰度"这个混杂驱动的？
  L2 的 n_eff 主要随细胞类型丰度变化，而类型本身有不同的生物学与噪声。
  做法：在**每个细胞类型内部**拟合（此时 n 的变异只来自供体/文库），
        看斜率是否仍集中在 −0.5。Ulm 的类型内中位是 −0.530。

用法（几秒，复用已存伪批量，不重跑 1–2 小时的扫描）：
  python3 p5_14_L2_null_and_stratified.py --pb results/seaad_pseudobulk_L2.npz
  python3 p5_14_L2_null_and_stratified.py --pb ... --pairs results/seaad_L2_pairs_genelevel.csv
"""

import argparse
import os

import numpy as np
import pandas as pd


def fit(x, y, groups=None, n_boot=1000, seed=0):
    lx, ly = np.log(np.asarray(x, float)), np.log(np.asarray(y, float))
    A = np.vstack([lx, np.ones_like(lx)]).T
    beta = np.linalg.lstsq(A, ly, rcond=None)[0]
    if groups is None or n_boot == 0:
        return beta, None
    rng = np.random.default_rng(seed)
    uq, bs = np.unique(groups), []
    for _ in range(n_boot):
        pick = rng.choice(uq, len(uq), replace=True)
        idx = np.concatenate([np.flatnonzero(groups == g) for g in pick])
        if len(idx) < 5:
            continue
        Ab = np.vstack([lx[idx], np.ones_like(idx, dtype=float)]).T
        bs.append(np.linalg.lstsq(Ab, ly[idx], rcond=None)[0])
    bs = np.array(bs)
    return beta, np.percentile(bs, [2.5, 97.5], axis=0)


def floor_from_counts(c1, c2, min_cpm=0.0):
    """与 p5_04 / p5_13 逐字一致的下限定义。"""
    t1, t2 = c1.sum(), c2.sum()
    if t1 <= 0 or t2 <= 0:
        return np.nan
    p1, p2 = c1 / t1 * 1e6, c2 / t2 * 1e6
    keep = (p1 > min_cpm) | (p2 > min_cpm)
    if keep.sum() < 200:
        return np.nan
    return float(np.median(np.abs(np.log2(p1[keep] + 1) - np.log2(p2[keep] + 1))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pb", required=True, help="results/seaad_pseudobulk_L2.npz")
    ap.add_argument("--out", default="results")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--min-cpm", type=float, default=0.0)
    args = ap.parse_args()

    z = np.load(args.pb, allow_pickle=True)
    acc = z["acc"]
    meta = pd.DataFrame(z["meta"], columns=list(z["meta_cols"]))
    for c in ("gid", "n_cells"):
        meta[c] = meta[c].astype(int)
    print(f"伪批量 {acc.shape[0]:,} 组 × {acc.shape[1]:,} 基因")

    rng = np.random.default_rng(args.seed)
    obs_rows, sim_rows = [], []
    for (don, ct), s in meta.groupby(["donor", "ct"]):
        if len(s) != 2:
            continue
        (i1, n1), (i2, n2) = [(r.gid, r.n_cells) for r in s.itertuples()]
        if min(n1, n2) < 5:
            continue
        c1, c2 = acc[i1].astype(np.float64), acc[i2].astype(np.float64)
        t1, t2 = c1.sum(), c2.sum()
        if t1 <= 0 or t2 <= 0:
            continue
        ne = 2 / (1 / n1 + 1 / n2)

        f = floor_from_counts(c1, c2, args.min_cpm)
        if np.isfinite(f) and f > 0:
            obs_rows.append(dict(donor=don, ct=ct, floor=f, n_eff=ne))

        # ---- 零模型：共同谱 + 各自实际深度的多项式重抽样 ----
        pool = c1 + c2
        p = pool / pool.sum()
        nz = p > 0
        s1 = np.zeros_like(c1); s2 = np.zeros_like(c2)
        s1[nz] = rng.multinomial(int(t1), p[nz])
        s2[nz] = rng.multinomial(int(t2), p[nz])
        fs = floor_from_counts(s1, s2, args.min_cpm)
        if np.isfinite(fs) and fs > 0:
            sim_rows.append(dict(donor=don, ct=ct, floor=fs, n_eff=ne))

    o, m = pd.DataFrame(obs_rows), pd.DataFrame(sim_rows)
    print(f"实测 {len(o):,} 对 / 零模型 {len(m):,} 对\n")

    # ================================================= 检验 1：零模型
    bo, cio = fit(o.n_eff, o.floor, o.donor.values)
    bm, cim = fit(m.n_eff, m.floor, m.donor.values)
    print("=" * 70)
    print("检验 1 · 零模型（纯多项式抽样，共同谱，实际深度）")
    print("=" * 70)
    print(f"{'':<10}{'斜率 b':>10}{'95% CI':>24}{'系数 a':>10}")
    print(f"{'实测':<10}{bo[0]:>10.3f}"
          f"{f'[{cio[0,0]:.3f}, {cio[1,0]:.3f}]':>24}{np.exp(bo[1]):>10.2f}")
    print(f"{'零模型':<10}{bm[0]:>10.3f}"
          f"{f'[{cim[0,0]:.3f}, {cim[1,0]:.3f}]':>24}{np.exp(bm[1]):>10.2f}")
    ratio = np.exp(bo[1]) / np.exp(bm[1])
    print(f"\n⭐ 实测/零模型 的系数比 = {ratio:.2f}×")
    # 逐对比值更稳健
    j = o.merge(m, on=["donor", "ct"], suffixes=("_o", "_m"))
    pr = (j.floor_o / j.floor_m)
    print(f"   逐对 floor 比值：中位 {pr.median():.2f}×  "
          f"IQR [{pr.quantile(.25):.2f}, {pr.quantile(.75):.2f}]")
    print("\n判读：")
    if pr.median() < 1.15:
        print("  ⚠️ **实测≈零模型** → L2 的 b=−0.5 是统计量的数学必然，")
        print("     '跨队列复现'只是确认零模型成立，不是发现。")
        print("     → 必须改写：L2 是**阴性对照臂**，其价值全在于反衬 L1 的偏离；")
        print("        不能单独作为卖点，摘要里也不应写成'标度律'。")
    else:
        print(f"  ✅ **实测比纯抽样嘈杂 {pr.median():.2f}×** → 超出部分是真实技术噪声，")
        print("     系数 a 有信息量 → 支持'一对重复标定 a + b=−0.5 外推'的协议。")
        print("     但 b 本身仍是理论预期，写作时必须明说，不可包装成新发现。")

    # ============================================= 检验 2：类型内分层
    print("\n" + "=" * 70)
    print("检验 2 · 细胞类型内分层（排除'类型丰度'混杂）")
    print("=" * 70)
    rows = []
    for ct, s in o.groupby("ct"):
        if len(s) < 8 or s.n_eff.max() / s.n_eff.min() < 3:
            continue
        b, _ = fit(s.n_eff, s.floor, None, 0)
        rows.append(dict(ct=ct, n=len(s), span=s.n_eff.max() / s.n_eff.min(),
                         b=b[0]))
    st = pd.DataFrame(rows).sort_values("b")
    if st.empty:
        print("  ⚠️ 无类型满足 (n≥8 且 n_eff 跨度≥3×) —— 无法分层检验")
    else:
        print(f"{'细胞类型':<24}{'n':>5}{'跨度':>8}{'b':>9}")
        for r in st.itertuples():
            print(f"{str(r.ct)[:23]:<24}{r.n:>5}{r.span:>8.1f}{r.b:>9.3f}")
        print(f"\n⭐ 类型内斜率：中位 {st.b.median():.3f}  "
              f"IQR [{st.b.quantile(.25):.3f}, {st.b.quantile(.75):.3f}]  "
              f"({len(st)} 个类型)     （Ulm 类型内中位 −0.530）")
        near = ((st.b > -0.7) & (st.b < -0.3)).mean()
        print(f"   落在 (−0.7, −0.3) 的比例：{near:.0%}")
        if abs(st.b.median() + 0.5) < 0.12:
            print("   ✅ 类型内仍集中于 −0.5 → 斜率**不是**类型丰度混杂造成的")
        else:
            print("   ⚠️ 类型内中位偏离 −0.5 → 合并斜率主要由跨类型变异驱动，")
            print("      两层对照的'同一回归'表述必须收回")

    os.makedirs(args.out, exist_ok=True)
    o.to_csv(os.path.join(args.out, "seaad_L2_obs.csv"), index=False)
    m.to_csv(os.path.join(args.out, "seaad_L2_null.csv"), index=False)
    if not st.empty:
        st.to_csv(os.path.join(args.out, "seaad_L2_by_celltype.csv"), index=False)
    print(f"\n[done] {args.out}/")


if __name__ == "__main__":
    main()
