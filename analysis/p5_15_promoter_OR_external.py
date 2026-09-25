#!/usr/bin/env python3
"""
P5 关键检验 · 启动子富集 OR 的阳性对照（cellranger-arc 官方 link 集）
=====================================================================
本项目最重要的一个检验。**它决定论文的上限。**

我们在 150 细胞的 link 集上测到 OR = 0.663 [0.568, 0.768]（显著**耗竭**），
并据此主张"细胞层 peak–gene 推断在可及细胞数下不可行"。
但这个主张有一个致命的未验证前提：

    ⚠️ **这个诊断真的能检出信号吗？**
       如果它在任何 link 集上都给 ~0.66，那 0.663 反映的是我们实现的问题，
       整个 L3 结论作废。

阳性对照：10x 官方 cellranger-arc 在 3,233 个核上产出的 link 集
（human_brain_3k，134,081 个 peak，CC BY 4.0）。
若此处 OR ≫ 1，诊断成立；若 OR ≈ 1，诊断作废。

⚠️⚠️ 口径必须与我们的分析逐字一致，否则就是第 6 次"拿不同口径的量比较"：
  · 检验集 = 所有 (peak, TSS) 在 ±WINDOW 内的配对（WINDOW 从 bedpe 实测距离取）
  · 启动子邻近 = |peak 中点 − TSS| <= 3 kb
  · TSS 用**同一份** GENCODE v32 表（results/p5_tss_gencode_v32.csv）
  · OR 用 Haldane–Anscombe 校正
"""

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

PROMOTER = 3_000


def load_bedpe(path):
    d = pd.read_csv(path, sep="\t", header=None, comment="#",
                    names=["c1", "s1", "e1", "c2", "s2", "e2", "name",
                           "corr", "x", "y", "sig", "dist", "type"])
    pg = d[d.type.isin(["peak-gene", "gene-peak"])].copy()
    nm = pg.name.str.extract(r"^<([^>]*)><([^>]*)>")
    # peak 锚点：peak-gene 时是 A，gene-peak 时是 B
    isA = pg.type.values == "peak-gene"
    pg["peak_mid"] = np.where(isA, (pg.s1 + pg.e1) // 2, (pg.s2 + pg.e2) // 2)
    pg["chrom"] = np.where(isA, pg.c1, pg.c2)
    gname = np.where(isA, nm[1], nm[0])
    # 去掉 10x 的后缀标记（_distal / _promoter / _intergenic）
    pg["gene"] = pd.Series(gname, index=pg.index).str.replace(
        r"_(distal|promoter|intergenic)$", "", regex=True)
    return pg


def build_tested(peaks, tss, window):
    """所有 (peak, TSS) 在 ±window 内的配对数，及其中启动子邻近的数量。"""
    n_tot = n_prox = 0
    for ch, pk in peaks.groupby("chrom"):
        ts = tss[tss.chrom == ch]
        if ts.empty:
            continue
        pm = np.sort(pk.mid.values)
        tv = np.sort(ts.tss.values)
        lo = np.searchsorted(pm, tv - window, "left")
        hi = np.searchsorted(pm, tv + window, "right")
        n_tot += int((hi - lo).sum())
        lo2 = np.searchsorted(pm, tv - PROMOTER, "left")
        hi2 = np.searchsorted(pm, tv + PROMOTER, "right")
        n_prox += int((hi2 - lo2).sum())
    return n_tot, n_prox


def odds_ratio(a, b, c, d, n_boot=2000, seed=0):
    """a/b = 检出中 近端/远端；c/d = 检验集中 近端/远端。Haldane–Anscombe。"""
    or_ = ((a + .5) * (d + .5)) / ((b + .5) * (c + .5))
    rng = np.random.default_rng(seed)
    n1, p1 = a + b, (a + .5) / (a + b + 1)
    bs = []
    for _ in range(n_boot):
        aa = rng.binomial(n1, p1)
        bs.append(((aa + .5) * (d + .5)) / ((n1 - aa + .5) * (c + .5)))
    return or_, np.percentile(bs, [2.5, 97.5])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bedpe", required=True)
    ap.add_argument("--peaks", required=True)
    ap.add_argument("--tss", required=True)
    ap.add_argument("--label", default="external")
    ap.add_argument("--out", default="results")
    args = ap.parse_args()

    tss = pd.read_csv(args.tss)
    pk = pd.read_csv(args.peaks, sep="\t", header=None, comment="#",
                     usecols=[0, 1, 2], names=["chrom", "start", "end"])
    pk = pk[pk.chrom.astype(str).str.startswith("chr")]
    pk["mid"] = (pk.start + pk.end) // 2
    print(f"{len(pk):,} 个 peak / {len(tss):,} 个 TSS")

    pg = load_bedpe(args.bedpe)
    print(f"{len(pg):,} 条 peak–gene link")

    # 窗口从实测距离取（cellranger-arc 默认 ±1 Mb）
    window = int(np.ceil(pg.dist.max() / 1e5) * 1e5)
    print(f"窗口（由 bedpe 实测最大距离定）: ±{window:,}")

    # 检出 link 的启动子邻近判定：用**同一份 TSS**
    t = tss.set_index("gene_name")
    hit = pg.gene.isin(t.index)
    print(f"  基因名可映射到 GENCODE v32: {hit.mean()*100:.1f}% "
          f"({int(hit.sum()):,}/{len(pg):,})")
    if hit.mean() < 0.5:
        sys.exit("[FATAL] 基因名映射率过低，检查 TSS 表与 bedpe 的命名体系")
    q = pg[hit].copy()
    q["tss"] = t.loc[q.gene, "tss"].values
    q["tchrom"] = t.loc[q.gene, "chrom"].values
    q = q[q.chrom == q.tchrom]
    q["d2tss"] = (q.peak_mid - q.tss).abs()

    a = int((q.d2tss <= PROMOTER).sum())      # 检出 & 近端
    b = len(q) - a                            # 检出 & 远端
    n_tot, n_prox = build_tested(pk, tss, window)
    c, d = n_prox, n_tot - n_prox             # 检验集

    print(f"\n检验集：{n_tot:,} 对，其中启动子邻近 {c:,} ({100*c/n_tot:.3f}%)")
    print(f"检出集：{len(q):,} 条，其中启动子邻近 {a:,} ({100*a/len(q):.3f}%)")

    or_, ci = odds_ratio(a, b, c, d)
    print("\n" + "=" * 66)
    print(f"启动子富集 OR = {or_:.3f}   95% CI [{ci[0]:.3f}, {ci[1]:.3f}]")
    print("=" * 66)
    print(f"  本项目（150 细胞，自建管线）: OR = 0.663 [0.568, 0.768]")
    print("\n判定：")
    if ci[0] > 1.5:
        print("  ⭐⭐ **阳性对照成立** —— 诊断能在有信号的 link 集上检出强富集。")
        print("      → 我们的 0.663 反映的是**数据**，不是实现问题。")
        print("      → L3 结论站得住，可以往下做文献普查。")
    elif ci[0] > 1.0:
        print("  ⭐ 显著富集但幅度有限，诊断可用，写作时须给出本对照的具体数值。")
    else:
        print("  ❌ **阳性对照失败** —— 诊断在官方流程的 link 集上也测不到富集。")
        print("     → 0.663 很可能是实现/口径问题，L3 的核心证据作废。")
        print("     → 立即停止 Nature Methods 方案。")

    os.makedirs(args.out, exist_ok=True)
    pd.DataFrame([dict(label=args.label, n_link=len(q), prox_link=a,
                       n_tested=n_tot, prox_tested=c, window=window,
                       OR=or_, lo=ci[0], hi=ci[1])]).to_csv(
        os.path.join(args.out, f"p5_promoterOR_{args.label}.csv"), index=False)


if __name__ == "__main__":
    main()
