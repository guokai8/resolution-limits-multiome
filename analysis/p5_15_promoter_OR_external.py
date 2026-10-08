#!/usr/bin/env python3
"""
Positive control for the promoter-enrichment diagnostic
======================================================
The single most important check in the project, because it sets a ceiling on
what the paper can claim.

At 150 nuclei the link set gives OR = 0.663 [0.568, 0.768], a significant
DEPLETION, and the paper reads that as evidence that nucleus-level peak-gene
inference is not supported at attainable nucleus numbers. That reading rests on
one unverified premise:

    Can this diagnostic detect signal at all?

If it returned ~0.66 on every link set, then 0.663 would be a property of our
implementation rather than of the data, and the whole Layer 3 conclusion would
be void.

The control: the link set 10x's own cellranger-arc produces from 3,233 nuclei
(human_brain_3k, 134,081 peaks, CC BY 4.0). An odds ratio well above 1 there
means the diagnostic works; an odds ratio near 1 means it does not, and the
Layer 3 result would have to be withdrawn.

The conventions must match the main analysis exactly, or this becomes another
comparison between quantities computed different ways:
  * background = every (peak, TSS) pair within +/-WINDOW
  * promoter proximity = |peak midpoint - TSS| <= 3 kb
  * TSS from the SAME GENCODE v32 table (results/p5_tss_gencode_v32.csv)
  * Haldane-Anscombe corrected odds ratio
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
    # The peak anchor is side A for a peak-gene row and side B for a gene-peak
    # row, so which columns to read depends on the row's orientation.
    isA = pg.type.values == "peak-gene"
    pg["peak_mid"] = np.where(isA, (pg.s1 + pg.e1) // 2, (pg.s2 + pg.e2) // 2)
    pg["chrom"] = np.where(isA, pg.c1, pg.c2)
    gname = np.where(isA, nm[1], nm[0])
    # Strip 10x's suffix tags (_distal / _promoter / _intergenic) from the names
    pg["gene"] = pd.Series(gname, index=pg.index).str.replace(
        r"_(distal|promoter|intergenic)$", "", regex=True)
    return pg


def build_tested(peaks, tss, window):
    """Count all (peak, TSS) pairs within +/-window, and how many are proximal.

    This is the background the odds ratio is computed against.
    """
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
    """Haldane-Anscombe corrected odds ratio.

    a/b is proximal/distal among detected links; c/d is the same ratio in the
    background. The +0.5 keeps the ratio finite when a cell is empty.
    """
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
    ap.add_argument("--window", type=int, default=0,
                    help="±窗口（bp）。0 = 沿用旧行为，由 bedpe 实测最大距离定。"
                         "给定时**同时**收紧背景集与检出集，两者必须一致。")
    args = ap.parse_args()

    tss = pd.read_csv(args.tss)
    pk = pd.read_csv(args.peaks, sep="\t", header=None, comment="#",
                     usecols=[0, 1, 2], names=["chrom", "start", "end"])
    pk = pk[pk.chrom.astype(str).str.startswith("chr")]
    pk["mid"] = (pk.start + pk.end) // 2
    print(f"{len(pk):,} 个 peak / {len(tss):,} 个 TSS")

    pg = load_bedpe(args.bedpe)
    print(f"{len(pg):,} 条 peak–gene link")

    # Window: taken from the observed distances by default (cellranger-arc uses
    # +/-1 Mb), or given explicitly.
    #
    # The window defines BOTH the background (build_tested) and the detected
    # set. An earlier version applied it only to the background, which happened
    # to be harmless while the window equalled the largest observed distance --
    # every link was inside it anyway. The moment the window is tightened
    # explicitly, the detected links must be filtered to match, or the 2x2 table
    # is internally inconsistent and the odds ratio is wrong.
    if args.window > 0:
        window = args.window
        print(f"窗口（显式指定）: ±{window:,}")
    else:
        window = int(np.ceil(pg.dist.max() / 1e5) * 1e5)
        print(f"窗口（由 bedpe 实测最大距离定）: ±{window:,}")

    # Proximity of detected links judged against the SAME TSS table
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

    # The detected set must use the same window as the background
    n_before = len(q)
    q = q[q.d2tss <= window]
    if n_before != len(q):
        print(f"  窗口过滤：检出 link {n_before:,} -> {len(q):,} "
              f"（去掉 {n_before - len(q):,} 条超出 ±{window:,} 的）")

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
