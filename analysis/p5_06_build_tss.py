#!/usr/bin/env python3
"""
P5 步骤 0.7 后处理 · 从 GENCODE GTF 建 TSS 表，并做 peak↔gene 配对预检
=====================================================================
在 `p5_05_download_refs.sh` 下完 GTF 之后跑这个。它做三件事：

  1. 解析 GTF → 每个基因一行的 TSS 表（按链取正确的 5' 端）
  2. **硬校验**：TSS 表能覆盖多少 Multiome 的 35,367 个基因名？
     覆盖率过低说明版本或命名空间不对，必须停下而不是将就往下走
  3. 预统计每个基因 ±500 kb 内有多少个 ATAC peak
     —— 这是 A2 `redundancy` 的**分母上界**，也是判断该指标有没有动态范围的前提

为什么必须做第 3 步：如果几乎所有基因附近的 peak 数都差不多，
`redundancy` 就没有类型间变异可言，命题 A 在技术上就不成立 —— 这一点
**在算任何相关性之前**就能知道，属于预注册 §4.1 之前的可行性检查。

用法：
  python3 p5_06_build_tss.py --ref refs --data ~/Desktop/ResearchD --out results
"""

import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd

# Ulm multiome 用 10x GRCh38-2020-A 参考 = GENCODE v32 / Ensembl 98
EXPECTED = "gencode.v32.annotation.gtf.gz"
MIN_COVERAGE = 0.85          # 低于此覆盖率即中止


def parse_gtf_tss(path):
    """解析 GTF 的 gene 行 → DataFrame[gene_name, gene_id, chrom, tss, strand, gene_type]"""
    rows = []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            chrom, start, end, strand, attr = f[0], int(f[3]), int(f[4]), f[6], f[8]
            # TSS = 按链取 5' 端。写反了会让所有 promoter/distal 划分失效。
            tss = start if strand == "+" else end
            d = {}
            for kv in attr.strip().split(";"):
                kv = kv.strip()
                if not kv:
                    continue
                k, _, v = kv.partition(" ")
                d[k] = v.strip('"')
            rows.append((d.get("gene_name"), d.get("gene_id"), chrom, tss,
                         strand, d.get("gene_type")))
    g = pd.DataFrame(rows, columns=["gene_name", "gene_id", "chrom", "tss",
                                    "strand", "gene_type"])
    return g.dropna(subset=["gene_name"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", default="refs")
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", default="results")
    ap.add_argument("--window", type=int, default=500_000,
                    help="peak↔gene 配对窗口（bp，单侧）")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    gtf = os.path.join(args.ref, EXPECTED)
    if not os.path.exists(gtf):
        sys.exit(f"[FATAL] 找不到 {gtf}\n先跑： bash code/p5_05_download_refs.sh")

    print("[1/4] 解析 GTF …")
    g = parse_gtf_tss(gtf)
    print(f"  {len(g):,} 个基因；gene_type 前 5: "
          f"{g.gene_type.value_counts().head(5).to_dict()}")

    print("[2/4] 与 Multiome 基因名对齐（硬校验）…")
    feats = pd.read_csv(os.path.join(args.data,
                        "Multiome_Dataset_RNA_features.tsv"), header=None)[0]
    # GTF 里同名基因可能多条（PAR 区、多版本）；每个名字保留一条主记录
    g1 = (g.sort_values(["gene_name", "chrom", "tss"])
            .drop_duplicates("gene_name", keep="first")
            .set_index("gene_name"))
    hit = feats.isin(g1.index)
    cov = hit.mean()
    print(f"  Multiome 基因 {len(feats):,}，匹配上 {int(hit.sum()):,} "
          f"（覆盖率 {cov:.1%}）")
    if cov < MIN_COVERAGE:
        sys.exit(f"[FATAL] 覆盖率 {cov:.1%} < {MIN_COVERAGE:.0%}。"
                 "很可能 GENCODE 版本不对（应为 v32）或命名空间不匹配。停止。")
    miss = feats[~hit]
    print(f"  未匹配示例: {list(miss[:8])}")

    tss = g1.loc[feats[hit]].reset_index()[
        ["gene_name", "gene_id", "chrom", "tss", "strand", "gene_type"]]
    f_tss = os.path.join(args.out, "p5_tss_gencode_v32.csv")
    tss.to_csv(f_tss, index=False)

    print("[3/4] 读 ATAC peak 坐标 …")
    pk = pd.read_csv(os.path.join(args.data,
                     "Multiome_Dataset_ATAC_features.tsv"), header=None)[0]
    # 形如 chr1-9908-10537
    parts = pk.str.rsplit("-", n=2, expand=True)
    peaks = pd.DataFrame({"peak": pk, "chrom": parts[0],
                          "start": parts[1].astype(int),
                          "end": parts[2].astype(int)})
    peaks["mid"] = (peaks.start + peaks.end) // 2
    print(f"  {len(peaks):,} 个 peak；染色体数 {peaks.chrom.nunique()}")
    bad = set(peaks.chrom) - set(tss.chrom)
    if bad:
        print(f"  ⚠️ peak 中有 {len(bad)} 个染色体名不在 GTF 里（将被忽略）: "
              f"{sorted(bad)[:5]}")

    print(f"[4/4] 统计每个基因 ±{args.window//1000} kb 内的 peak 数 …")
    counts, prom = [], []
    for chrom, sub in tss.groupby("chrom"):
        pm = np.sort(peaks.loc[peaks.chrom == chrom, "mid"].values)
        if len(pm) == 0:
            counts += [0] * len(sub); prom += [0] * len(sub); continue
        lo = np.searchsorted(pm, sub.tss.values - args.window)
        hi = np.searchsorted(pm, sub.tss.values + args.window)
        counts += list(hi - lo)
        plo = np.searchsorted(pm, sub.tss.values - 3000)
        phi = np.searchsorted(pm, sub.tss.values + 3000)
        prom += list(phi - plo)
    order = tss.groupby("chrom").cumcount().index  # 保持 groupby 顺序
    tss2 = pd.concat([sub for _, sub in tss.groupby("chrom")])
    tss2["n_peaks_window"] = counts
    tss2["n_peaks_promoter"] = prom
    tss2 = tss2.sort_values("gene_name")
    tss2.to_csv(os.path.join(args.out, "p5_gene_peak_window_counts.csv"),
                index=False)

    n = tss2["n_peaks_window"]
    print(f"\n  每基因窗口内 peak 数： 中位 {n.median():.0f}  "
          f"IQR [{n.quantile(.25):.0f}, {n.quantile(.75):.0f}]  "
          f"范围 {n.min()}–{n.max()}")
    print(f"  有启动子 peak 的基因: {(tss2.n_peaks_promoter>0).mean():.1%}")
    print(f"  窗口内 0 个 peak 的基因: {(n==0).mean():.1%}")

    print("\n" + "=" * 72)
    print("可行性判定（在算任何相关性之前）")
    print("=" * 72)
    cv = n.std() / max(n.mean(), 1e-9)
    print(f"  每基因 peak 数的变异系数 CV = {cv:.2f}")
    if cv < 0.3:
        print("  ❌ 变异太小 —— redundancy 几乎没有动态范围，命题 A 在技术上难以成立。")
    else:
        print("  ✅ 有足够动态范围，可以进入 A2。")
    print("\n  ⚠️ 注意：这里数的是**窗口内 peak 总数**（分母上界），")
    print("     不是 A2 要的 `redundancy`（显著相关且去共线后的 link 数）。")
    print("     窗口内 peak 数强烈依赖基因密度，**必须**按预注册 §4.1 做")
    print("     深度降采样，并检验 partial_r(redundancy, depth) < 0.2。")
    print(f"\n[done] {f_tss}")
    print(f"[done] {os.path.join(args.out,'p5_gene_peak_window_counts.csv')}")


if __name__ == "__main__":
    main()
