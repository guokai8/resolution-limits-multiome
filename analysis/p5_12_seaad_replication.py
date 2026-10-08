#!/usr/bin/env python3
"""
External validation in SEA-AD (84 donors, 1.18M nuclei, middle temporal gyrus)
=============================================================================
Two independent jobs sharing one 1 GB metadata CSV. Neither needs the 31 GB
h5ad except where stated.

  --mode l1        composition-layer replication; metadata CSV only
  --mode scaling   composition and expression layers together. The expression
                   layer needs --h5ad; see the warning on mode_scaling below.
                   80 donors have two or more library preps with at least 50
                   nuclei each. Expectation: the composition exponent is
                   clearly flatter than -0.5, the expression exponent is
                   indistinguishable from -0.5, and their difference is
                   clearly above zero.

  --mode l3sat     the Layer 3 saturation curve. This is NOT a technical
                   replicate analysis. The 28 multiome libraries are one per
                   donor, so there are no replicates; what they give is link
                   count against nucleus count. Sequencing depth has to be
                   controlled first, or an exponent near 1 may be nothing but
                   a depth effect.

Why SEA-AD is a good validation set:
  * 80 replicated donors against the primary cohort's 26
  * median nucleus-count ratio between a donor's two libraries is 1.21, against
    2.67 in the primary cohort, so these are genuinely independent libraries
    rather than one library split
  * 24 subclasses against 12

Two differences from the primary cohort that must be stated in the Methods:
  1. SEA-AD is mainly snRNA-seq (10x v3.1); the primary cohort is multiome. The
     composition and expression definitions do not involve ATAC so they remain
     comparable, but the platform difference belongs in the limitations.
  2. The batch unit differs: chip here (one chip, several wells) against
     library prep there (one preparation). Both are a single independent
     loading and library construction, so they are the same level of technical
     unit.

Usage:
  python3 p5_12_seaad_replication.py --meta SEAAD_..._metadata.csv --mode scaling
  python3 p5_12_seaad_replication.py --meta SEAAD_..._metadata.csv --mode l3sat
"""

import argparse
import itertools
import os

import numpy as np
import pandas as pd

DONOR, BATCH, CTYPE = "Donor ID", "library_prep", "Subclass"


# ---------------------------------------------------------------- shared
def fit_loglog(x, y, groups, n_boot=1500, seed=0):
    """log-log OLS with a bootstrap clustered BY DONOR, since the several cell
    types of one donor are not independent."""
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


def partial_corr(x, y, z):
    """Partial correlation of x and y controlling for z, all on the log scale."""
    lx, ly, lz = np.log(x), np.log(y), np.log(z)
    def resid(v):
        A = np.vstack([lz, np.ones_like(lz)]).T
        return v - A @ np.linalg.lstsq(A, v, rcond=None)[0]
    rx, ry = resid(lx), resid(ly)
    return float(np.corrcoef(rx, ry)[0, 1])


# ---------------------------------------------------------------- main modes
def mode_scaling(d, n_boot, allow_umi_proxy=False):
    """
    Composition and expression scaling, side by side.

    A limitation caught by the smoke test, recorded so it is not repeated.

    The expression floor is the median |delta log2 CPM| ACROSS GENES within a
    cell type, which requires an expression matrix. The SEA-AD metadata CSV
    carries no gene-level expression at all -- only total UMI per nucleus.
    Using |delta log2(median UMI)| as a stand-in gives a systematically
    shallower slope: checked against the primary cohort's own data, the true
    value is -0.507 and the UMI proxy returns -0.248. That is more than enough
    to call a successful replication a failure.

    So the expression layer must come from the 31 GB h5ad, or only the
    composition layer should be reported. This function refuses to print the
    expression-layer comparison unless --allow-umi-proxy is given explicitly.
    """
    d = d[d["Number of UMIs"].notna()].copy()
    d["Number of UMIs"] = pd.to_numeric(d["Number of UMIs"], errors="coerce")
    d = d.dropna(subset=["Number of UMIs", DONOR, BATCH, CTYPE])

    # Total nuclei per (donor, batch), and nuclei plus expression proxy per
    # (donor, batch, cell type)
    tot = d.groupby([DONOR, BATCH]).size().rename("tot")
    per = d.groupby([DONOR, BATCH, CTYPE]).agg(
        n=("Number of UMIs", "size"),
        umi=("Number of UMIs", "median")).reset_index()
    per = per.merge(tot, on=[DONOR, BATCH])
    per["f"] = per.n / per.tot

    rows1, rows2 = [], []
    for don, sub in per.groupby(DONOR):
        bats = sorted(sub[BATCH].unique())
        if len(bats) < 2:
            continue
        for b1, b2 in itertools.combinations(bats, 2):
            a = sub[sub[BATCH] == b1].set_index(CTYPE)
            b = sub[sub[BATCH] == b2].set_index(CTYPE)
            if a.empty or b.empty:
                continue
            T1, T2 = a.tot.iloc[0], b.tot.iloc[0]
            if min(T1, T2) < 50:
                continue
            # ---- Composition: standardised proportion difference. n_eff uses
            # the harmonic mean of TOTAL nuclei, since a proportion's
            # denominator is the whole library.
            for ct in a.index.union(b.index):
                f1 = a.f.get(ct, 0.0); f2 = b.f.get(ct, 0.0)
                p = (f1 * T1 + f2 * T2) / (T1 + T2)
                if not (0 < p < 1):
                    continue
                z = abs(f1 - f2) / np.sqrt(p * (1 - p))
                if z > 0:
                    rows1.append(dict(donor=don, ct=ct, z=z,
                                      n_eff=2 / (1 / T1 + 1 / T2)))
            # ---- Expression: |delta log2(median UMI)|, n_eff from this cell
            # type's own nuclei. This is the UMI proxy; see the docstring.
            for ct in a.index.intersection(b.index):
                n1, n2 = a.n.get(ct, 0), b.n.get(ct, 0)
                if min(n1, n2) < 5:
                    continue
                u1, u2 = a.umi.get(ct), b.umi.get(ct)
                fl = abs(np.log2(max(u1, 1) / max(u2, 1)))
                if fl > 0:
                    rows2.append(dict(donor=don, ct=ct, floor=fl,
                                      n_eff=2 / (1 / n1 + 1 / n2)))

    c, e = pd.DataFrame(rows1), pd.DataFrame(rows2)
    print(f"L1 组成：{len(c):,} 观测 / {c.donor.nunique()} 供体 / {c.ct.nunique()} 类型")
    print(f"L2 表达：{len(e):,} 观测 / {e.donor.nunique()} 供体 / {e.ct.nunique()} 类型")
    print(f"  n_eff 跨度  L1 {c.n_eff.min():.0f}–{c.n_eff.max():.0f}"
          f"（{c.n_eff.max()/c.n_eff.min():.0f}×）"
          f"  L2 {e.n_eff.min():.0f}–{e.n_eff.max():.0f}"
          f"（{e.n_eff.max()/e.n_eff.min():.0f}×）")

    b1, ci1, bs1 = fit_loglog(c.n_eff.values, c.z.values, c.donor.values, n_boot)
    print(f"\n{'='*68}\nL1 组成层（元数据即可，与 Ulm 口径一致）\n{'='*68}")
    print(f"  b = {b1[0]:+.3f}  95% CI [{ci1[0,0]:+.3f}, {ci1[1,0]:+.3f}]"
          f"   （Ulm: −0.152 [−0.340, +0.111]）")
    print(f"  排除纯抽样(b≤−0.5 的 bootstrap 比例) = "
          f"{float((bs1[:,0]<=-0.5).mean()):.3f}"
          f"{'  ✅ 复现' if float((bs1[:,0]<=-0.5).mean())<0.025 else '  ❌ 未复现'}")
    if not allow_umi_proxy:
        print("\n⚠️ L2 已跳过：元数据只有每细胞总 UMI，无基因层表达。")
        print("   用 UMI 代理会系统性低估斜率（Ulm 自测：真值 −0.507 → 代理 −0.248），")
        print("   足以把'复现'误判为'未复现'。")
        print("   → 要复现 L2，需下载 31 GB 的 h5ad 并用基因层 log-CPM 重算。")
        return c, e
    print("\n⚠️⚠️ 以下 L2 使用 UMI 代理，**不可与 Ulm 的 −0.507 直接比较**")
    b2, ci2, bs2 = fit_loglog(e.n_eff.values, e.floor.values, e.donor.values, n_boot)

    print("\n" + "=" * 68)
    print("SEA-AD 独立复现：两层标度对照")
    print("=" * 68)
    print(f"{'层':<12}{'斜率 b':>10}{'95% CI':>22}{'Ulm 原值':>12}")
    print(f"{'L1 组成':<12}{b1[0]:>10.3f}"
          f"{f'[{ci1[0,0]:.3f}, {ci1[1,0]:.3f}]':>22}{-0.152:>12.3f}")
    print(f"{'L2 表达':<12}{b2[0]:>10.3f}"
          f"{f'[{ci2[0,0]:.3f}, {ci2[1,0]:.3f}]':>22}{-0.507:>12.3f}")

    k = min(len(bs1), len(bs2))
    dif = bs1[:k, 0] - bs2[:k, 0]
    lo, hi = np.percentile(dif, [2.5, 97.5])
    print(f"\n斜率之差 (L1−L2) = {b1[0]-b2[0]:+.3f}  95% CI [{lo:+.3f}, {hi:+.3f}]"
          f"   （Ulm: +0.355 [+0.158, +0.633]）")
    print(f"  两层显著不同: {'✅ 是' if lo > 0 else '❌ 否'}")
    print(f"  L1 排除纯抽样(b≤−0.5 比例): {float((bs1[:,0]<=-0.5).mean()):.3f}")
    print(f"  L2 与纯抽样一致: {'✅' if ci2[0,0] <= -0.5 <= ci2[1,0] else '⚠️ 不含 −0.5'}")

    print("\n判定：")
    if lo > 0 and float((bs1[:, 0] <= -0.5).mean()) < 0.025:
        print("  ⭐ **复现成功** —— 两层噪声结构分离在独立数据集中重现。")
        print("     主张从'一套数据的观察'升级为'可复现的结构'。")
    else:
        print("  ⚠️ **未复现** —— 同样重要：说明两层分离可能是数据集特异的。")
        print("     主张须收窄为'该结构需逐数据集测量'，目标期刊回到 Genome Biology 层。")
    return c, e


# ---------------------------------------------------------------- L3
def mode_l3sat(d, n_boot):
    """Layer 3 saturation: link count against nucleus count, controlling depth."""
    m = d[d.method == "10xMulti"].copy()
    num = ["Multiome_Feature_linkages_detected", "Multiome_Linked_peaks",
           "Multiome_Linked_genes", "ATAC_Mean_raw_read_pairs_per_cell",
           "GEX_Mean_raw_reads_per_cell", "ATAC_Median_high_quality_fragments_per_cell"]
    for c in num:
        if c in m.columns:
            m[c] = pd.to_numeric(m[c], errors="coerce")

    lib = m.groupby(BATCH).agg(
        n_cells=(CTYPE, "size"),
        donor=(DONOR, "first"),
        link=("Multiome_Feature_linkages_detected", "first"),
        peaks=("Multiome_Linked_peaks", "first"),
        genes=("Multiome_Linked_genes", "first"),
        atac_depth=("ATAC_Mean_raw_read_pairs_per_cell", "first"),
        atac_frag=("ATAC_Median_high_quality_fragments_per_cell", "first"),
        gex_depth=("GEX_Mean_raw_reads_per_cell", "first")).dropna(
        subset=["n_cells", "link"])
    lib = lib[(lib.link > 0) & (lib.n_cells > 0)]
    print(f"{len(lib)} 个 Multiome 文库")
    print(f"  细胞数 {lib.n_cells.min():,}–{lib.n_cells.max():,}"
          f"（{lib.n_cells.max()/lib.n_cells.min():.0f}×）")
    print(f"  linkage {lib.link.min():,.0f}–{lib.link.max():,.0f}")

    b, ci, bs = fit_loglog(lib.n_cells.values, lib.link.values,
                           lib.donor.values, n_boot)
    print(f"\n⭐ log(linkage) ~ log(细胞数)  b = {b[0]:+.3f}  "
          f"95% CI [{ci[0,0]:+.3f}, {ci[1,0]:+.3f}]")
    print("   b≈1 = 线性增长、**完全未饱和**；b≈0 = 已饱和（达到生物学真值）")

    # Critical: rule out the depth confound before reading the exponent
    print("\n=== 混杂排除：控制测序深度后的偏相关 ===")
    for dcol, lab in [("atac_depth", "ATAC 每细胞原始读对"),
                      ("atac_frag", "ATAC 每细胞高质量片段"),
                      ("gex_depth", "GEX 每细胞原始读数")]:
        s = lib[lib[dcol].notna() & (lib[dcol] > 0)]
        if len(s) < 8:
            print(f"  {lab:<24} 样本不足，跳过"); continue
        r0 = np.corrcoef(np.log(s.n_cells), np.log(s.link))[0, 1]
        rp = partial_corr(s.n_cells.values, s.link.values, s[dcol].values)
        rd = np.corrcoef(np.log(s[dcol]), np.log(s.link))[0, 1]
        print(f"  {lab:<24} 深度↔linkage r={rd:+.3f} | "
              f"细胞数↔linkage 原始 r={r0:+.3f} → 偏相关 {rp:+.3f}"
              f"{'  ✅ 稳健' if abs(rp) > 0.4 else '  ⚠️ 大幅衰减，可能是深度效应'}")

    print("\n⚠️ 解读时必须说明的一点：")
    print("   cellranger-arc 的 linkage 数很可能是**未过滤**输出（本数据最大 2.5M）。")
    print("   若如此，b≈1 部分反映假阳性随细胞数线性增长——")
    print("   但这**同样支持** L3 结论：一个随细胞数线性增长、不收敛的量，")
    print("   不是在测量生物学。两种解释导向同一个结论。")
    return lib


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta", required=True)
    ap.add_argument("--mode", choices=["l1", "scaling", "l3sat"], default="l1")
    ap.add_argument("--allow-umi-proxy", action="store_true",
                    help="允许用中位 UMI 作 L2 代理（结果不可与 Ulm 直接比较）")
    ap.add_argument("--out", default="results")
    ap.add_argument("--n-boot", type=int, default=1500)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    use = [DONOR, BATCH, CTYPE, "method", "Number of UMIs", "Class",
           "Multiome_Feature_linkages_detected", "Multiome_Linked_peaks",
           "Multiome_Linked_genes", "ATAC_Mean_raw_read_pairs_per_cell",
           "GEX_Mean_raw_reads_per_cell",
           "ATAC_Median_high_quality_fragments_per_cell"]
    d = pd.read_csv(args.meta, usecols=lambda c: c in use, low_memory=False)
    print(f"读入 {len(d):,} 核\n")

    if args.mode in ("l1", "scaling"):
        c, e = mode_scaling(d, args.n_boot,
                            allow_umi_proxy=(args.mode == "scaling"
                                             and args.allow_umi_proxy))
        c.to_csv(os.path.join(args.out, "seaad_L1_pairs.csv"), index=False)
        e.to_csv(os.path.join(args.out, "seaad_L2_pairs.csv"), index=False)
    else:
        lib = mode_l3sat(d, args.n_boot)
        lib.to_csv(os.path.join(args.out, "seaad_L3_libraries.csv"))
    print(f"\n[done] 输出写入 {args.out}/")


if __name__ == "__main__":
    main()
