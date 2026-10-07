#!/usr/bin/env python3
"""原队列核数梯队 · C1 对照：把 single-cell 相关换成 ArchR 式 KNN 聚合相关。

与 `primary_ladder_run.py` 的唯一差别是 peak–gene 相关在**元细胞**（metacell）
之间计算，而不是在单核之间计算。细胞抽取、深度降采样、±500 kb 窗口、3 kb 启动子
定义、BH FDR 0.05（含同样的 m_extra 校正）、Haldane–Anscombe 校正的 OR 与
bootstrap 区间全部逐字保持一致——这是对照的全部意义。

聚合方式照 ArchR `addPeak2GeneLinks`：在 ATAC 的 LSI 空间里取 KNN 聚合体，
按 overlapCutoff 去掉高度重叠者，组内 raw counts 求和后 log1p。

`--match-tested` 把被检验的 gene×peak 集合锁定为 single-cell 臂的集合，
使两臂共享分母；不加时用聚合臂自己的检出规则（即真跑 ArchR 会得到的集合）。

两个零假设：
  --null permute  把基因向量在观测之间打乱。对聚合臂**不是**合适的零假设：
                  它同时破坏了 y 在 KNN 图上的平滑性，而 x 仍然平滑，
                  于是低估假阳性。保留它只为显示这一点。
  --null trans    把每个基因配到**另一条染色体**上另一个基因的窗口 peak 集合。
                  两臂的图平滑性都原样保留，只打断配对，因此这是诊断
                  「聚合把多少个链接变成假阳性」的正确零假设。窗口大小与
                  近端/远端比例的分布也一并保留，OR 的分母因此可比。

`--linker single` 用单核相关跑同一套代码，使两臂除聚合外无任何差别。

用法：
  python3 primary_ladder_aggregated.py --dir <d> --celltype X --n 150 --k 25
  python3 primary_ladder_aggregated.py --dir <d> --celltype X --n 150 --linker single
"""
from __future__ import annotations

import argparse
import json
import logging
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import sparse, stats as st
from scipy.spatial import cKDTree

WINDOW, PROMOTER, FDR = 500_000, 3_000, 0.05
# 距离带。启动子二分法只用上检出链接里 ~5% 的那一小撮，带状分解把全部链接都用上，
# 因此位置信息的检验力高得多。
DIST_BANDS = (0, 3_000, 10_000, 50_000, 100_000, 250_000, 500_001)
DEPTH_ATAC, DEPTH_RNA = 5564, 5265
N_SVD = 30
OVERLAP_CUTOFF = 0.8
MIN_DET_FRAC = 0.10
# ArchR `getPeak2GeneLinks` 的实际默认阈值（见 R/IntegrativeAnalysis.R）。
# 只用 BH FDR 0.05、不设相关系数下限，会放出远多于 ArchR 真实默认的链接，
# 那是在打稻草人；`--archr-defaults` 照抄它自己的阈值。
ARCHR_COR_CUTOFF = 0.45
ARCHR_FDR_CUTOFF = 1e-4
ARCHR_VAR_QUANTILE = 0.25
ARCHR_SCALE_TO = 1e4
ARCHR_MAXDIST = 250_000

logger = logging.getLogger("ladder_agg")


def bh_threshold(p: np.ndarray, q: float, m_extra: int = 0) -> float:
    """BH 阈值。m_extra 计入被跳过（无变异）的检验，与 single-cell 臂一致。"""
    p = np.sort(p[np.isfinite(p)])
    m = p.size + int(m_extra)
    if p.size == 0:
        return -1.0
    ok = p <= q * np.arange(1, p.size + 1) / m
    return float(p[ok][-1]) if ok.any() else -1.0


def downsample(M: sparse.spmatrix, target: int, rng: np.random.Generator) -> sparse.csc_matrix:
    """按列多项降采样到 target 计数。与 primary_ladder_run.py 逐字相同。"""
    M = M.tocsc().astype(np.int64)
    out = M.copy()
    for j in range(M.shape[1]):
        s, e = M.indptr[j], M.indptr[j + 1]
        v = M.data[s:e]
        tot = v.sum()
        if tot <= target:
            continue
        out.data[s:e] = rng.multinomial(int(target), v / tot)
    out.eliminate_zeros()
    return out.tocsc()


def lsi_embedding(A: sparse.csc_matrix, n_comp: int, rng: np.random.Generator) -> np.ndarray:
    """TF-IDF + SVD，丢掉第 1 个成分（惯例：与测序深度共线）。A 为 细胞 × peak。"""
    counts = A.tocsr().astype(np.float64)
    n_cells = counts.shape[0]
    peak_sums = np.asarray(counts.sum(0)).ravel()
    keep = peak_sums > 0
    counts = counts[:, keep]
    peak_sums = peak_sums[keep]
    cell_sums = np.asarray(counts.sum(1)).ravel()
    cell_sums[cell_sums == 0] = 1.0
    # TF-IDF：tf = count / cell_total，idf = n_cells / peak_total
    tf = counts.multiply(1.0 / cell_sums[:, None]).tocsc()
    idf = n_cells / peak_sums
    X = tf.multiply(idf[None, :]).tocsr()
    X.data = np.log1p(X.data * 1e4)
    k = min(n_comp + 1, min(X.shape) - 1)
    if k < 2:
        raise SystemExit("[FATAL] 矩阵太小，无法做 LSI")
    rs = int(rng.integers(0, 2**31 - 1))
    U, S, _ = sparse.linalg.svds(X, k=k, random_state=rs)
    order = np.argsort(-S)
    emb = (U * S)[:, order]
    return emb[:, 1:]  # 丢第 1 成分


def knn_aggregates(emb: np.ndarray, k: int, n_seeds: int,
                   overlap_cutoff: float, rng: np.random.Generator) -> List[np.ndarray]:
    """ArchR 式聚合：随机种子细胞 → 其 k 近邻为一个聚合体 → 去掉高度重叠者。"""
    n = emb.shape[0]
    z = (emb - emb.mean(0)) / np.maximum(emb.std(0), 1e-12)
    seeds = rng.permutation(n)[:min(n_seeds, n)]
    # 只对被选中的种子查 k 近邻，避免 n×n 距离矩阵
    tree = cKDTree(z)
    _, nbr_seeds = tree.query(z[seeds], k=k, workers=-1)
    nbr_seeds = np.atleast_2d(nbr_seeds)
    kept: List[np.ndarray] = []
    kept_sets: List[set] = []
    for si in range(len(seeds)):
        cand = nbr_seeds[si]
        cs = set(cand.tolist())
        if any(len(cs & ks) / k > overlap_cutoff for ks in kept_sets):
            continue
        kept.append(cand)
        kept_sets.append(cs)
    return kept


def aggregate_counts(M: sparse.csc_matrix, groups: List[np.ndarray]) -> sparse.csr_matrix:
    """把 细胞 × 特征 的 counts 按聚合体求和，返回 聚合体 × 特征。"""
    n_cells = M.shape[0]
    rows, cols = [], []
    for gi, members in enumerate(groups):
        rows.extend([gi] * len(members))
        cols.extend(members.tolist())
    S = sparse.csr_matrix(
        (np.ones(len(rows), dtype=np.float32), (rows, cols)),
        shape=(len(groups), n_cells))
    return (S @ M.tocsr()).tocsr()


def haldane_or(a_: int, b_: int, c_: int, d_: int,
               rng: np.random.Generator, n_boot: int = 2000) -> Tuple[float, float, float]:
    """Haldane–Anscombe 校正的 OR，区间由检出链接数的二项 bootstrap 给出。"""
    o = ((a_ + .5) * (d_ + .5)) / ((b_ + .5) * (c_ + .5))
    if a_ + b_ == 0:
        return float(o), float("nan"), float("nan")
    bs = rng.binomial(a_ + b_, (a_ + .5) / (a_ + b_ + 1), n_boot)
    ci = np.percentile(((bs + .5) * (d_ + .5)) / ((a_ + b_ - bs + .5) * (c_ + .5)), [2.5, 97.5])
    return float(o), float(ci[0]), float(ci[1])


def load_features(features_dir: Path) -> Tuple[pd.DataFrame, pd.Series]:
    peaks = pd.read_csv(features_dir / "Multiome_Dataset_ATAC_features.tsv",
                        sep="\t", header=None)[0].astype(str)
    gpath = features_dir / "Multiome_Dataset_RNA_features.tsv"
    if not gpath.exists():
        gpath = features_dir / "Multiome_Dataset_RNA_barcodes.tsv"
    genes = pd.read_csv(gpath, sep="\t", header=None)[0].astype(str)
    pv = pd.DataFrame({"raw": peaks})
    spl = pv.raw.str.extract(r"^(chr[^:\-]+)[:\-](\d+)[\-](\d+)$")
    pv["chrom"] = spl[0]
    pv["mid"] = (spl[1].astype(float) + spl[2].astype(float)) // 2
    pv["i"] = np.arange(len(pv))
    return pv.dropna(), genes


def trans_pairing(windows: Dict[str, Tuple[np.ndarray, np.ndarray, str]],
                  rng: np.random.Generator
                  ) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """把每个基因重新配到另一条染色体上另一个基因的窗口 peak 集合。

    窗口大小与近端/远端比例的分布原样保留（只是换了持有者），两臂的图平滑性
    也不受影响，所以这是诊断聚合带来的假阳性的正确零假设。
    """
    names = list(windows)
    out: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    for g in names:
        own_chrom = windows[g][2]
        for _ in range(50):
            h = names[int(rng.integers(0, len(names)))]
            if windows[h][2] != own_chrom:
                out[g] = (windows[h][0], windows[h][1])
                break
    return out


def correlate(Rd: np.ndarray, A: sparse.csc_matrix, pv: pd.DataFrame,
              tss: pd.DataFrame, gi: Dict[str, int], n_obs: int,
              min_det: float, restrict: Optional[Dict[str, np.ndarray]],
              rng: Optional[np.random.Generator],
              null: str = "none", window: int = WINDOW
              ) -> Tuple[List[Tuple[np.ndarray, ...]], int, int, Dict[str, np.ndarray]]:
    """在 n_obs 个观测（单核或聚合体）之间算窗口内 peak–gene 相关。

    restrict 非 None 时只检验给定的 gene→peak 列集合，使两臂共享分母。
    null="permute" 打乱基因向量；null="trans" 换成异染色体基因的窗口。
    """
    pmean = np.asarray(A.sum(0)).ravel() / n_obs
    psq = np.asarray(A.multiply(A).sum(0)).ravel() / n_obs
    pstd = np.sqrt(np.maximum(psq - pmean ** 2, 0.0))
    # 方差在全部 peak / 基因上的分位秩（ArchR 的 varCutOffATAC / varCutOffRNA）
    pvar = pstd ** 2
    varq_peak = st.rankdata(pvar, method="average") / pvar.size
    gvar = Rd.var(axis=0)
    varq_gene = st.rankdata(gvar, method="average") / max(gvar.size, 1)

    # ---- 第一遍：确定每个基因的窗口（cols, dist, chrom） ----
    # n_win/prox_win 是窗口内**全部**对的计数（未按 peak 变异过滤），
    # 供 m_extra 与 all-peaks 分母使用，与 single-cell 臂的定义一致。
    windows: Dict[str, Tuple[np.ndarray, np.ndarray, str]] = {}
    n_win = prox_win = 0
    for ch, tg in tss.groupby("chrom"):
        pc = pv[pv.chrom == ch].sort_values("mid")
        if pc.empty:
            continue
        mids = pc.mid.to_numpy()
        pidx = pc.i.to_numpy()
        for gname, t in zip(tg.gene_name, tg.tss):
            if restrict is not None and gname not in restrict:
                continue
            lo, hi = np.searchsorted(mids, [t - window, t + window])
            if hi - lo < 2:
                continue
            y = Rd[:, gi[gname]]
            if restrict is None and (y > 0).sum() < min_det:
                continue
            cols = pidx[lo:hi]
            dist = np.abs(mids[lo:hi] - t)
            if restrict is not None:
                sel = np.isin(cols, restrict[gname])
                cols, dist = cols[sel], dist[sel]
            keep = pstd[cols] > 0
            if keep.sum() < 2:
                continue
            n_win += dist.size
            prox_win += int((dist <= PROMOTER).sum())
            windows[gname] = (cols[keep], dist[keep], str(ch))

    if not windows:
        return [], 0, 0, {}

    # ---- 零假设：换配对 ----
    swapped = trans_pairing(windows, rng) if null == "trans" and rng is not None else None

    # ---- 第二遍：相关 ----
    rows: List[Tuple[np.ndarray, ...]] = []
    tested: Dict[str, np.ndarray] = {}
    for gname, (cols, dist, _ch) in windows.items():
        if swapped is not None:
            if gname not in swapped:
                continue
            cols, dist = swapped[gname]
        y = Rd[:, gi[gname]].astype(np.float64)
        if null == "permute" and rng is not None:
            y = y[rng.permutation(y.size)]
        sy = y.std()
        if sy <= 0:
            continue
        xty = np.asarray(A[:, cols].T @ y).ravel()
        cov = xty / n_obs - pmean[cols] * y.mean()
        r = cov / (pstd[cols] * sy) * (n_obs / (n_obs - 1))
        rows.append((cols, r, dist, varq_peak[cols],
                     float(varq_gene[gi[gname]]), _ch))
        tested[gname] = cols
    return rows, n_win, prox_win, tested


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="results_methods/primary_ladder")
    ap.add_argument("--celltype", required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--linker", choices=("aggregate", "single"), default="aggregate")
    ap.add_argument("--k", type=int, default=0, help="每个聚合体的细胞数（--linker aggregate 必填）")
    ap.add_argument("--n-seeds", type=int, default=500, help="聚合体种子数（ArchR knnIteration）")
    ap.add_argument("--match-tested", action="store_true",
                    help="锁定为 single-cell 臂检验过的 gene×peak 集合")
    ap.add_argument("--null", choices=("none", "permute", "trans"), default="none")
    ap.add_argument("--trans-rep", type=int, default=0,
                    help="trans 配对的重复编号。改变它只改变配对的随机抽取，"
                         "不改变细胞抽样，用来检验位置 OR 对单次配对是否敏感")
    ap.add_argument("--archr-defaults", action="store_true",
                    help="照 ArchR getPeak2GeneLinks 的默认阈值取链接："
                         "r>=0.45、FDR<=1e-4、两侧方差分位>0.25、scaleTo 1e4 + log2(x+1)")
    ap.add_argument("--depth-mult", type=float, default=1.0,
                    help="把等深度目标乘以该倍数（核数 × 深度响应面用）")
    ap.add_argument("--depth-mult-atac", type=float, default=0.0,
                    help="只缩放 ATAC 深度的倍数；默认 0 表示沿用 --depth-mult。"
                         "已发表的深度序列缩放的正是 ATAC 而把 RNA 固定在 5,265，"
                         "所以要复现那个操作必须用这个参数而不是 --depth-mult")
    ap.add_argument("--depth-mult-rna", type=float, default=0.0,
                    help="只缩放 RNA 深度的倍数；默认 0 表示沿用 --depth-mult")
    ap.add_argument("--pool-mult", type=float, default=0.0,
                    help="资格门槛用这个倍数，而抽到的核再降采样到 --depth-mult。"
                         "默认 0 表示与 --depth-mult 相同。把它固定在高倍数、让 "
                         "--depth-mult 变化，就能在同一批核上做深度序列，"
                         "从而把深度本身与深核是被挑出来的这一混淆分开")
    ap.add_argument("--require-depth", action="store_true",
                    help="只从总计数达到深度目标的核里抽样。不加时沿用已发表梯队的"
                         "做法（随机抽样后降采样，未达标的核保持原深度）；在 >1x "
                         "的深度上必须加，否则等深度会被静默破坏")
    ap.add_argument("--window", type=int, default=WINDOW,
                    help=f"±窗口，默认 {WINDOW}（正文口径）；ArchR 自己的默认是 {ARCHR_MAXDIST}")
    ap.add_argument("--tss", default="Resolution_limit/data/derived_results/p5_tss_gencode_v32.csv")
    ap.add_argument("--features", default="Data/Other_Datasets/Multiome_Dataset/files")
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S")
    D = Path(a.dir)
    t0 = time.time()

    # ---- 与 single-cell 臂逐字相同的前半段，rng 消耗顺序也相同 ----
    rng = np.random.default_rng(a.seed)
    cells = pd.read_csv(D / "cells.csv")
    idx = np.flatnonzero((cells.WNN_L4 == a.celltype).to_numpy())
    if idx.size < a.n:
        raise SystemExit(f"[FATAL] {a.celltype}: 仅 {idx.size} 核，不足 {a.n}")
    if a.linker == "aggregate":
        if a.k <= 0:
            raise SystemExit("[FATAL] --linker aggregate 需要 --k")
        if a.k > a.n // 3:
            raise SystemExit(f"[FATAL] k={a.k} 对 n={a.n} 过大（要求 k <= n/3），聚合体会退化")
    mult_a = a.depth_mult_atac if a.depth_mult_atac > 0 else a.depth_mult
    mult_r = a.depth_mult_rna if a.depth_mult_rna > 0 else a.depth_mult
    target_a = int(round(DEPTH_ATAC * mult_a))
    target_r = int(round(DEPTH_RNA * mult_r))
    if a.require_depth:
        pool_mult = a.pool_mult if a.pool_mult > 0 else max(mult_a, mult_r)
        if pool_mult < max(mult_a, mult_r):
            raise SystemExit("[FATAL] --pool-mult 不能小于所要求的深度倍数")
        gate_a = int(round(DEPTH_ATAC * pool_mult))
        gate_r = int(round(DEPTH_RNA * pool_mult))
        Afull = sparse.load_npz(D / "atac.npz").tocsc()
        Rfull = sparse.load_npz(D / "rna.npz").tocsc()
        tot_a = np.asarray(Afull.sum(0)).ravel()
        tot_r = np.asarray(Rfull.sum(0)).ravel()
        ok = idx[(tot_a[idx] >= gate_a) & (tot_r[idx] >= gate_r)]
        if ok.size < a.n:
            raise SystemExit(
                f"[FATAL] {a.celltype} 在 {pool_mult}x 资格门槛（ATAC {gate_a} / "
                f"RNA {gate_r}）下只有 {ok.size} 个达标核，不足 {a.n}")
        logger.info("达标核 %d / %d（资格 %.1fx，降采样到 %.1fx）",
                    ok.size, idx.size, pool_mult, a.depth_mult)
        pick = np.sort(rng.choice(ok, size=a.n, replace=False))
        A = Afull[:, pick]
        R = Rfull[:, pick]
        del Afull, Rfull
    else:
        pick = np.sort(rng.choice(idx, size=a.n, replace=False))
        A = sparse.load_npz(D / "atac.npz").tocsc()[:, pick]
        R = sparse.load_npz(D / "rna.npz").tocsc()[:, pick]
    A = downsample(A, target_a, rng)
    R = downsample(R, target_r, rng)
    A = A.T.tocsc()
    R = R.T.tocsr()                                  # 细胞 × 特征
    logger.info("矩阵就绪 RNA %s / ATAC %s (%.0fs)", R.shape, A.shape, time.time() - t0)

    pv, genes = load_features(Path(a.features))
    tss = pd.read_csv(a.tss)
    gi = {g: i for i, g in enumerate(genes)}
    tss = tss[tss.gene_name.isin(gi)].drop_duplicates("gene_name")

    # ---- single-cell 臂：始终先算一遍，用来锁定共享的被检验集合 ----
    A_sc = A.astype(np.float32)
    A_sc.data = np.log1p(A_sc.data)
    Rd_sc = np.log1p(np.asarray(R.todense(), dtype=np.float32))
    _, n_all_sc, prox_all_sc, tested_sc = correlate(
        Rd_sc, A_sc, pv, tss, gi, a.n, MIN_DET_FRAC * a.n, None, None)
    logger.info("single-cell 臂：%d 基因 / %d 对", len(tested_sc), n_all_sc)

    nrng = (np.random.default_rng(20_000 + a.seed + 1_000 * a.trans_rep)
            if a.null != "none" else None)
    restrict = tested_sc if a.match_tested else None

    if a.linker == "aggregate":
        arng = np.random.default_rng(10_000 + a.seed)
        emb = lsi_embedding(A, N_SVD, arng)
        groups = knn_aggregates(emb, a.k, a.n_seeds, OVERLAP_CUTOFF, arng)
        n_agg = len(groups)
        if n_agg < 10:
            raise SystemExit(f"[FATAL] 只剩 {n_agg} 个聚合体，无法做相关")
        Aa_raw = aggregate_counts(A, groups).tocsc()
        Ra_raw = np.asarray(aggregate_counts(R, groups).todense(), dtype=np.float64)
        if a.archr_defaults:
            # ArchR: t(t(M)/colSums(M)) * scaleTo，然后 log2(M+1)。
            # 本研究的核已做等深度，每个聚合体恰好 k 个核，所以 colSums 本来就相等；
            # 仍然照抄这一步，免得归一化成为对照的第二个差异。
            acs = np.asarray(Aa_raw.sum(1)).ravel()
            acs[acs == 0] = 1.0
            Aa = Aa_raw.multiply((ARCHR_SCALE_TO / acs)[:, None]).tocsc().astype(np.float32)
            Aa.data = np.log2(Aa.data + 1.0)
            rcs = Ra_raw.sum(1)
            rcs[rcs == 0] = 1.0
            Rd_use = np.log2(Ra_raw * (ARCHR_SCALE_TO / rcs)[:, None] + 1.0).astype(np.float32)
        else:
            Aa = Aa_raw.astype(np.float32)
            Aa.data = np.log1p(Aa.data)
            Rd_use = np.log1p(Ra_raw).astype(np.float32)
        A_use, n_obs = Aa, n_agg
        logger.info("聚合体 %d 个 × 每个 %d 核（重叠上限 %.2f）", n_agg, a.k, OVERLAP_CUTOFF)
    else:
        A_use, Rd_use, n_obs, n_agg = A_sc, Rd_sc, a.n, 0

    rows, n_all, prox_all, _ = correlate(
        Rd_use, A_use, pv, tss, gi, n_obs, MIN_DET_FRAC * n_obs, restrict, nrng,
        a.null, a.window)
    if not rows:
        raise SystemExit("[FATAL] 没有可检验的对")

    allr = np.concatenate([x[1] for x in rows])
    df = n_obs - 2
    tstat = allr * np.sqrt(df / np.maximum(1 - allr ** 2, 1e-12))
    pvals = 2 * st.t.sf(np.abs(tstat), df)
    prox_tested = int((np.concatenate([x[2] for x in rows]) <= PROMOTER).sum())

    if a.archr_defaults:
        # ArchR 的 FDR 是在全部被检验对上做 BH，没有 m_extra 的概念
        order = np.argsort(pvals)
        m = pvals.size
        fdr = np.empty(m)
        fdr[order] = np.minimum.accumulate(
            (pvals[order] * m / np.arange(1, m + 1))[::-1])[::-1]
        thr = FDR  # 占位，ArchR 分支不用 BH 阈值
    else:
        thr = bh_threshold(pvals, FDR, m_extra=n_all - allr.size)
        fdr = None

    a_ = b_ = 0
    off = 0
    link_dist: List[np.ndarray] = []
    # 每条染色体：[近端链接, 远端链接, 近端被检验, 远端被检验]
    per_chrom: Dict[str, List[int]] = {}
    for cols, _r, dist, _vqa, _vqr, chrom in rows:
        cb = per_chrom.setdefault(chrom, [0, 0, 0, 0])
        cb[2] += int((dist <= PROMOTER).sum())
        cb[3] += int((dist > PROMOTER).sum())
    for cols, r, dist, vqa, vqr, chrom in rows:
        k = r.size
        sl = slice(off, off + k)
        off += k
        if a.archr_defaults:
            # r >= 0.45（只要正相关）且 FDR <= 1e-4 且两侧方差分位 > 0.25
            sig = np.flatnonzero(
                (r >= ARCHR_COR_CUTOFF) & (fdr[sl] <= ARCHR_FDR_CUTOFF)
                & (vqa > ARCHR_VAR_QUANTILE) & (vqr > ARCHR_VAR_QUANTILE))
        else:
            if thr < 0:
                continue
            sig = np.flatnonzero(pvals[sl] <= thr)
        if sig.size == 0:
            continue
        np_ = int((dist[sig] <= PROMOTER).sum())
        nd_ = int(sig.size) - np_
        a_ += np_
        b_ += nd_
        cb = per_chrom.setdefault(chrom, [0, 0, 0, 0])
        cb[0] += np_
        cb[1] += nd_
        link_dist.append(dist[sig])
    dl = np.concatenate(link_dist) if link_dist else np.zeros(0)
    tested_dist = np.concatenate([x[2] for x in rows])
    band_links = np.histogram(dl, bins=DIST_BANDS)[0].tolist()
    band_tested = np.histogram(tested_dist, bins=DIST_BANDS)[0].tolist()

    brng = np.random.default_rng(a.seed)
    oA = haldane_or(a_, b_, prox_all, n_all - prox_all, brng)
    oB = haldane_or(a_, b_, prox_tested, allr.size - prox_tested, brng)
    res = dict(cohort="primary", celltype=a.celltype, n=a.n, seed=a.seed,
               linker=("knn_aggregate" if a.linker == "aggregate" else "single_cell"),
               k=int(a.k), n_aggregates=int(n_agg), n_obs=int(n_obs),
               match_tested=bool(a.match_tested), null=a.null,
               archr_defaults=bool(a.archr_defaults), window=int(a.window),
               depth_atac=target_a, depth_rna=target_r,
               depth_mult=float(a.depth_mult), depth_mult_atac=float(mult_a),
               depth_mult_rna=float(mult_r), require_depth=bool(a.require_depth),
               pool_mult=float(a.pool_mult if a.pool_mult > 0 else max(mult_a, mult_r)),
               n_tested_all=int(n_all), prox_tested_all=int(prox_all),
               n_tested_var=int(allr.size), prox_tested_var=int(prox_tested),
               n_links=int(a_ + b_), prox_links=int(a_),
               trans_rep=int(a.trans_rep), per_chrom=per_chrom,
               band_edges=list(DIST_BANDS), band_links=band_links,
               band_tested=band_tested,
               median_link_dist=(float(np.median(dl)) if dl.size else float("nan")),
               OR=oA[0], lo=oA[1], hi=oA[2],
               OR_varpeaks=oB[0], lo_varpeaks=oB[1], hi_varpeaks=oB[2],
               sc_n_tested_all=int(n_all_sc), sc_prox_tested_all=int(prox_all_sc))
    tag = (f"k{a.k}" if a.linker == "aggregate" else "sc")
    tag += ("_archr" if a.archr_defaults else "")
    tag += ("" if a.window == WINDOW else f"_w{a.window // 1000}kb")
    tag += ("" if mult_a == 1.0 else f"_a{mult_a:g}x")
    tag += ("" if mult_r == 1.0 else f"_r{mult_r:g}x")
    tag += ("" if a.pool_mult in (0.0, a.depth_mult) else f"_pool{a.pool_mult:g}x")
    tag += ("_matched" if a.match_tested else "") + ("" if a.null == "none" else f"_{a.null}")
    tag += ("" if a.trans_rep == 0 else f"_rep{a.trans_rep}")
    out = Path(a.out) if a.out else D / f"c1_{a.celltype}_n{a.n}_s{a.seed}_{tag}.json"
    out.write_text(json.dumps(res, indent=2))
    logger.info("结果: %s", json.dumps(res))


if __name__ == "__main__":
    main()
