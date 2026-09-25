#!/usr/bin/env python3
"""Layer 3 的外部复现：NABEC/HBCC 前额叶多组学（Catching et al. Cell Rep 2026）。

口径照 p5_08 / p5_15：
  等化基质：每种细胞类型取 n 个核，均超过共同深度阈值，再多项式降采样到该深度
  连接：基因 TSS ±500 kb 内的 peak，Pearson 相关，细胞类型内统一 BH FDR<0.05
  启动子邻近：|peak 中点 − TSS| ≤ 3 kb
  富集 OR：检出 link 中近端/远端 的比数，比检验集中同一比数，Haldane–Anscombe 校正

⚠️ 与原队列的关键差别：这里每种细胞类型有数万到数十万个核，因此可以【扫 n】，
   把「150 个核时失败」推广为「到多少个核才不失败」。
"""
from __future__ import annotations
import argparse, json, logging, time
from pathlib import Path
import numpy as np, pandas as pd, h5py, scipy.sparse as sp

logger = logging.getLogger("L3")
WINDOW, PROMOTER, FDR, COLLIN = 500_000, 3_000, 0.05, 0.7

def obs_col(f, n):
    g = f["obs"][n]
    if isinstance(g, h5py.Group):
        return np.asarray(g["categories"][:]).astype(str)[g["codes"][:]]
    return np.asarray(g[:])

def var_index(f):
    v = f["var"]; return np.asarray(v[v.attrs.get("_index", "_index")][:]).astype(str)

def read_rows(path, key, rows, ncol):
    """按行号取稀疏子矩阵（行号需升序）。"""
    with h5py.File(path, "r") as f:
        ip = f[f"{key}/indptr"]; dat = f[f"{key}/data"]; idx = f[f"{key}/indices"]
        ptr = ip[:]
        parts = []
        for r in rows:
            s, e = int(ptr[r]), int(ptr[r+1])
            parts.append((dat[s:e].astype(np.float32), idx[s:e]))
    indptr = np.zeros(len(rows)+1, np.int64)
    indptr[1:] = np.cumsum([len(d) for d, _ in parts])
    return sp.csr_matrix((np.concatenate([d for d, _ in parts]),
                          np.concatenate([i for _, i in parts]), indptr),
                         shape=(len(rows), ncol))

def downsample(M, target, rng):
    """按行多项式降采样到 target 总计数。"""
    M = M.tocsr().astype(np.float64)
    for i in range(M.shape[0]):
        s, e = M.indptr[i], M.indptr[i+1]
        v = M.data[s:e]; tot = v.sum()
        if tot > target:
            M.data[s:e] = rng.multinomial(int(target), v/tot).astype(np.float64)
    M.eliminate_zeros(); return M

def bh_threshold(p, q, m_extra=0):
    """BH 阈值。m_extra 为额外的 p=1 检验数（零方差 peak），只进分母。"""
    p = np.sort(p[np.isfinite(p)]); m = p.size + int(m_extra)
    if p.size == 0: return -1.0
    ok = p <= q*np.arange(1, p.size+1)/m
    return float(p[ok][-1]) if ok.any() else -1.0

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rna", type=Path, default=Path("Data/final_rna_data.h5ad"))
    ap.add_argument("--atac", type=Path, default=Path("Data/final_atac_data.h5ad"))
    ap.add_argument("--tss", type=Path,
                    default=Path.home()/"Desktop/P5_VulnerableEpigenome/results/p5_tss_gencode_v32.csv")
    ap.add_argument("--celltype", required=True)
    ap.add_argument("--label-col", default="cell_type",
                    help="obs 中的注释列；用 leiden_2 可做细亚型粒度对照")
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--depth-atac", type=float, default=None,
                    help="强制 ATAC 深度（片段）；缺省由 n 决定共同阈值")
    ap.add_argument("--depth-rna", type=float, default=None, help="强制 RNA 深度（counts）")
    ap.add_argument("--min-det-frac", type=float, default=0.0,
                    help="与 p5_08 一致取 0.10：只测在该比例细胞中检出的基因")
    ap.add_argument("--fdr-all-pairs", action="store_true",
                    help="与 p5_08 一致：BH 分母含窗口内零方差 peak（p=1）")
    ap.add_argument("--log1p", action="store_true",
                    help="与 p5_08 一致：相关前对 RNA/ATAC 做 log1p")
    ap.add_argument("--tag", default="", help="输出文件名后缀")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=Path("results_methods/nabec_L3"))
    a = ap.parse_args(); a.out.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S")
    rng = np.random.default_rng(a.seed); t0 = time.time()

    with h5py.File(a.rna, "r") as fr:
        genes = var_index(fr); rna_ct = obs_col(fr, a.label_col)
        rna_idx = obs_col(fr, "_index").astype(str)
        rna_tot = obs_col(fr, "total_counts").astype(float)
        NG = int(fr["X"].attrs["shape"][1])
    with h5py.File(a.atac, "r") as fa:
        peaks = var_index(fa)
        atac_idx = np.array([x.rsplit("_", 1)[0] if x.count("_") >= 2 and
                             x.rsplit("_",1)[1] == x.rsplit("_",2)[1] else x
                             for x in obs_col(fa, "_index").astype(str)])
        atac_frag = obs_col(fa, "Unique_nr_frag").astype(float)
        NP = int(fa["X"].attrs["shape"][1])

    pos = {b: i for i, b in enumerate(atac_idx)}
    amap = np.array([pos.get(b, -1) for b in rna_idx])
    sel = np.flatnonzero((rna_ct == a.celltype) & (amap >= 0))
    logger.info("%s: %d 个核可用", a.celltype, sel.size)

    # 共同深度阈值：取该类型中第 n 大的深度（ATAC 与 RNA 各自），选满足两者的核
    af = atac_frag[amap[sel]]; rt = rna_tot[sel]
    # 两个阈值必须【联合】满足：按同一分位数放松，取仍能凑够 n 个核的最严分位
    lo, hi = 0.0, 1.0
    for _ in range(40):
        q = (lo + hi) / 2
        ta_, tr_ = np.quantile(af, 1-q), np.quantile(rt, 1-q)
        if int(((af >= ta_) & (rt >= tr_)).sum()) >= a.n: hi = q
        else: lo = q
    ta, tr = np.quantile(af, 1-hi), np.quantile(rt, 1-hi)
    if a.depth_atac or a.depth_rna:          # 强制深度：只要核超过该深度即可，再降采样到该深度
        ta = a.depth_atac if a.depth_atac else ta
        tr = a.depth_rna if a.depth_rna else tr
    elig = sel[(af >= ta) & (rt >= tr)]
    if elig.size < a.n:
        raise SystemExit(f"[FATAL] {a.celltype}: 阈值 ATAC {ta:.0f}/RNA {tr:.0f} 下仅 "
                         f"{elig.size} 个核，不足 {a.n}")
    pick = np.sort(rng.choice(elig, size=min(a.n, elig.size), replace=False))
    logger.info("  深度阈值 ATAC %.0f / RNA %.0f；选中 %d 个核 (%.0fs)", ta, tr, pick.size, time.time()-t0)

    R = downsample(read_rows(a.rna, "X", pick, NG), tr, rng)
    A = downsample(read_rows(a.atac, "X", np.sort(amap[pick]), NP), ta, rng)
    logger.info("  矩阵就绪 RNA %s / ATAC %s (%.0fs)", R.shape, A.shape, time.time()-t0)
    if a.log1p:                      # log1p(0)=0，可直接作用于 .data，稀疏性不变
        R = R.copy(); A = A.copy()
        R.data = np.log1p(R.data); A.data = np.log1p(A.data)
        logger.info("  已做 log1p 变换（与 p5_08 一致）")

    tss = pd.read_csv(a.tss)
    gi = {g: i for i, g in enumerate(genes)}
    tss = tss[tss.gene_name.isin(gi)].drop_duplicates("gene_name")
    pv = pd.DataFrame({"raw": peaks})
    spl = pv.raw.str.extract(r"^(chr[^:]+):(\d+)-(\d+)$")
    pv["chrom"] = spl[0]; pv["mid"] = ((spl[1].astype(float)+spl[2].astype(float))//2)
    pv["i"] = np.arange(len(pv)); pv = pv.dropna()

    # 相关用稀疏矩阵向量积算，避免把 n × 521k 的 ATAC 展开成稠密（大 n 会爆内存）
    n = pick.size
    A = A.tocsc()
    pmean = np.asarray(A.sum(0)).ravel() / n
    psq = np.asarray(A.multiply(A).sum(0)).ravel() / n
    pstd = np.sqrt(np.maximum(psq - pmean**2, 0.0))
    Rd = np.asarray(R.todense(), dtype=np.float32)          # n × 38k，可承受
    logger.info("  peak 非常数比例 %.1f%% (%.0fs)", 100*np.mean(pstd > 0), time.time()-t0)

    min_det = a.min_det_frac * n
    rows = []; n_all = 0; prox_all = 0
    for ch, tg in tss.groupby("chrom"):
        pc = pv[pv.chrom == ch].sort_values("mid")
        if pc.empty: continue
        mids = pc.mid.to_numpy(); pidx = pc.i.to_numpy()
        for gname, t in zip(tg.gene_name, tg.tss):
            lo, hi = np.searchsorted(mids, [t-WINDOW, t+WINDOW])
            if hi-lo < 2: continue
            y = Rd[:, gi[gname]].astype(np.float64)
            if (y > 0).sum() < min_det: continue
            sy = y.std()
            if sy <= 0: continue
            cols = pidx[lo:hi]
            ok = pstd[cols] > 0
            if ok.sum() < 2: continue
            dist_all = np.abs(mids[lo:hi]-t)
            n_all += dist_all.size; prox_all += int((dist_all <= PROMOTER).sum())
            cols = cols[ok]
            xty = np.asarray(A[:, cols].T @ y).ravel()      # 稀疏 matvec
            cov = xty/n - pmean[cols]*y.mean()
            r = cov / (pstd[cols]*sy) * (n/(n-1))
            rows.append((gname, cols, r, dist_all[ok]))
    logger.info("  相关计算完成：%d 个基因 (%.0fs)", len(rows), time.time()-t0)

    allr = np.concatenate([x[2] for x in rows])
    df = n - 2
    tstat = allr*np.sqrt(df/np.maximum(1-allr**2, 1e-12))
    from scipy import stats as st
    pvals = 2*st.t.sf(np.abs(tstat), df)
    n_const = (n_all - allr.size) if a.fdr_all_pairs else 0   # p5_08 把零方差对也计入分母
    thr = bh_threshold(pvals, FDR, m_extra=n_const)
    n_tested = allr.size
    prox_tested = int(np.concatenate([x[3] for x in rows]).__le__(PROMOTER).sum())
    if thr < 0:
        res = dict(celltype=a.celltype, n=int(n), n_tested=int(n_tested),
                   prox_tested=prox_tested, n_links=0, note="no link passes FDR")
    else:
        off = 0; a_=b_=0; ndec=0
        for gname, pidx, r, dist in rows:
            k = r.size; pg = pvals[off:off+k]; off += k
            sig = np.flatnonzero(pg <= thr)
            if sig.size == 0: continue
            order = sig[np.argsort(-np.abs(r[sig]))]
            keep = []
            for j in order[:200]:               # 每个基因最多考察 200 个显著 peak
                cj = np.asarray(A[:, pidx[j]].todense()).ravel()
                if all(abs(np.corrcoef(cj, np.asarray(A[:, pidx[m]].todense()).ravel())[0,1])
                       <= COLLIN for m in keep):
                    keep.append(j)
            ndec += len(keep)
            a_ += int((dist[sig] <= PROMOTER).sum()); b_ += int((dist[sig] > PROMOTER).sum())
        def OR(c_, d_):
            o = ((a_+.5)*(d_+.5))/((b_+.5)*(c_+.5))
            bs = rng.binomial(a_+b_, (a_+.5)/(a_+b_+1), 2000)
            ci = np.percentile(((bs+.5)*(d_+.5))/((a_+b_-bs+.5)*(c_+.5)), [2.5, 97.5])
            return float(o), float(ci[0]), float(ci[1])
        # 口径 A（与原队列 p5_15 一致）：检验集 = 窗口内全部 (peak, TSS) 对
        oA = OR(prox_all, n_all - prox_all)
        # 口径 B：仅计入非零方差 peak，即理论上可被检出的对
        oB = OR(prox_tested, n_tested - prox_tested)
        res = dict(celltype=a.celltype, n=int(n), depth_atac=float(ta), depth_rna=float(tr),
                   n_tested_all=int(n_all), prox_tested_all=int(prox_all),
                   n_tested_var=int(n_tested), prox_tested_var=int(prox_tested),
                   n_links=int(a_+b_), prox_links=int(a_), n_independent=int(ndec),
                   OR=oA[0], lo=oA[1], hi=oA[2],
                   OR_varpeaks=oB[0], lo_varpeaks=oB[1], hi_varpeaks=oB[2])
    json.dump(res, open(a.out/f"{a.celltype}_n{a.n}{a.tag}_seed{a.seed}.json", "w"), indent=2)
    logger.info("结果: %s", json.dumps(res, ensure_ascii=False))
    logger.info("用时 %.1f 分", (time.time()-t0)/60)

if __name__ == "__main__":
    main()
