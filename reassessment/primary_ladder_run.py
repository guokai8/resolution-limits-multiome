#!/usr/bin/env python3
"""Primary-cohort nucleus ladder, step 2: run the Layer 3 analysis.

Operates on the sub-matrices written by primary_ladder_extract.py, using
exactly the parameters of p5_08 and the external NABEC run, so that a rung of
this ladder is comparable with both.

The sequence is: draw n nuclei, downsample both modalities to the common
depths, correlate every gene against every peak within the window, apply
Benjamini-Hochberg across the cell type, then ask what fraction of the
surviving links are promoter-proximal relative to the background.

Only the nucleus count varies between rungs. Everything else -- window,
promoter definition, FDR, collinearity cutoff, depth -- is frozen at module
level so that no rung can differ in any other respect.
"""
import argparse, json, logging, time
from pathlib import Path
import numpy as np, pandas as pd
from scipy import sparse, stats as st

WINDOW, PROMOTER, FDR, COLLIN = 500_000, 3_000, 0.05, 0.7
logger = logging.getLogger("ladder")

def bh_threshold(p, q, m_extra=0):
    """Largest p-value passing Benjamini-Hochberg at level q.

    m_extra adds pairs that were never tested -- peaks with no variance at this
    depth -- to the denominator. Including them is the stricter convention and
    is the one the paper uses throughout; it matters because the number of
    untestable peaks itself depends on depth.

    Returns -1.0 when nothing passes, which the caller treats as "no links".
    """
    p = np.sort(p[np.isfinite(p)]); m = p.size + int(m_extra)
    if p.size == 0: return -1.0
    ok = p <= q*np.arange(1, p.size+1)/m
    return float(p[ok][-1]) if ok.any() else -1.0

def downsample(M, target, rng):
    """Multinomially downsample each column to `target` total counts.

    Per column, so every nucleus ends at the same depth. Columns already below
    target are left alone rather than being discarded or scaled up, which is
    why the extraction step pre-filters on the same thresholds.
    """
    M = M.tocsc().astype(np.int64); out = M.copy()
    for j in range(M.shape[1]):
        s, e = M.indptr[j], M.indptr[j+1]
        v = M.data[s:e]; tot = v.sum()
        if tot <= target: continue
        out.data[s:e] = rng.multinomial(int(target), v/tot)
    out.eliminate_zeros(); return out.tocsc()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="results_methods/primary_ladder")
    ap.add_argument("--celltype", required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--tss", default="Resolution_limit/data/derived_results/p5_tss_gencode_v32.csv")
    ap.add_argument("--features", default="Data/Other_Datasets/Multiome_Dataset/files")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S")
    D = Path(a.dir); rng = np.random.default_rng(a.seed); t0 = time.time()

    cells = pd.read_csv(D/"cells.csv")
    idx = np.flatnonzero((cells.WNN_L4 == a.celltype).to_numpy())
    if idx.size < a.n: raise SystemExit(f"[FATAL] {a.celltype}: 仅 {idx.size} 核，不足 {a.n}")
    pick = np.sort(rng.choice(idx, size=a.n, replace=False))

    A = sparse.load_npz(D/"atac.npz").tocsc()[:, pick]
    R = sparse.load_npz(D/"rna.npz").tocsc()[:, pick]
    A = downsample(A, 5564, rng); R = downsample(R, 5265, rng)
    A = A.T.tocsc(); R = R.T.tocsr()      # transpose to cells x features
    logger.info("矩阵就绪 RNA %s / ATAC %s (%.0fs)", R.shape, A.shape, time.time()-t0)

    # log1p on both modalities. The depths are already equalised, so this is
    # for variance stabilisation rather than normalisation.
    A = A.astype(np.float32); A.data = np.log1p(A.data)
    Rd = np.asarray(R.todense(), dtype=np.float32); Rd = np.log1p(Rd)

    F = Path(a.features)
    peaks = pd.read_csv(F/"Multiome_Dataset_ATAC_features.tsv", sep="\t", header=None)[0].astype(str)
    genes = pd.read_csv(F/"Multiome_Dataset_RNA_features.tsv", sep="\t", header=None)[0].astype(str) \
            if (F/"Multiome_Dataset_RNA_features.tsv").exists() else None
    if genes is None:
        genes = pd.read_csv(F/"Multiome_Dataset_RNA_barcodes.tsv", sep="\t", header=None)[0].astype(str)
    pv = pd.DataFrame({"raw": peaks})
    spl = pv.raw.str.extract(r"^(chr[^:\-]+)[:\-](\d+)[\-](\d+)$")
    pv["chrom"] = spl[0]; pv["mid"] = ((spl[1].astype(float)+spl[2].astype(float))//2)
    pv["i"] = np.arange(len(pv)); pv = pv.dropna()
    logger.info("peak 可解析 %d / %d", len(pv), len(peaks))

    tss = pd.read_csv(a.tss); gi = {g: i for i, g in enumerate(genes)}
    tss = tss[tss.gene_name.isin(gi)].drop_duplicates("gene_name")

    # Peak means and standard deviations, computed once over all peaks from the
    # sparse sums rather than densifying: E[x^2] - E[x]^2, clipped at 0 for
    # floating-point noise. A gene must be detected in at least 10% of nuclei.
    n = a.n; min_det = 0.10*n
    pmean = np.asarray(A.sum(0)).ravel()/n
    psq = np.asarray(A.multiply(A).sum(0)).ravel()/n
    pstd = np.sqrt(np.maximum(psq - pmean**2, 0.0))
    rows=[]; n_all=0; prox_all=0
    for ch, tg in tss.groupby("chrom"):
        pc = pv[pv.chrom == ch].sort_values("mid")
        if pc.empty: continue
        mids = pc.mid.to_numpy(); pidx = pc.i.to_numpy()
        for gname, t in zip(tg.gene_name, tg.tss):
            # Peaks are sorted by midpoint per chromosome, so the window is a
            # pair of binary searches rather than a scan.
            lo, hi = np.searchsorted(mids, [t-WINDOW, t+WINDOW])
            if hi-lo < 2: continue
            y = Rd[:, gi[gname]].astype(np.float64)
            if (y > 0).sum() < min_det: continue
            sy = y.std()
            if sy <= 0: continue
            cols = pidx[lo:hi]; ok = pstd[cols] > 0
            if ok.sum() < 2: continue
            dist = np.abs(mids[lo:hi]-t)
            # Two backgrounds are tracked. n_all counts every eligible pair in
            # the window; the var-peak counts below drop peaks with zero
            # variance, which cannot be tested at this depth. Reporting both is
            # what makes the convention's effect on the answer visible.
            n_all += dist.size; prox_all += int((dist <= PROMOTER).sum())
            cols = cols[ok]
            # Pearson r from precomputed moments: one sparse matvec per gene
            # instead of densifying the peak block. n/(n-1) converts the
            # population covariance to the sample one.
            xty = np.asarray(A[:, cols].T @ y).ravel()
            cov = xty/n - pmean[cols]*y.mean()
            r = cov/(pstd[cols]*sy)*(n/(n-1))
            rows.append((cols, r, dist[ok]))
    allr = np.concatenate([x[1] for x in rows]); df = n-2
    tstat = allr*np.sqrt(df/np.maximum(1-allr**2, 1e-12))
    pvals = 2*st.t.sf(np.abs(tstat), df)
    thr = bh_threshold(pvals, FDR, m_extra=n_all-allr.size)
    prox_tested = int((np.concatenate([x[2] for x in rows]) <= PROMOTER).sum())
    a_=b_=0; off=0
    if thr >= 0:
        for cols, r, dist in rows:
            k=r.size; pg=pvals[off:off+k]; off+=k
            sig=np.flatnonzero(pg <= thr)
            if sig.size==0: continue
            a_ += int((dist[sig] <= PROMOTER).sum()); b_ += int((dist[sig] > PROMOTER).sum())
    def OR(c_, d_):
        """Haldane-Anscombe corrected odds ratio with a bootstrap interval.

        The +0.5 keeps the ratio finite when a cell is empty, which happens at
        the low rungs. The interval resamples the promoter-proximal count among
        detected links binomially; the background is treated as fixed, since it
        rests on millions of pairs and contributes negligible uncertainty.
        """
        o=((a_+.5)*(d_+.5))/((b_+.5)*(c_+.5))
        bs=rng.binomial(a_+b_, (a_+.5)/(a_+b_+1), 2000)
        ci=np.percentile(((bs+.5)*(d_+.5))/((a_+b_-bs+.5)*(c_+.5)), [2.5,97.5])
        return float(o), float(ci[0]), float(ci[1])
    oA=OR(prox_all, n_all-prox_all); oB=OR(prox_tested, allr.size-prox_tested)
    res=dict(cohort="primary", celltype=a.celltype, n=n, seed=a.seed,
             depth_atac=5564, depth_rna=5265,
             n_tested_all=int(n_all), prox_tested_all=int(prox_all),
             n_tested_var=int(allr.size), prox_tested_var=int(prox_tested),
             n_links=int(a_+b_), prox_links=int(a_),
             OR=oA[0], lo=oA[1], hi=oA[2],
             OR_varpeaks=oB[0], lo_varpeaks=oB[1], hi_varpeaks=oB[2])
    out=D/f"{a.celltype}_n{n}_s{a.seed}.json"; out.write_text(json.dumps(res, indent=2))
    logger.info("结果: %s", json.dumps(res))

if __name__ == "__main__":
    main()
