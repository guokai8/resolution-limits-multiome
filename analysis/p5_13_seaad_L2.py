#!/usr/bin/env python3
"""
The expression floor in SEA-AD, computed exactly as in the primary cohort
========================================================================
p5_12 --mode l1 can work from metadata alone; the expression layer cannot,
because it is defined across genes. This script builds (donor x library x cell
type) pseudobulk by STREAMING the 31 GB h5ad, then applies the identical
statistic:

    floor = median over genes of |log2 CPM_lib1 - log2 CPM_lib2|
    n_eff = harmonic mean of this cell type's nuclei in the two libraries

Why the metadata's total-UMI proxy will not do, demonstrated on the primary
cohort's own data: the true exponent is -0.507 and the proxy returns -0.248 --
enough to turn a replication into a non-replication.

Memory strategy for a 3 GB machine:
  * float32 accumulator of shape (n_group x n_gene). 80 donors x 2 libraries x
    24 types is about 3,840 groups; at 36k genes that is roughly 560 MB. If
    that still exhausts memory, --donor-chunks K splits it into K passes, each
    accumulating some donors, at the cost of reading the file K times.
  * read the CSR data and indices block by block; X is never loaded whole.

Usage:
  # 1) inspect first (seconds, and mandatory): is X raw counts or normalised?
  python3 p5_13_seaad_L2.py --h5ad SEAAD_...h5ad --inspect

  # 2) then run, choosing --layer from what the inspection reported
  python3 p5_13_seaad_L2.py --h5ad SEAAD_...h5ad --layer raw --out results
"""

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

try:
    import h5py
except ImportError:
    sys.exit("需要 h5py：/opt/homebrew/bin/python -m pip install h5py")

DONOR, BATCH, CTYPE = "Donor ID", "library_prep", "Subclass"


# --------------------------------------------------------------- reading obs
def read_obs(f, cols):
    """Read columns from an h5ad /obs, decoding categoricals correctly."""
    g = f["obs"]
    out = {}
    for c in cols:
        if c not in g:
            continue
        n = g[c]
        if isinstance(n, h5py.Group) and "codes" in n:      # categorical
            codes = n["codes"][:]
            cats = n["categories"][:]
            cats = np.array([x.decode() if isinstance(x, bytes) else x
                             for x in cats], dtype=object)
            v = np.where(codes >= 0, cats[np.clip(codes, 0, None)], None)
        else:
            v = n[:]
            if v.dtype.kind == "S":
                v = np.array([x.decode() for x in v], dtype=object)
        out[c] = v
    return pd.DataFrame(out)


def inspect(path):
    with h5py.File(path, "r") as f:
        print("顶层:", list(f.keys()))
        for key in ("X", "raw/X", "layers"):
            if key in f:
                n = f[key]
                if isinstance(n, h5py.Group) and "data" in n:
                    d = n["data"]
                    enc = n.attrs.get("encoding-type", b"?")
                    print(f"\n{key}: encoding={enc}  nnz={d.shape[0]:,}  "
                          f"dtype={d.dtype}  shape={n.attrs.get('shape')}")
                    s = d[:200000]
                    intish = np.allclose(s, np.round(s))
                    print(f"   前 20 万个值：min={s.min():.4g} max={s.max():.4g} "
                          f"是否整数={intish}  "
                          f"{'✅ 像原始计数' if intish else '⚠️ 已归一化，不能直接用'}")
                elif isinstance(n, h5py.Group):
                    print(f"\n{key}: 子项 {list(n.keys())}")
        print("\nobs 列:", list(f["obs"].keys())[:60])
        if "var" in f:
            print("var 列:", list(f["var"].keys())[:20])


# --------------------------------------------------------- streaming pseudobulk
def build_pseudobulk(path, xkey, gid, n_group, n_gene, chunk=8_000_000,
                     ckpt=None, ckpt_min=20.0):
    """
    One streaming pass over the CSR, accumulating (group x gene) counts.

    Cells with gid < 0 are skipped.

    SEA-AD has 6.6 billion non-zeros, about 79 GB decompressed, so a single
    pass takes one to two hours and must be resumable. Following the lesson
    from p5_01, the checkpoint is written to the system temp directory -- not
    to iCloud, where a 750 MB checkpoint stalls.
    """
    acc = np.zeros((n_group, n_gene), dtype=np.float32)
    start0 = 0
    if ckpt and os.path.exists(ckpt):
        z = np.load(ckpt)
        if int(z["n_group"]) == n_group and int(z["n_gene"]) == n_gene:
            acc, start0 = z["acc"], int(z["pos"])
            print(f"  ↻ 从断点续跑：已完成 {start0:,} 个非零", flush=True)
        else:
            print("  ⚠️ 断点形状不符，忽略并重跑")
    t_ck = time.time()
    with h5py.File(path, "r") as f:
        X = f[xkey]
        indptr, data, indices = X["indptr"], X["data"], X["indices"]
        n_cell = indptr.shape[0] - 1
        ip = indptr[:]                      # (n_cell+1,) int64，约 10 MB
        nnz = int(ip[-1])
        print(f"  {n_cell:,} 细胞 / {nnz:,} 非零 / {n_gene:,} 基因", flush=True)
        t0, start = time.time(), start0
        while start < nnz:
            stop = min(start + chunk, nnz)
            # The cell rows this nnz span covers
            lo = int(np.searchsorted(ip, start, "right")) - 1
            hi = int(np.searchsorted(ip, stop, "left"))
            stop = int(ip[hi]) if hi <= n_cell else nnz    # 对齐到行边界
            if stop <= start:
                stop = min(start + chunk, nnz); hi = n_cell
            d = data[start:stop]
            j = indices[start:stop]
            rows = np.repeat(np.arange(lo, hi),
                             np.diff(ip[lo:hi + 1]).astype(np.int64))
            g = gid[rows]
            m = g >= 0
            if m.any():
                # Merge duplicate (group, gene) keys before scattering, so no
                # temporary the size of the accumulator is ever allocated
                key = g[m].astype(np.int64) * n_gene + j[m].astype(np.int64)
                uk, inv = np.unique(key, return_inverse=True)
                s = np.bincount(inv, weights=d[m].astype(np.float64))
                acc.reshape(-1)[uk] += s.astype(np.float32)
            start = stop
            if hi >= n_cell:
                break
            el = (time.time() - t0) / 60
            done = start - start0
            eta = el / max(done, 1) * (nnz - start)
            print(f"    {start:,}/{nnz:,} ({100*start/nnz:.1f}%) "
                  f"{el:.1f} min，剩余约 {eta:.0f} min", flush=True)
            if ckpt and (time.time() - t_ck) / 60 >= ckpt_min:
                np.savez(ckpt + ".tmp", acc=acc, pos=start,
                         n_group=n_group, n_gene=n_gene)
                os.replace(ckpt + ".tmp", ckpt)
                t_ck = time.time()
                print(f"    💾 断点已存 {ckpt}", flush=True)
    if ckpt and os.path.exists(ckpt):
        try:
            os.remove(ckpt)          # p5_01 教训：删断点失败不能中断结果写出
        except OSError as e:
            print(f"    (断点删除失败，无妨: {e})")
    return acc


# ------------------------------------------------------------------ L2
def compute_L2(acc, meta, n_group, min_cpm):
    """
    Compute the expression floor from accumulated counts.

    acc is (group x gene) counts; meta carries donor, batch, cell type and
    nucleus count per group.

    The gene filter must match the primary cohort's exactly:

        keep = (x > 0) | (y > 0)      # union, not intersection

    The first version here wrote (cpm1 >= 1) & (cpm2 >= 1), an intersection
    with a threshold. That removes precisely the genes present in one library
    and absent from the other -- the genes with the largest |delta log2 CPM|.
    The floor was pushed down systematically (a = 1.40 against the primary
    cohort's 4.79), and because the fraction removed varies with n, the slope
    was contaminated too (b = -0.371). That was not a failed replication; it
    was a different statistic.

    This was the fifth instance of the same error in this project, and the rule
    it enforces: a floor can only be compared with a floor computed the same
    way.
    """
    tot = acc.sum(1, keepdims=True)
    ok = tot[:, 0] > 0
    cpm = np.zeros_like(acc)
    cpm[ok] = acc[ok] / tot[ok] * 1e6
    lg = np.log2(cpm + 1.0)

    rows = []
    for (don, ct), s in meta[meta.gid >= 0].groupby(["donor", "ct"]):
        if len(s) != 2:
            continue
        (i1, n1), (i2, n2) = [(r.gid, r.n_cells) for r in s.itertuples()]
        if min(n1, n2) < 5 or not (ok[i1] and ok[i2]):
            continue
        keep = (cpm[i1] > min_cpm) | (cpm[i2] > min_cpm)   # 与 p5_04 一致：或
        if keep.sum() < 200:
            continue
        fl = float(np.median(np.abs(lg[i1][keep] - lg[i2][keep])))
        if fl > 0:
            rows.append(dict(donor=don, ct=ct, floor=fl,
                             n_eff=2 / (1 / n1 + 1 / n2), n_genes=int(keep.sum())))
    return pd.DataFrame(rows)


def fit(x, y, groups, n_boot=1500, seed=0):
    lx, ly = np.log(x), np.log(y)
    A = np.vstack([lx, np.ones_like(lx)]).T
    beta = np.linalg.lstsq(A, ly, rcond=None)[0]
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
    return beta, np.percentile(bs, [2.5, 97.5], axis=0), bs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--inspect", action="store_true")
    ap.add_argument("--layer", default="layers/UMIs",
                    help="SEA-AD 已确认：X 是归一化值，原始计数在 layers/UMIs")
    ap.add_argument("--ckpt-dir", default="/tmp",
                    help="断点目录。⚠️ 不要指向 iCloud")
    ap.add_argument("--chunk", type=int, default=8_000_000)
    ap.add_argument("--out", default="results")
    ap.add_argument("--min-cells", type=int, default=50, help="每文库最少核数")
    ap.add_argument("--min-cpm", type=float, default=0.0,
                    help="与 p5_04 一致：`(x>min) | (y>min)`，默认 0 即'任一非零'")
    ap.add_argument("--reuse-pseudobulk", action="store_true",
                    help="复用已存的伪批量，跳过 1–2 小时的流式扫描")
    ap.add_argument("--n-boot", type=int, default=1500)
    args = ap.parse_args()

    if args.inspect:
        inspect(args.h5ad); return
    os.makedirs(args.out, exist_ok=True)

    with h5py.File(args.h5ad, "r") as f:
        xkey = args.layer
        if xkey == "auto":
            xkey = "layers/UMIs" if "layers/UMIs" in f else "X"
        if xkey not in f:
            sys.exit(f"[FATAL] {xkey} 不存在。先跑 --inspect")
        # Integrality check: computing CPM from already-normalised values
        # yields a wrong floor, silently
        s = f[xkey]["data"][:200000]
        if not np.allclose(s, np.round(s)):
            sys.exit(f"[FATAL] {xkey} 不是整数计数（已归一化）。"
                     "用 --layer 指定原始计数所在处；先跑 --inspect 看清楚。")
        n_gene = int(f[xkey].attrs["shape"][1]) if "shape" in f[xkey].attrs \
            else int(f["var/_index"].shape[0])
        obs = read_obs(f, [DONOR, BATCH, CTYPE])
        print(f"[1/3] obs {len(obs):,} 行，{xkey} 为原始计数 ✅")

    obs.columns = ["donor", "batch", "ct"][:len(obs.columns)]
    obs = obs.dropna()
    # Keep only (donor, library) pairs with at least min_cells nuclei, from
    # donors having at least two such libraries
    lib = obs.groupby(["donor", "batch"]).size()
    lib = lib[lib >= args.min_cells]
    nlib = lib.reset_index().groupby("donor").batch.nunique()
    keep_don = set(nlib[nlib >= 2].index)
    keep_lib = set(map(tuple, lib.reset_index()[["donor", "batch"]].values))
    print(f"  {len(keep_don)} 个供体有 ≥2 个 ≥{args.min_cells} 核的文库")

    key = list(zip(obs.donor, obs.batch, obs.ct))
    valid = np.array([(d in keep_don) and ((d, b) in keep_lib)
                      for d, b, c in key])
    uk = pd.unique(pd.Series([k for k, v in zip(key, valid) if v]))
    kmap = {k: i for i, k in enumerate(uk)}
    gid = np.array([kmap.get(k, -1) if v else -1
                    for k, v in zip(key, valid)], dtype=np.int32)
    n_group = len(uk)
    cnt = pd.Series(gid[gid >= 0]).value_counts()
    meta = pd.DataFrame([dict(gid=i, donor=k[0], batch=k[1], ct=k[2],
                              n_cells=int(cnt.get(i, 0))) for k, i in kmap.items()])
    print(f"  {n_group:,} 个 (供体×文库×类型) 组")

    print(f"[2/3] 流式构建伪批量（累加器 {n_group*n_gene*4/1e9:.2f} GB）…")
    pb = os.path.join(args.out, "seaad_pseudobulk_L2.npz")
    if args.reuse_pseudobulk and os.path.exists(pb):
        z = np.load(pb, allow_pickle=True)
        acc = z["acc"]
        meta = pd.DataFrame(z["meta"], columns=z["meta_cols"])
        for c in ("gid", "n_cells"):
            meta[c] = meta[c].astype(int)
        print(f"  ↻ 复用已存伪批量 {pb}（跳过流式扫描）")
    else:
        ck = os.path.join(args.ckpt_dir, "p5_13_seaad_L2.ckpt.npz")
        acc = build_pseudobulk(args.h5ad, xkey, gid, n_group, n_gene,
                               chunk=args.chunk, ckpt=ck)
        # Store the pseudobulk, so changing the gene-filtering convention does
        # not mean another 1-2 hour scan
        np.savez_compressed(pb, acc=acc, meta=meta.values,
                            meta_cols=np.array(meta.columns, dtype=object))
        print(f"  💾 伪批量已存 {pb}（下次加 --reuse-pseudobulk 秒出）")

    print("[3/3] 计算 L2 下限…")
    e = compute_L2(acc, meta, n_group, args.min_cpm)
    del acc
    e.to_csv(os.path.join(args.out, "seaad_L2_pairs_genelevel.csv"), index=False)
    print(f"  {len(e):,} 个 (供体×类型) 配对 / {e.donor.nunique()} 供体 / "
          f"{e.ct.nunique()} 类型；每对用基因数中位 {e.n_genes.median():.0f}")
    print(f"  n_eff 跨度 {e.n_eff.min():.0f}–{e.n_eff.max():.0f} "
          f"({e.n_eff.max()/e.n_eff.min():.0f}×)"
          f"{'  ⚠️ <20× ，斜率估计不稳' if e.n_eff.max()/e.n_eff.min()<20 else ''}")

    b, ci, bs = fit(e.n_eff.values, e.floor.values, e.donor.values, args.n_boot)
    print("\n" + "=" * 66)
    print("SEA-AD · L2 表达层（基因层，与 Ulm 同口径）")
    print("=" * 66)
    print(f"  b = {b[0]:+.3f}  95% CI [{ci[0,0]:+.3f}, {ci[1,0]:+.3f}]"
          f"    （Ulm: −0.507 [−0.569, −0.441]）")
    inside = ci[0, 0] <= -0.5 <= ci[1, 0]
    print(f"  与纯抽样 (−0.5) 一致: {'✅ 是' if inside else '❌ 否'}")
    print(f"  系数 a = {np.exp(b[1]):.2f}    （Ulm: 4.79）")
    print("\n判定：")
    if inside:
        print("  ⭐ **L2 复现** —— 表达层标度律跨队列稳定。")
        print("     可写：三层里唯一能用公式外推的是表达层；")
        print("     组成层的噪声性质随队列而变，调控层在可及细胞数下不可行。")
    else:
        print("  ⚠️ **L2 也不复现** —— 则结论收窄为：三层下限均不可迁移，")
        print("     必须逐数据集测量。这仍可发表，但公式不能外推。")


if __name__ == "__main__":
    main()
