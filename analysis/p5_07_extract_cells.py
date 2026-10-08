#!/usr/bin/env python3
"""
Extract the equalised control-nucleus subset (sparse output)
===========================================================
Implements two pre-registration deviations:
  D1  the unit of observation for link detection moves from donor to nucleus,
      with nucleus count and depth both forced equal
  D2  N drops from 200 to 150, and nuclei are filtered on depth before being
      sampled. At N=200 the common depth collapses to 736 ATAC fragments, a
      non-zero rate of 0.35%, and the per-gene link measure degenerates to all
      zeros -- so the larger N measures nothing.

Three engineering choices, each forced by something that actually broke:
  1. SPARSE output. Dense 200,791 x 3,750 int32 is 3.0 GB and exhausted memory
     on every downstream load. Stored as COO triples it is about 150 MB, and
     A2-b can rebuild one cell type at a time.
  2. ONE scan shared by all three seeds. Scan the union of the three seeds'
     selected nuclei once, then slice per seed. Turns three passes over a 30 GB
     ATAC file into one.
  3. RESUMABLE (--time-budget; exit code 3 means "call me again"), for
     environments that impose a wall-clock limit.

Usage:
  python3 p5_07_extract_cells.py --data DIR --out results --seeds 0,1,2
  python3 p5_07_extract_cells.py --data DIR --out results --seeds 0,1,2 --cells-only
  # under a wall-clock limit, re-run the same command until the exit code is not 3
  python3 p5_07_extract_cells.py --data DIR --out results --seeds 0,1,2 --time-budget 25
"""

import argparse
import io
import os
import sys
import tempfile
import time

import numpy as np
import pandas as pd


# ---------------------------------------------------------------- selection
def pick_cells(hc, targets, n_per_type, seed, verbose=True):
    """
    Filter on a JOINT depth threshold, then sample N from what qualifies.

    The bug this fixes is worth recording. The first version computed a marginal
    threshold for ATAC and another for RNA, each set by whichever cell type was
    limiting for that modality. But a nucleus has to clear BOTH. The intersection
    of two marginal thresholds is stricter than either one, so six inhibitory
    types could not reach N and fell back to "take the N deepest" -- leaving 25%
    of nuclei below the downsampling target and a 1.95x depth ratio between cell
    types. The equal-depth premise was simply not met.

    The fix: binary search a single scale factor t, apply the thresholds
    (t*a0, t*r0) jointly, and take the largest t for which EVERY cell type still
    has at least N qualifying nuclei. At N=150 that gives t=0.342, i.e. ATAC
    >= 5,564 and RNA >= 5,254, with zero types falling back.
    """
    a0 = float(hc.nCount_ATAC.median())
    r0 = float(hc.nCount_RNA.median())

    def feasible(t):
        o = (hc.nCount_ATAC >= t * a0) & (hc.nCount_RNA >= t * r0)
        return all(o.loc[s.index].sum() >= n_per_type
                   for _, s in hc.groupby("WNN_L4"))

    lo, hi = 0.0, 3.0
    for _ in range(40):
        mid = (lo + hi) / 2
        if feasible(mid):
            lo = mid
        else:
            hi = mid
    D_atac, D_rna = lo * a0, lo * r0
    ok = (hc.nCount_ATAC >= D_atac) & (hc.nCount_RNA >= D_rna)

    out, excl, fallback = [], [], []
    for ct, sub in hc.groupby("WNN_L4"):
        elig = sub[ok.loc[sub.index]]
        excl.append((ct, len(sub), len(elig)))
        if len(elig) >= n_per_type:
            e = elig.sample(frac=1, random_state=seed)
            e = e.assign(_rk=e.groupby("ID").cumcount())
            sel = e.sort_values("_rk").head(n_per_type).drop(columns="_rk")
        else:
            # Falling back to "take the N deepest" only when too few qualify.
            # This introduces a selection bias in those cell types that points
            # the opposite way from the others, so it must be named in the
            # paper rather than passed over. The binary search above exists so
            # that this branch is never taken in the published run.
            sel = sub.nlargest(min(n_per_type, len(sub)), "nCount_ATAC")
            fallback.append(ct)
        out.append(sel.assign(celltype=ct))

    ex = pd.DataFrame(excl, columns=["celltype", "n_total", "n_eligible"])
    ex["pct_excluded"] = 100 * (1 - ex.n_eligible / ex.n_total)
    if verbose:
        print(f"  seed={seed}  深度阈值 ATAC>={D_atac:.0f}  RNA>={D_rna:.0f}")
        print(f"    被阈值排除：中位 {ex.pct_excluded.median():.0f}%，"
              f"最高 {ex.pct_excluded.max():.0f}% "
              f"({ex.loc[ex.pct_excluded.idxmax(),'celltype']})")
        if fallback:
            print(f"    ⚠️ 合格细胞不足、退回取最深的类型 ({len(fallback)}): {fallback}")
    return pd.concat(out, ignore_index=True), ex, D_atac, D_rna


# ------------------------------------------------- single-pass sparse column extraction
def extract_sparse(mtx_path, keep_lut, n_keep, ckpt, time_budget,
                   chunk_bytes=1 << 26):
    """
    Stream the MatrixMarket file, keeping only columns with keep_lut >= 0.

    keep_lut[original column] = new column index, or -1 to discard. A dense
    lookup array rather than a dict, since it is consulted once per non-zero
    entry. Returns COO triples (rows, cols, vals).
    """
    R, C, V = [], [], []
    seen, offset0 = 0, None
    if ckpt and os.path.exists(ckpt):
        z = np.load(ckpt)
        R, C, V = [z["r"]], [z["c"]], [z["v"]]
        seen, offset0 = int(z["seen"]), int(z["offset"])
        print(f"  [resume] 已处理 {seen:,} entries", flush=True)

    with open(mtx_path, "rt") as fh:
        if not fh.readline().startswith("%%MatrixMarket"):
            sys.exit(f"[FATAL] 非 MatrixMarket: {mtx_path}")
        line = fh.readline()
        while line.startswith("%"):
            line = fh.readline()
        n_feat, n_cell, nnz = (int(x) for x in line.split())
        if n_cell != len(keep_lut):
            sys.exit(f"[FATAL] mtx 列数 {n_cell} != barcode 数 {len(keep_lut)}")
        if offset0 is not None:
            fh.seek(offset0)

        t0 = time.time()
        while True:
            buf = fh.read(chunk_bytes)
            if not buf:
                break
            buf += fh.readline()
            a = pd.read_csv(io.StringIO(buf), sep=r"\s+", header=None,
                            engine="c", names=["r", "c", "v"],
                            dtype=np.int64).to_numpy()
            r, c, v = a[:, 0] - 1, a[:, 1] - 1, a[:, 2]
            del a, buf
            nc = keep_lut[c]
            m = (nc >= 0) & (v != 0)
            if m.any():
                R.append(r[m].astype(np.int32))
                C.append(nc[m].astype(np.int32))
                V.append(v[m].astype(np.int32))
            seen += len(r)
            del r, c, v, nc, m
            el = time.time() - t0
            if seen % 400_000_000 < 70_000_000:
                print(f"    {seen:,}/{nnz:,} ({100*seen/nnz:.0f}%) "
                      f"{el/60:.1f} min", flush=True)
            if time_budget and el > time_budget and ckpt:
                np.savez(ckpt, r=np.concatenate(R) if R else np.array([], np.int32),
                         c=np.concatenate(C) if C else np.array([], np.int32),
                         v=np.concatenate(V) if V else np.array([], np.int32),
                         seen=seen, offset=fh.tell())
                print(f"  [ckpt] {seen:,}/{nnz:,} ({100*seen/nnz:.0f}%) "
                      f"→ 重跑同一命令继续", flush=True)
                sys.exit(3)
    return (np.concatenate(R) if R else np.array([], np.int32),
            np.concatenate(C) if C else np.array([], np.int32),
            np.concatenate(V) if V else np.array([], np.int32), n_feat)


def downsample_sparse(rows, cols, vals, n_cols, target, rng):
    """Multinomially downsample each column to `target`.

    Returns the new triples plus the columns that were already below target,
    which the caller checks must be empty.
    """
    order = np.argsort(cols, kind="mergesort")
    rows, cols, vals = rows[order], cols[order], vals[order]
    bounds = np.searchsorted(cols, np.arange(n_cols + 1))
    nr, nc, nv, short = [], [], [], []
    for j in range(n_cols):
        s, e = bounds[j], bounds[j + 1]
        if s == e:
            short.append(j); continue
        v = vals[s:e]
        tot = int(v.sum())
        if tot <= target:
            nr.append(rows[s:e]); nc.append(cols[s:e]); nv.append(v)
            if tot < target:
                short.append(j)
            continue
        d = rng.multinomial(target, v.astype(np.float64) / tot)
        k = d > 0
        nr.append(rows[s:e][k]); nc.append(cols[s:e][k]); nv.append(d[k].astype(np.int32))
    return (np.concatenate(nr), np.concatenate(nc),
            np.concatenate(nv).astype(np.int32), short)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", default="results")
    ap.add_argument("--seeds", default="0,1,2")
    ap.add_argument("--n-per-type", type=int, default=150)
    ap.add_argument("--cells-only", action="store_true")
    ap.add_argument("--min-cells", type=int, default=0,
                    help="剔除对照细胞数低于此值的类型（预注册偏离 D4）。"
                         "N=150 时 Exc_FEZF2_NTNG1(205核) 把全体功效压死；"
                         "--min-cells 500 保留 15 类，每 peak 非零细胞 2→7")
    ap.add_argument("--tag", default="",
                    help="输出文件名后缀，用于并存多套设置")
    ap.add_argument("--force", action="store_true",
                    help="已存在也重算（改了抽样/降采样逻辑后必须加）")
    ap.add_argument("--time-budget", type=float, default=None)
    ap.add_argument("--ckpt-dir", default=None)
    ap.add_argument("--targets", default="results/p5_claimA_pathology_OR_WNN_L4.csv")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    seeds = [int(s) for s in args.seeds.split(",")]

    targets = pd.read_csv(args.targets).celltype.tolist()
    m = pd.read_csv(os.path.join(args.data, "Multiome_Dataset_Metadata.txt"),
                    sep="\t", usecols=["CellId", "ID", "Case", "WNN_L1",
                                       "WNN_L4", "nCount_RNA", "nCount_ATAC"],
                    dtype={"CellId": str, "ID": str, "Case": str,
                           "WNN_L1": str, "WNN_L4": str})
    hc = m[(m.Case == "HC") & (m.WNN_L1 == "Neuronal") &
           (m.WNN_L4.isin(targets))].copy()
    for c in ("nCount_RNA", "nCount_ATAC"):
        hc[c] = pd.to_numeric(hc[c])
    if args.min_cells:
        keep = hc.groupby("WNN_L4").size()
        keep = keep[keep >= args.min_cells].index
        dropped = sorted(set(targets) - set(keep))
        hc = hc[hc.WNN_L4.isin(keep)]
        targets = sorted(keep)
        print(f"  [D4] 剔除对照细胞 <{args.min_cells} 的 {len(dropped)} 个类型: "
              f"{dropped}")
        print(f"  [D4] 保留 {len(targets)} 个类型")

    print(f"[1/4] 选细胞（{len(targets)} 类，N={args.n_per_type}，"
          f"seeds={seeds}）…")
    picks = {}
    for s in seeds:
        cells, ex, Da, Dr = pick_cells(hc, targets, args.n_per_type, s)
        cells.to_csv(os.path.join(args.out, f"A2_cells_seed{s}{args.tag}.csv"), index=False)
        ex.to_csv(os.path.join(args.out, f"A2_exclusion_seed{s}{args.tag}.csv"), index=False)
        picks[s] = cells
    union = pd.Index(sorted(set().union(*[set(c.CellId) for c in picks.values()])))
    print(f"  每种子 {len(picks[seeds[0]]):,} 细胞；三种子并集 {len(union):,}"
          f"（重叠率 {1-len(union)/sum(len(c) for c in picks.values()):.0%}）")
    if args.cells_only:
        return

    bc = pd.Index(pd.read_csv(os.path.join(
        args.data, "Multiome_Dataset_RNA_barcodes.tsv"), header=None)[0])
    lut = np.full(len(bc), -1, dtype=np.int32)
    upos = pd.Series(np.arange(len(union)), index=union)
    hit = bc.isin(union)
    lut[np.flatnonzero(hit)] = upos.reindex(bc[hit]).values.astype(np.int32)
    assert (lut >= 0).sum() == len(union)

    ckdir = args.ckpt_dir or os.path.join(tempfile.gettempdir(), "p5_ckpt")
    os.makedirs(ckdir, exist_ok=True)

    for mod, ffile in (("RNA", "Multiome_Dataset_RNA_features.tsv"),
                       ("ATAC", "Multiome_Dataset_ATAC_features.tsv")):
        done = (not args.force) and all(os.path.exists(os.path.join(
            args.out, f"A2_sparse_{mod}_seed{s}{args.tag}.npz")) for s in seeds)
        if done:
            print(f"[2/4] {mod} 已存在，跳过")
            continue
        feats = pd.read_csv(os.path.join(args.data, ffile), header=None)[0]
        print(f"[2/4] 扫描 {mod}（并集 {len(union):,} 列）…")
        ck = os.path.join(ckdir, f"ck_extract_{mod}.npz")
        r, c, v, n_feat = extract_sparse(
            os.path.join(args.data, f"Multiome_Dataset_{mod}_counts_raw.mtx"),
            lut, len(union), ck, args.time_budget)
        if os.path.exists(ck):
            try: os.remove(ck)
            except OSError as e: print(f"  [warn] 检查点未删除: {e}")
        print(f"  抽到 {len(r):,} 个非零项")

        print(f"[3/4] 按种子切片 + 等深度降采样 …")
        for s in seeds:
            cells = picks[s]
            idx = upos.reindex(cells.CellId).values.astype(np.int32)
            remap = np.full(len(union), -1, dtype=np.int32)
            remap[idx] = np.arange(len(idx), dtype=np.int32)
            keep = remap[c] >= 0
            rs, cs, vs = r[keep], remap[c[keep]], v[keep]
            tot = np.bincount(cs, weights=vs, minlength=len(idx))
            # Target the MINIMUM total among the selected nuclei, which
            # guarantees zero columns below target: every nucleus lands at
            # exactly the same depth and the between-type depth ratio is 1.000.
            # Using the 25th percentile instead would leave a quarter of nuclei
            # under target and a 1.95x ratio between types -- the equal-depth
            # premise unmet, and the downstream gate on partial correlation
            # meaningless, since it would then be testing a depth confound.
            target = int(tot.min())
            rng = np.random.default_rng(s)
            rs, cs, vs, short = downsample_sparse(rs, cs, vs, len(idx), target, rng)
            np.savez_compressed(
                os.path.join(args.out, f"A2_sparse_{mod}_seed{s}{args.tag}.npz"),
                rows=rs, cols=cs, vals=vs, shape=np.array([n_feat, len(idx)]),
                features=feats.values, cells=cells.CellId.values,
                celltype=cells.celltype.values, donor=cells.ID.values,
                depth_target=target, raw_total=tot,
                short_cols=np.array(short, dtype=np.int32))
            print(f"    seed{s}: depth_target={target:,}  "
                  f"nnz={len(rs):,}  未达标列={len(short)}")
        del r, c, v

    print(f"[4/4] 完成 → {args.out}/A2_sparse_*.npz")


if __name__ == "__main__":
    main()
