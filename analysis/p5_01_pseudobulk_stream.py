#!/usr/bin/env python3
"""
Streaming pseudobulk over the 30 GB ATAC and 9.4 GB RNA matrices
================================================================
One sequential pass over each MatrixMarket file, accumulating by
(donor x cell type), without ever holding the full sparse matrix in memory.

Memory is features x groups x 8 bytes:
    ATAC  200,791 peaks x ~948 groups = 1.5 GB
    RNA    35,367 genes x ~948 groups = 0.27 GB
plus 2-3 GB of chunk buffers. Under 6 GB in total, which an ordinary
workstation can run.

The RNA matrix carries about 13% EXPLICIT ZERO entries, an artefact of how it
was written. They are discarded here and the count is logged: if that fraction
does not match expectation, the file format has been misread and the run should
stop rather than produce a quietly wrong answer.

The RNA and ATAC barcodes.tsv files were verified byte-identical, so one
cell-to-group mapping serves both column indices. The script re-checks this
anyway and aborts on any mismatch.

Usage:
  # smoke test: first 5 million lines, about 20 seconds
  python3 p5_01_pseudobulk_stream.py --data DIR --out out/ --modality RNA --max-lines 5000000

  # full run
  python3 p5_01_pseudobulk_stream.py --data DIR --out out/ --modality RNA
  python3 p5_01_pseudobulk_stream.py --data DIR --out out/ --modality ATAC

Outputs:
  out/pseudobulk_{RNA,ATAC}_{level}.npz    counts (features x groups), with group
                                           and feature names
  out/pseudobulk_{...}_groupinfo.csv       donor, cell type, nuclei and total
                                           counts per group
"""

import argparse
import gzip
import io
import os
import sys
import tempfile
import time

import numpy as np
import pandas as pd

CHUNK_ROWS = 20_000_000          # 每次读入的 mtx 行数


def load_cell_groups(data_dir, level, modality, group_extra=None):
    """Returns (group_of_col[int32, n_cells], group_names[list], donors, celltypes)."""
    # Use pandas Index for every set operation here. np.setdiff1d and
    # np.unique fall onto a pathologically slow path for object-dtype string
    # arrays -- slow enough to look like a hang.
    bc_file = os.path.join(data_dir, f"Multiome_Dataset_{modality}_barcodes.tsv")
    barcodes = pd.Index(pd.read_csv(bc_file, header=None)[0])

    # Cross-check: the two modalities' barcodes must match exactly, or the
    # column indices cannot share one mapping.
    other = "ATAC" if modality == "RNA" else "RNA"
    other_file = os.path.join(data_dir, f"Multiome_Dataset_{other}_barcodes.tsv")
    if os.path.exists(other_file):
        ob = pd.Index(pd.read_csv(other_file, header=None)[0])
        if not barcodes.equals(ob):
            sys.exit("[FATAL] RNA 与 ATAC 的 barcode 顺序不一致 —— "
                     "不能假设同核配对。停止。")

    meta = pd.read_csv(os.path.join(data_dir, "Multiome_Dataset_Metadata.txt"),
                       sep="\t",
                       usecols=["CellId", "ID", level] +
                               ([group_extra] if group_extra else []),
                       dtype=str)
    meta = meta.set_index("CellId")
    if not meta.index.is_unique:
        sys.exit("[FATAL] 元数据 CellId 有重复。停止。")

    missing = barcodes.difference(meta.index)
    if len(missing):
        sys.exit(f"[FATAL] {len(missing)} 个 barcode 在元数据中找不到。停止。")

    meta = meta.reindex(barcodes)
    # Group key: donor || [extra keys] || cell type
    key = meta["ID"].astype(str)
    if group_extra:
        key = key + "||" + meta[group_extra].astype(str)
    key = key + "||" + meta[level].astype(str)
    codes, uniques = pd.factorize(key, sort=True)
    return (codes.astype(np.int32), list(uniques),
            meta["ID"].to_numpy(), meta[level].to_numpy())


def stream_pseudobulk(mtx_path, group_of_col, n_groups, max_lines=None,
                      drop_zeros=True, log_every=100_000_000,
                      group_subset=None, chunk_rows=CHUNK_ROWS,
                      ckpt_path=None, time_budget=None):
    """
    Scan the mtx sequentially, returning a (features x groups) int32 matrix.

    The memory choices here were all forced by an actual out-of-memory kill on
    a 3 GB machine:

      * the accumulator is INT32, not int64. ATAC totals about 2e6 counts per
        group, far below 2^31.
      * each chunk merges duplicate (feature, group) keys with np.unique before
        scattering. Temporaries are then CHUNK-sized rather than
        (n_feat x n_groups)-sized. Calling np.bincount over the full length
        instead allocates a float64 temporary as large as the accumulator on
        every chunk, doubling peak memory -- which is exactly what the OOM
        killer caught. Pure numpy, no scipy dependency.
      * group_subset accumulates only part of the groups, trading extra passes
        over the file for lower memory.

    ATAC at full size: 200,791 x 935 x 4 B = 0.75 GB accumulator, a dense
    temporary of the same size, and about 0.5 GB of chunk buffers -- roughly
    2.0-2.5 GB. Below 4 GB of RAM, split into N passes with --group-chunks N.
    """
    opener = gzip.open if mtx_path.endswith(".gz") else open
    with opener(mtx_path, "rt") as fh:
        header = fh.readline()
        if not header.startswith("%%MatrixMarket"):
            sys.exit(f"[FATAL] 不是 MatrixMarket 文件: {mtx_path}")
        line = fh.readline()
        while line.startswith("%"):
            line = fh.readline()
        n_feat, n_cell, nnz = (int(x) for x in line.split())
        print(f"  header: {n_feat} features × {n_cell} cells, {nnz:,} stored entries")
        if n_cell != len(group_of_col):
            sys.exit(f"[FATAL] mtx 列数 {n_cell} ≠ barcode 数 {len(group_of_col)}")

        # group_subset: accumulate only these global group ids this pass
        if group_subset is None:
            keep_mask = None
            n_out_groups = n_groups
            local_of_global = None
        else:
            keep_mask = np.zeros(n_groups, dtype=bool)
            keep_mask[group_subset] = True
            local_of_global = np.full(n_groups, -1, dtype=np.int32)
            local_of_global[group_subset] = np.arange(len(group_subset),
                                                      dtype=np.int32)
            n_out_groups = len(group_subset)

        data_start = fh.tell()

        # ---- Resumable -------------------------------------------------
        # Read BYTE blocks by hand and hand them to pandas' C parser, rather
        # than iterating pd.read_csv(chunksize=). That keeps fh.tell() reliable,
        # so a byte offset can go into the checkpoint and the next run can seek
        # back to it. Iterating read_csv buffers internally and tell() is then
        # meaningless.
        out = np.zeros((n_feat, n_out_groups), dtype=np.int32)
        seen = 0
        n_zero_dropped = 0
        if ckpt_path and os.path.exists(ckpt_path):
            z = np.load(ckpt_path)
            if int(z["n_out_groups"]) != n_out_groups or int(z["n_feat"]) != n_feat:
                sys.exit("[FATAL] 检查点维度与当前参数不符。删掉它或换 --out。")
            out = z["out"]
            seen = int(z["seen"]); n_zero_dropped = int(z["n_zero_dropped"])
            fh.seek(int(z["offset"]))
            print(f"  [resume] 从检查点恢复：已处理 {seen:,} entries "
                  f"({100*seen/nnz:.1f}%)", flush=True)

        t0 = time.time()
        block_bytes = max(chunk_rows * 12, 1 << 20)   # 每行约 12 字节

        def blocks():
            while True:
                buf = fh.read(block_bytes)
                if not buf:
                    return
                tail = fh.readline()       # 补齐被截断的最后一行
                yield buf + tail

        for text in blocks():
            chunk = pd.read_csv(io.StringIO(text), sep=r"\s+", header=None,
                                engine="c", names=["r", "c", "v"],
                                dtype=np.int64)
            a = chunk.to_numpy()
            del chunk
            if max_lines is not None and seen + len(a) > max_lines:
                a = a[: max(max_lines - seen, 0)]
            seen_this = len(a)
            r, c, v = a[:, 0] - 1, a[:, 1] - 1, a[:, 2]
            del a
            if drop_zeros:
                nz = v != 0
                n_zero_dropped += int((~nz).sum())
                r, c, v = r[nz], c[nz], v[nz]
            g = group_of_col[c]
            if keep_mask is not None:
                sel = keep_mask[g]
                r, v, g = r[sel], v[sel], local_of_global[g[sel]]
            # Merge duplicate (feature, group) keys within the chunk, then
            # scatter. Temporaries stay chunk-sized, independent of the
            # accumulator.
            flat = r.astype(np.int64) * n_out_groups + g.astype(np.int64)
            uniq, inv = np.unique(flat, return_inverse=True)
            sums = np.bincount(inv, weights=v.astype(np.float64),
                               minlength=len(uniq))
            out_flat = out.reshape(-1)
            out_flat[uniq] += sums.astype(np.int32)
            del flat, uniq, inv, sums, r, c, v, g
            seen += seen_this
            el = time.time() - t0
            if seen % log_every < chunk_rows:
                print(f"    {seen:,}/{nnz:,} entries  "
                      f"({100*seen/nnz:5.1f}%)  {el/60:.1f} min", flush=True)
            if max_lines is not None and seen >= max_lines:
                break
            # Budget exhausted: checkpoint and exit 3, meaning "call me again"
            if time_budget and el > time_budget and ckpt_path:
                np.savez(ckpt_path, out=out, seen=seen, offset=fh.tell(),
                         n_zero_dropped=n_zero_dropped,
                         n_out_groups=n_out_groups, n_feat=n_feat)
                print(f"  [ckpt] 已存检查点：{seen:,}/{nnz:,} "
                      f"({100*seen/nnz:.1f}%) → 再次运行同一命令即可继续",
                      flush=True)
                sys.exit(3)

    # Nuclei per group, computed from the mapping alone -- no matrix needed
    n_cells_per_group = np.bincount(group_of_col, minlength=n_groups)

    print(f"  扫描完成: {seen:,} entries, 丢弃显式 0 共 {n_zero_dropped:,} "
          f"({100*n_zero_dropped/max(seen,1):.1f}%), 用时 {(time.time()-t0)/60:.1f} min")
    if drop_zeros and seen > 1_000_000:
        frac = n_zero_dropped / seen
        if frac > 0.30:
            print(f"  ⚠️  显式 0 占比 {frac:.1%} 异常偏高，请核对文件格式再往下走。")
    return out, n_cells_per_group


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", default="out")
    ap.add_argument("--modality", choices=["RNA", "ATAC"], required=True)
    ap.add_argument("--level", default="WNN_L2.5",
                    help="细胞类型层级列名 (WNN_L1/L1.5/L2/L2.5/L3/L4)")
    ap.add_argument("--max-lines", type=int, default=None,
                    help="冒烟测试用：只读前 N 行")
    ap.add_argument("--group-chunks", type=int, default=1,
                    help="内存不足时把组分成 N 趟累加（代价是扫描 N 遍文件）")
    ap.add_argument("--chunk-rows", type=int, default=CHUNK_ROWS,
                    help="每次读入的 mtx 行数；内存紧张就调小")
    ap.add_argument("--group-extra", default=None,
                    help="在 供体||细胞类型 之外再加一个分组键（如 Well），"
                         "用于跨芯片技术重复标定（步骤 0.4 的 L2 层）")
    ap.add_argument("--ckpt-dir", default=None,
                    help="检查点目录。默认放系统临时目录，**不要放 iCloud/网络盘**："
                         "检查点可达数百 MB，云同步会拖垮写入甚至写不完")
    ap.add_argument("--time-budget", type=float, default=None,
                    help="单次运行秒数上限；到点存检查点并以退出码 3 退出，"
                         "重跑同一命令自动续跑（适合有超时限制的环境）")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    print(f"[1/3] 构建 cell→group 映射 (level={args.level}) …")
    gcol, gnames, donors, ctypes = load_cell_groups(args.data, args.level,
                                                    args.modality,
                                                    args.group_extra)
    print(f"  {len(gnames)} 组 = {len(set(donors))} 供体 × {len(set(ctypes))} 细胞类型")

    feat_file = os.path.join(args.data,
                             f"Multiome_Dataset_{args.modality}_features.tsv")
    features = pd.read_csv(feat_file, header=None)[0].to_numpy()

    print(f"[2/3] 流式扫描 {args.modality} 矩阵 "
          f"({args.group_chunks} 趟) …")
    mtx = os.path.join(args.data, f"Multiome_Dataset_{args.modality}_counts_raw.mtx")
    if args.group_chunks <= 1:
        ckdir = args.ckpt_dir or os.path.join(
            tempfile.gettempdir(), "p5_ckpt")
        os.makedirs(ckdir, exist_ok=True)
        ck = os.path.join(ckdir, f"ckpt_{args.modality}_"
                          f"{args.level.replace('.','')}.npz")
        mat, _ = stream_pseudobulk(mtx, gcol, len(gnames),
                                   max_lines=args.max_lines,
                                   chunk_rows=args.chunk_rows,
                                   ckpt_path=ck, time_budget=args.time_budget)
        # Clean up the checkpoint. On iCloud or a network drive os.remove can
        # raise PermissionError -- at which point the results are already
        # computed and written. Failing the whole job over a failed cleanup
        # would invite someone to re-run a ten-minute scan for nothing.
        if os.path.exists(ck):
            try:
                os.remove(ck)
            except OSError as e:
                print(f"  [warn] 检查点未能删除（{e}）。结果不受影响，"
                      f"下次重跑前请手动删除：{ck}")
    else:
        parts, order = [], []
        for k, sub in enumerate(np.array_split(np.arange(len(gnames)),
                                               args.group_chunks)):
            print(f"  -- 第 {k+1}/{args.group_chunks} 趟，{len(sub)} 组")
            p, _ = stream_pseudobulk(mtx, gcol, len(gnames),
                                     max_lines=args.max_lines,
                                     group_subset=sub,
                                     chunk_rows=args.chunk_rows)
            parts.append(p); order.append(sub)
        mat = np.concatenate(parts, axis=1)
        mat = mat[:, np.argsort(np.concatenate(order))]
        del parts
    ncells = np.bincount(gcol, minlength=len(gnames))
    if mat.shape[0] != len(features):
        sys.exit(f"[FATAL] 特征数不符: mtx {mat.shape[0]} vs tsv {len(features)}")

    print("[3/3] 写出 …")
    tag = f"{args.modality}_{args.level.replace('.', '')}"
    if args.group_extra:
        tag += f"_by{args.group_extra}"
    if args.max_lines:
        tag += "_SMOKE"
    np.savez_compressed(os.path.join(args.out, f"pseudobulk_{tag}.npz"),
                        counts=mat, groups=np.array(gnames),
                        features=features, n_cells=ncells)

    info = pd.DataFrame({
        "group": gnames,
        "donor": [g.split("||")[0] for g in gnames],
        "celltype": [g.split("||")[-1] for g in gnames],
        "n_cells": ncells,
        "total_counts": mat.sum(0),
    })
    # total_counts is what every later power-balancing step rests on. The
    # source paper concedes that DEG counts correlate strongly with read counts,
    # so any comparison of effect size across cell types has to match or
    # downsample on this column first.
    info.to_csv(os.path.join(args.out, f"pseudobulk_{tag}_groupinfo.csv"),
                index=False)
    print(f"[done] {args.out}/pseudobulk_{tag}.npz  "
          f"shape={mat.shape}  总计数={mat.sum():,}")


if __name__ == "__main__":
    main()
