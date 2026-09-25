#!/usr/bin/env python3
"""
P5 基础对象 · 30 GB ATAC / 9.4 GB RNA 矩阵的流式伪批量
=======================================================
单次顺序扫描 MatrixMarket 文件，按 (供体 × 细胞类型) 分组累加，
不把完整稀疏矩阵读进内存。

内存占用 = features × groups × 8 bytes：
    ATAC  200,791 peaks × ~948 groups ≈ 1.5 GB
    RNA    35,367 genes × ~948 groups ≈ 0.27 GB
外加分块缓冲约 2–3 GB。整体 < 6 GB，普通工作站可跑。

⚠️  RNA 矩阵含约 13% 的**显式 0** 条目（writeMM 的产物）。本脚本会丢弃它们，
    并在日志中报告丢弃数——若这个比例与预期不符，说明文件格式理解有误，必须停下。

⚠️  RNA 与 ATAC 的 barcodes.tsv 已核实逐字节相同，因此列索引可共用同一份
    cell→group 映射。脚本仍会重新校验，不一致即中止。

用法：
  # 冒烟测试（只读前 500 万行，约 20 秒）
  python3 p5_01_pseudobulk_stream.py --data DIR --out out/ --modality RNA --max-lines 5000000

  # 正式运行
  python3 p5_01_pseudobulk_stream.py --data DIR --out out/ --modality RNA
  python3 p5_01_pseudobulk_stream.py --data DIR --out out/ --modality ATAC

输出：
  out/pseudobulk_{RNA,ATAC}_{level}.npz    counts (features × groups), 附 group/feature 名
  out/pseudobulk_{...}_groupinfo.csv       每组的 供体/细胞类型/细胞数/总计数
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
    """返回 (group_of_col[int32, n_cells], group_names[list], donors, celltypes)"""
    # ⚠️ 全程用 pandas Index 做集合运算。不要用 np.setdiff1d / np.unique：
    #    在 object-dtype 的字符串数组上它们会退化到极慢的路径（实测挂死）。
    bc_file = os.path.join(data_dir, f"Multiome_Dataset_{modality}_barcodes.tsv")
    barcodes = pd.Index(pd.read_csv(bc_file, header=None)[0])

    # 交叉校验：两个模态的 barcode 必须完全一致，否则列索引不可共用
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
    # 分组键：供体 || [额外键] || 细胞类型
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
    顺序扫描 mtx，返回 (features × groups) int32 计数矩阵。

    内存要点（血泪教训：3 GB 机器上 np.bincount 路线直接 OOM 被杀）：
      * 累加器用 **int32** 而非 int64 —— ATAC 每组总计数约 2e6，远低于 2^31
      * 每个分块先 np.unique 合并重复的 (feature, group) 键，再散射累加。
        临时数组只有**分块大小**，而不是 (n_feat × n_groups) 全长。
        对整个全长调用 np.bincount 会每块分配一个与累加器等大的 float64
        临时数组，直接翻倍内存 —— 实测在 3 GB 机器上被 OOM killer 干掉。
        （纯 numpy 实现，不依赖 scipy）
      * group_subset 允许一次只累加一部分组，用多趟扫描换内存

    ATAC 全量所需内存：200,791 × 935 × 4 B ≈ 0.75 GB（累加器）
                        + 同等大小的稠密临时 ≈ 0.75 GB
                        + 分块缓冲 ≈ 0.5 GB   →  约 2.0–2.5 GB
    若机器内存 < 4 GB，用 --group-chunks N 分 N 趟。
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

        # group_subset: 本趟只累加这些全局组编号；其余行丢弃
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

        # ---- 断点续跑 ----------------------------------------------------
        # 自己按**字节块**读，再交给 pandas 的 C 解析器。这样 fh.tell() 是可靠的，
        # 可以把字节偏移量存进检查点，下次 seek 回来继续。
        # （直接用 pd.read_csv(chunksize=) 迭代时 tell() 因缓冲而不可靠。）
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
            # 合并分块内重复的 (feature, group) 键后散射累加。
            # 临时数组均为分块大小，与累加器规模无关。
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
            # 时间预算用尽 → 存检查点后退出（退出码 3 = "请再调用我一次"）
            if time_budget and el > time_budget and ckpt_path:
                np.savez(ckpt_path, out=out, seen=seen, offset=fh.tell(),
                         n_zero_dropped=n_zero_dropped,
                         n_out_groups=n_out_groups, n_feat=n_feat)
                print(f"  [ckpt] 已存检查点：{seen:,}/{nnz:,} "
                      f"({100*seen/nnz:.1f}%) → 再次运行同一命令即可继续",
                      flush=True)
                sys.exit(3)

    # 每组细胞数（与矩阵无关，直接从映射算）
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
        # 完成后清理检查点。⚠️ 在 iCloud / 网络盘上 os.remove 可能抛
        # PermissionError —— 那时结果其实已经算完写好了，绝不能因为清理失败
        # 而让整个任务以非零码退出（否则会诱使人重跑 10 分钟的扫描）。
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
    # ⚠️ total_counts 是后续所有"功效均衡"的依据：论文自己承认 DEG 数与
    #    reads 数强相关，任何跨细胞类型的效应量比较都必须先在这一列上匹配或降采样。
    info.to_csv(os.path.join(args.out, f"pseudobulk_{tag}_groupinfo.csv"),
                index=False)
    print(f"[done] {args.out}/pseudobulk_{tag}.npz  "
          f"shape={mat.shape}  总计数={mat.sum():,}")


if __name__ == "__main__":
    main()
