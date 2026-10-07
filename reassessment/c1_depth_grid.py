#!/usr/bin/env python3
"""审稿意见 V3：核数 × 测序深度的二维响应面。

正文只做了两条线——固定深度变核数，和固定 150 核变深度——所以"限制资源是核数而
不是读长"只在那一个操作点上成立。这个脚本把两者铺成网格，用来区分信息量是否
近似只依赖 N × depth 的乘积，还是存在不可用深度替代的 N 效应。

网格受数据本身限制，而这个限制本身就是答案的一半。要求 ATAC 与 RNA **同时**达到
倍数目标，`Exc_LINC00507_FREM3` 的达标核数是 9,429（1x）、6,226（2x）、3,850（3x）、
2,275（4x）、631（6x）、169（8x）；寡树突胶质细胞在 2x 以上就只剩 283 个。也就是说
在一个固定文库里，核数与每核深度**不能各自独立地买到**：可行前沿近似满足
N × depth ≈ 一万，直到 4 倍以上深核的尾部耗尽。原梯队的深度序列停在 150 核正是
这个缘故。

因此这里不铺满析因网格，而是沿**等乘积线**取点：如果信息量只依赖 N × depth 的
乘积，沿一条等乘积线启动子富集应当不变；如果存在不可用深度替代的 N 效应，沿线
富集应当随 N 上升。这正是审稿意见要区分的两件事，而且比满网格省三分之二的算力。

  等乘积 ≈ 1200：(150, 8x) (300, 4x) (600, 2x) (1200, 1x)
  等乘积 ≈ 2400：(600, 4x) (1200, 2x) (2400, 1x)
  纯核数线（1x）：150 / 400 / 900 / 2400
  纯深度线（150 核）：1x / 2x / 4x / 8x

所有格子都用 --require-depth 抽样（只从达标核里抽），否则降采样会放过未达标的核，
等深度被静默破坏。代价是 1 倍那一行与已发表梯队相差一个核（9,429/9,430）。

用法：
  python3 code/reassessment/c1_depth_grid.py --dir <d> --out depth_grid.csv
"""
from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
import time
from pathlib import Path
from typing import List, Optional, Tuple

import pandas as pd

CELLTYPE = "Exc_LINC00507_FREM3"
# (核数, 深度倍数) —— 只列数据支持的格子
CELLS: List[Tuple[int, float]] = [
    # 纯深度线（150 核）与等乘积 1200 线的端点
    (150, 1.0), (150, 2.0), (150, 4.0), (150, 8.0),
    # 等乘积 ≈ 1200
    (300, 4.0), (600, 2.0), (1200, 1.0),
    # 等乘积 ≈ 2400
    (600, 4.0), (1200, 2.0), (2400, 1.0),
    # 纯核数线（1x）补齐
    (400, 1.0), (900, 1.0),
]
SEEDS = (0, 1, 2)
NULLS = ("none", "trans")
RUNNER = Path(__file__).with_name("primary_ladder_aggregated.py")
logger = logging.getLogger("v3")


def run_one(d: Path, n: int, mult: float, seed: int, null: str,
            tss: str, feats: str) -> Optional[dict]:
    tag = "sc_matched" + ("" if null == "none" else f"_{null}")
    tag += "" if mult == 1.0 else f"_d{mult:g}x"
    out = d / f"v3_{CELLTYPE}_n{n}_s{seed}_{tag}.json"
    if out.exists():
        return json.loads(out.read_text())
    cmd = [sys.executable, str(RUNNER), "--dir", str(d), "--celltype", CELLTYPE,
           "--n", str(n), "--seed", str(seed), "--linker", "single",
           "--null", null, "--depth-mult", str(mult), "--require-depth",
           "--match-tested", "--tss", tss, "--features", feats,
           "--out", str(out)]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        logger.warning("失败 n=%d %.0fx s=%d /%s: %s", n, mult, seed, null,
                       (r.stderr or "").strip().splitlines()[-1:])
        return None
    logger.info("  ok n=%-5d %.0fx s=%d %-5s (%.0fs)", n, mult, seed, null,
                time.time() - t0)
    return json.loads(out.read_text())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tss", default="Resolution_limit/data/derived_results/p5_tss_gencode_v32.csv")
    ap.add_argument("--features", default="Data/Other_Datasets/Multiome_Dataset/files")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S")
    d = Path(a.dir)

    jobs = [(n, m, s, nl) for (n, m) in CELLS for s in SEEDS for nl in NULLS]
    logger.info("共 %d 个配置（%d 个格子 × %d 种子 × %d 零假设）",
                len(jobs), len(CELLS), len(SEEDS), len(NULLS))
    rows = []
    for i, (n, m, s, nl) in enumerate(jobs, 1):
        logger.info("[%d/%d]", i, len(jobs))
        r = run_one(d, n, m, s, nl, a.tss, a.features)
        if r is not None:
            rows.append(r)
    pd.DataFrame(rows).to_csv(a.out, index=False)
    logger.info("写出 %s (%d 行)", a.out, len(rows))


if __name__ == "__main__":
    main()
