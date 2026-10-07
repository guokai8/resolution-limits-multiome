#!/usr/bin/env python3
"""C1 驱动脚本：把 single-cell 与 KNN 聚合两种 linker 在同一梯队上跑满，收成一张表。

每个配置都跑两次：真实配对与 trans 零假设（基因配到异染色体基因的窗口）。
诊断量是 OR_ratio = OR(观测) / OR(trans)：它衡量链接集合里**位置信息**的含量。
已存在的 JSON 会跳过，所以脚本可以中断后继续。

用法：
  python3 c1_run_grid.py --dir results_methods/primary_ladder --out c1_grid.csv
"""
from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

LADDER: Dict[str, List[int]] = {
    "Exc_LINC00507_FREM3": [150, 400, 900, 2400, 5000, 7500],
    "Oligodendrocytes": [150, 400, 900, 2400],
}
SEEDS = (0, 1, 2)
K_MAIN = 25
K_SWEEP = (10, 50, 100)
NULLS = ("none", "trans")
RUNNER = Path(__file__).with_name("primary_ladder_aggregated.py")

logger = logging.getLogger("c1")


def tag_for(linker: str, k: int, matched: bool, null: str) -> str:
    t = f"k{k}" if linker == "aggregate" else "sc"
    return t + ("_matched" if matched else "") + ("" if null == "none" else f"_{null}")


def run_one(d: Path, ct: str, n: int, seed: int, linker: str, k: int,
            matched: bool, null: str, tss: str, features: str) -> Optional[dict]:
    out = d / f"c1_{ct}_n{n}_s{seed}_{tag_for(linker, k, matched, null)}.json"
    if out.exists():
        return json.loads(out.read_text())
    cmd = [sys.executable, str(RUNNER), "--dir", str(d), "--celltype", ct,
           "--n", str(n), "--seed", str(seed), "--linker", linker,
           "--null", null, "--tss", tss, "--features", features]
    if linker == "aggregate":
        cmd += ["--k", str(k)]
    if matched:
        cmd += ["--match-tested"]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        logger.warning("失败 %s n=%d s=%d %s/%s: %s", ct, n, seed, linker, null,
                       (r.stderr or "").strip().splitlines()[-1:] or "?")
        return None
    logger.info("  ok %-22s n=%-5d s=%d %-10s k=%-3s %-5s (%.0fs)",
                ct, n, seed, linker, k if linker == "aggregate" else "-", null,
                time.time() - t0)
    return json.loads(out.read_text())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="results_methods/primary_ladder")
    ap.add_argument("--out", default="results_methods/primary_ladder/c1_grid.csv")
    ap.add_argument("--tss", default="Resolution_limit/data/derived_results/p5_tss_gencode_v32.csv")
    ap.add_argument("--features", default="Data/Other_Datasets/Multiome_Dataset/files")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S")
    d = Path(a.dir)

    jobs = []
    # 主网格：两臂 × 共享被检验集合 × 两个零假设
    for ct, ns in LADDER.items():
        for n in ns:
            for seed in SEEDS:
                for null in NULLS:
                    jobs.append((ct, n, seed, "single", 0, True, null))
                    if K_MAIN <= n // 3:
                        jobs.append((ct, n, seed, "aggregate", K_MAIN, True, null))
    # 原生检出规则（真跑 ArchR 会得到的集合）
    for n in (150, 2400):
        for seed in SEEDS:
            for null in NULLS:
                jobs.append(("Exc_LINC00507_FREM3", n, seed, "aggregate", K_MAIN, False, null))
    # k 敏感性
    for k in K_SWEEP:
        if k <= 2400 // 3:
            for null in NULLS:
                jobs.append(("Exc_LINC00507_FREM3", 2400, 0, "aggregate", k, True, null))

    logger.info("共 %d 个配置", len(jobs))
    rows = []
    for i, (ct, n, seed, linker, k, matched, null) in enumerate(jobs, 1):
        logger.info("[%d/%d]", i, len(jobs))
        r = run_one(d, ct, n, seed, linker, k, matched, null, a.tss, a.features)
        if r is not None:
            rows.append(r)

    df = pd.DataFrame(rows)
    df.to_csv(a.out, index=False)
    logger.info("写出 %s (%d 行)", a.out, len(df))


if __name__ == "__main__":
    main()
