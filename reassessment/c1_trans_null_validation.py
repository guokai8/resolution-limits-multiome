#!/usr/bin/env python3
"""Validate the trans-pairing null.

The cross-procedure Layer 3 argument now rests on this null, so the null itself
has to be tested rather than assumed. Three questions:

1. Does it depend on one particular pairing? Ten independent trans pairings are
   drawn on the same nuclei at the same depth, and the spread of the positional
   odds ratio across them is measured. If one draw happens to pair a gene into
   a promoter-rich window, the conclusion must not move with it.

2. Is it driven by one chromosome? A leave-one-chromosome-out jackknife, using
   the per-chromosome decomposition.

3. Does it behave as expected where there is plenty of information? At the
   highest nucleus number the observed set should carry clear positional
   information and the trans set should not. A null that returns a large
   positional odds ratio there is a broken null.

rep 0 reuses the main grid's existing output by filename, so only reps 1-9 are
actually run here.

Usage:
  python3 c1_trans_null_validation.py --dir <d> --out trans_validation.csv
"""
from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
import time
from pathlib import Path
from typing import List, Optional

import pandas as pd

CELLTYPE = "Exc_LINC00507_FREM3"
LEVELS = (150, 900, 7500)
ARMS = (("single", 0, False), ("aggregate", 25, True))
N_REPS = 10
SEED = 0
RUNNER = Path(__file__).with_name("primary_ladder_aggregated.py")
logger = logging.getLogger("n1")


def path_for(d: Path, n: int, linker: str, k: int, archr: bool,
             null: str, rep: int) -> Path:
    tag = f"k{k}" if linker == "aggregate" else "sc"
    tag += "_archr" if archr else ""
    tag += "_matched"
    tag += "" if null == "none" else f"_{null}"
    tag += "" if rep == 0 else f"_rep{rep}"
    # rep 0 deliberately produces the same filename as the main grid, so that
    # run is reused rather than repeated.
    return d / f"c1_{CELLTYPE}_n{n}_s{SEED}_{tag}.json"


def run_one(d: Path, n: int, linker: str, k: int, archr: bool,
            null: str, rep: int, tss: str, feats: str) -> Optional[dict]:
    out = path_for(d, n, linker, k, archr, null, rep)
    if out.exists():
        return json.loads(out.read_text())
    cmd = [sys.executable, str(RUNNER), "--dir", str(d), "--celltype", CELLTYPE,
           "--n", str(n), "--seed", str(SEED), "--linker", linker,
           "--null", null, "--match-tested", "--trans-rep", str(rep),
           "--tss", tss, "--features", feats]
    if linker == "aggregate":
        cmd += ["--k", str(k)]
    if archr:
        cmd += ["--archr-defaults"]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        logger.warning("失败 n=%d %s rep=%d: %s", n, linker, rep,
                       (r.stderr or "").strip().splitlines()[-1:])
        return None
    logger.info("  ok n=%-5d %-10s archr=%-5s %-5s rep=%-2d (%.0fs)",
                n, linker, archr, null, rep, time.time() - t0)
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

    jobs = []
    for n in LEVELS:
        for linker, k, archr in ARMS:
            jobs.append((n, linker, k, archr, "none", 0))
            for rep in range(N_REPS):
                jobs.append((n, linker, k, archr, "trans", rep))
    logger.info("共 %d 个配置（已有的会跳过）", len(jobs))

    rows: List[dict] = []
    for i, j in enumerate(jobs, 1):
        r = run_one(d, *j, a.tss, a.features)
        if r is not None:
            rows.append(r)
    pd.DataFrame(rows).to_csv(a.out, index=False)
    logger.info("写出 %s (%d 行)", a.out, len(rows))


if __name__ == "__main__":
    main()
