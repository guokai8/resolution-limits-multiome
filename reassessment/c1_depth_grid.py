#!/usr/bin/env python3
"""The nucleus-number by sequencing-depth response surface.

The main text runs two lines only -- vary nuclei at fixed depth, and vary depth
at 150 nuclei -- so the claim that nuclei rather than reads are the limiting
resource holds at one operating point. This script spreads the two over a grid
to separate two possibilities: that information depends only on the product
N x depth, or that there is an N effect that depth cannot substitute for.

The grid is limited by the data, and that limit is half the answer. Requiring
BOTH ATAC and RNA to reach the depth multiple, Exc_LINC00507_FREM3 has 9,429
qualifying nuclei at 1x, 6,226 at 2x, 3,850 at 3x, 2,275 at 4x, 631 at 6x and
169 at 8x; oligodendrocytes fall to 283 beyond 2x. In other words, within one
fixed library nucleus count and per-nucleus depth CANNOT be bought
independently: the feasible frontier sits near N x depth = ten thousand until
the tail of deep nuclei runs out past 4x. That is also why the original depth
series stopped at 150 nuclei.

So rather than a full factorial, points are taken along LINES OF CONSTANT
PRODUCT. If information depends only on the product, promoter enrichment should
be flat along such a line; if there is an N effect depth cannot replace, it
should rise with N along the line. That is exactly the distinction at issue,
and it costs a third of a full grid.

  product ~ 1200:  (150, 8x) (300, 4x) (600, 2x) (1200, 1x)
  product ~ 2400:  (600, 4x) (1200, 2x) (2400, 1x)
  pure nucleus line (1x):     150 / 400 / 900 / 2400
  pure depth line (150 nuclei): 1x / 2x / 4x / 8x

Every cell samples with --require-depth, drawing only from nuclei that already
clear the target. Without it, downsampling would quietly admit nuclei below
target and the equal-depth guarantee would be broken. The cost is that the 1x
row differs from the published ladder by a single nucleus, 9,429 against 9,430.

Usage:
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
# (nuclei, depth multiple) -- only the cells the data can actually support
CELLS: List[Tuple[int, float]] = [
    # Pure depth line (150 nuclei), which also supplies the endpoint of the
    # constant-product-1200 line
    (150, 1.0), (150, 2.0), (150, 4.0), (150, 8.0),
    # constant product ~ 1200
    (300, 4.0), (600, 2.0), (1200, 1.0),
    # constant product ~ 2400
    (600, 4.0), (1200, 2.0), (2400, 1.0),
    # pure nucleus line at 1x, filling in the rungs not already covered
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
