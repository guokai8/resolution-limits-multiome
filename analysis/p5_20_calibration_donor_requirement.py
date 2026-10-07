#!/usr/bin/env python3
"""审稿意见 V2：标定表达层系数需要多少个重复供体？

正文说"从一对技术重复测出常数，再用 n^(−1/2) 外推"。一对重复是系数的极嘈杂
估计量，这个脚本把"需要多少供体"量化出来，从而把那句话换成可执行的数字。

做法：对供体做有放回重抽样，抽样规模 m = 1..全部，每次用论文自己的估计量
（`technical_floor.estimate_expression_floor`，同一条 ≥50 核准则、同一 log–log
拟合）重拟合，记录系数与指数的分布。报出使系数 95% 区间落在全数据估计
±10% / ±20% / ±30% 之内所需的最小 m。

用法：
    python3 code/analysis/p5_20_calibration_donor_requirement.py
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
for cand in (ROOT, ROOT / "code"):
    if (cand / "technical_floor").is_dir():
        sys.path.insert(0, str(cand))
        break
else:
    raise SystemExit("[FATAL] 找不到 technical_floor 包")

from technical_floor.config import FloorConfig          # noqa: E402
from technical_floor.io import apply_pair_criterion     # noqa: E402

N_BOOT = 2000
TOLERANCES = (0.10, 0.20, 0.30)
logger = logging.getLogger("v2")


def loglog_fit(n: np.ndarray, y: np.ndarray) -> tuple:
    b, log_a = np.polyfit(np.log(n), np.log(y), 1)
    return float(b), float(np.exp(log_a))


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    cfg = FloorConfig()
    src = ROOT / "Resolution_limit" / "data" / "derived_results" / "p5_floor_scaling_pairs.csv"
    if not src.exists():
        src = ROOT / "data" / "derived_results" / "p5_floor_scaling_pairs.csv"
    df = pd.read_csv(src).rename(columns={"median_abs_log2FC": "floor"})
    df = apply_pair_criterion(df, cfg.min_nuclei_per_pair)
    df = df[df.included_in_fit & (df.n_eff > 0) & (df.floor > 0)]

    n_all = df.n_eff.to_numpy()
    y_all = df.floor.to_numpy()
    donors = df.donor.to_numpy()
    uniq = np.unique(donors)
    idx_by_donor: Dict[str, np.ndarray] = {d: np.flatnonzero(donors == d) for d in uniq}

    b_full, a_full = loglog_fit(n_all, y_all)
    logger.info("全数据：指数 %.4f，系数 %.4f（%d 观测 / %d 供体）\n",
                b_full, a_full, len(df), len(uniq))

    rng = np.random.default_rng(cfg.seed)
    rows: List[Dict[str, float]] = []
    for m in range(1, len(uniq) + 1):
        a_draws, b_draws = [], []
        for _ in range(N_BOOT):
            pick = rng.choice(uniq, size=m, replace=True)
            idx = np.concatenate([idx_by_donor[d] for d in pick])
            if idx.size < 5 or np.unique(n_all[idx]).size < 2:
                continue
            b, a = loglog_fit(n_all[idx], y_all[idx])
            a_draws.append(a)
            b_draws.append(b)
        if len(a_draws) < 50:
            continue
        a_arr = np.asarray(a_draws)
        b_arr = np.asarray(b_draws)
        a_lo, a_hi = np.percentile(a_arr, [2.5, 97.5])
        rows.append(dict(
            n_donors=m, n_fits=len(a_draws),
            coef_median=float(np.median(a_arr)),
            coef_lo=float(a_lo), coef_hi=float(a_hi),
            coef_rel_halfwidth=float((a_hi - a_lo) / 2 / a_full),
            exponent_median=float(np.median(b_arr)),
            exponent_lo=float(np.percentile(b_arr, 2.5)),
            exponent_hi=float(np.percentile(b_arr, 97.5)),
        ))

    out = pd.DataFrame(rows)
    dest = src.parent / "p5_20_calibration_donor_requirement.csv"
    out.to_csv(dest, index=False)

    print(f"{'供体数':>6s} {'系数中位':>9s} {'系数 95% 区间':>22s} "
          f"{'相对半宽':>9s} {'指数 95% 区间':>22s}")
    for _, r in out.iterrows():
        print(f"{int(r.n_donors):>6d} {r.coef_median:>9.3f} "
              f"{r.coef_lo:>10.3f}–{r.coef_hi:<10.3f} {r.coef_rel_halfwidth:>9.1%} "
              f"{r.exponent_lo:>10.3f}–{r.exponent_hi:<10.3f}")

    print()
    for tol in TOLERANCES:
        ok = out[out.coef_rel_halfwidth <= tol]
        if len(ok):
            print(f"  系数落在 ±{tol:.0%} 之内所需最小供体数：{int(ok.n_donors.min())}")
        else:
            print(f"  系数落在 ±{tol:.0%} 之内：{len(uniq)} 个供体仍不够")
    print(f"\n写出 {dest}")


if __name__ == "__main__":
    main()
