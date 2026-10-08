#!/usr/bin/env python3
"""Does information depend only on N x depth, or is there an N effect?

For each line of constant product, regress log(promoter OR) on log(N):

  slope ~ 0   information depends only on the product; depth substitutes for
              nuclei, and a budget can be spent either way
  slope > 0   at the same product, more nuclei is better; nuclei cannot be
              replaced by depth

The pure nucleus line (vary N at 1x) and the pure depth line (vary depth at 150
nuclei) are reported alongside as a reference. Those are the two lines the main
text originally had, and on their own they cannot tell the two cases apart --
both rise, whichever is true.

Intervals come from bootstrapping the seeds. The seeds are repeated
computational draws of the same nuclei, so they are treated as repeated
measures rather than as independent strata.

Usage:
  python3 c1_depth_summary.py --grid v3_depth_grid.csv
"""
from __future__ import annotations

import argparse
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

NBOOT = 4000
SEED = 0


def mh(strata: List[Tuple[float, float, float, float]]) -> float:
    """Mantel-Haenszel pooled odds ratio over strata of (a, b, c, d).

    Point estimate only; the interval comes from the seed bootstrap instead,
    since the seeds are not independent strata.
    """
    R = S = 0.0
    for a, b, c, d in strata:
        n = a + b + c + d
        if n == 0:
            continue
        R += a * d / n
        S += b * c / n
    return R / S if S > 0 else float("nan")


def or_for(g: pd.DataFrame) -> float:
    """Pooled promoter odds ratio for one grid cell over its three seeds."""
    strata = []
    for _, r in g.iterrows():
        a = float(r.prox_links)
        b = float(r.n_links) - a
        c = float(r.prox_tested_all)
        d = float(r.n_tested_all) - c
        if min(a, b, c, d) == 0:
            a, b, c, d = a + .5, b + .5, c + .5, d + .5
        strata.append((a, b, c, d))
    return mh(strata)


def slope_ci(cells: Dict[int, pd.DataFrame]) -> Tuple[float, float, float, int]:
    """log(OR) ~ log(N) 的斜率，按种子重抽样给区间。"""
    ns = sorted(cells)
    if len(ns) < 2:
        return (float("nan"),) * 3 + (len(ns),)
    y = np.array([np.log(or_for(cells[n])) for n in ns])
    x = np.log(np.array(ns, dtype=float))
    ok = np.isfinite(y)
    if ok.sum() < 2:
        return (float("nan"),) * 3 + (int(ok.sum()),)
    point = float(np.polyfit(x[ok], y[ok], 1)[0])

    rng = np.random.default_rng(SEED)
    draws: List[float] = []
    for _ in range(NBOOT):
        yb = []
        for n in ns:
            g = cells[n]
            seeds = g.seed.to_numpy()
            pick = rng.choice(seeds, size=len(seeds), replace=True)
            gb = pd.concat([g[g.seed == s] for s in pick])
            yb.append(np.log(or_for(gb)))
        yb = np.asarray(yb)
        m = np.isfinite(yb)
        if m.sum() >= 2:
            draws.append(float(np.polyfit(x[m], yb[m], 1)[0]))
    if len(draws) < 50:
        return point, float("nan"), float("nan"), int(ok.sum())
    return (point, float(np.percentile(draws, 2.5)),
            float(np.percentile(draws, 97.5)), int(ok.sum()))


LINES = {
    "等乘积 ≈ 1200": [(150, 8.0), (300, 4.0), (600, 2.0), (1200, 1.0)],
    "等乘积 ≈ 2400": [(600, 4.0), (1200, 2.0), (2400, 1.0)],
    "纯核数线 (1x)": [(150, 1.0), (400, 1.0), (900, 1.0), (1200, 1.0), (2400, 1.0)],
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", required=True)
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    df = pd.read_csv(a.grid)
    obs = df[df.null == "none"]

    print("每个格子的启动子 OR（三个种子 MH 合并，相对基因组背景）")
    print(f"{'N':>6s} {'depth':>6s} {'N×depth':>8s} {'links':>7s} {'prox':>6s} {'OR':>8s}")
    rows = []
    for (n, m), g in obs.groupby(["n", "depth_mult"]):
        o = or_for(g)
        rows.append(dict(n=int(n), depth_mult=float(m), product=int(n * m),
                         links=int(g.n_links.mean()), prox=int(g.prox_links.mean()),
                         promoter_OR=o, seeds=len(g)))
        print(f"{int(n):>6d} {m:>5.0f}x {int(n*m):>8d} {int(g.n_links.mean()):>7,} "
              f"{int(g.prox_links.mean()):>6d} {o:>8.3f}")

    print("\nlog(OR) ~ log(N) 的斜率（种子 bootstrap 95% 区间）")
    print(f"{'线':<16s} {'点数':>4s} {'斜率':>8s} {'95% 区间':>20s}")
    out_lines = []
    for label, cells_spec in LINES.items():
        cells = {}
        for n, m in cells_spec:
            g = obs[(obs.n == n) & (obs.depth_mult == m)]
            if len(g):
                cells[n] = g
        s, lo, hi, k = slope_ci(cells)
        print(f"{label:<16s} {k:>4d} {s:>8.3f} {lo:>9.3f}–{hi:<9.3f}")
        out_lines.append(dict(line=label, n_points=k, slope=s, lo=lo, hi=hi))

    # Pure depth line: log(OR) on log(depth) at 150 nuclei
    cells = {}
    for m in (1.0, 2.0, 4.0, 8.0):
        g = obs[(obs.n == 150) & (obs.depth_mult == m)]
        if len(g):
            cells[int(m)] = g
    s, lo, hi, k = slope_ci(cells)
    print(f"{'纯深度线 (150核)':<16s} {k:>4d} {s:>8.3f} {lo:>9.3f}–{hi:<9.3f}"
          f"   ← 横轴是深度倍数，不是 N")
    out_lines.append(dict(line="纯深度线 (150核)", n_points=k, slope=s, lo=lo, hi=hi))

    print("\n读法：等乘积线斜率 ≈ 0 表示深度可替代核数；> 0 表示同样乘积下核数更值钱。")
    if a.out:
        pd.DataFrame(rows).to_csv(a.out, index=False)
        pd.DataFrame(out_lines).to_csv(a.out.replace('.csv', '_slopes.csv'), index=False)
        print(f"写出 {a.out}")


if __name__ == "__main__":
    main()
