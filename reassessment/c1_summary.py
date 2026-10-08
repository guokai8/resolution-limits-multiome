#!/usr/bin/env python3
"""C1 summary: score each arm's link set against the trans null.

The 2x2 behind the positional odds ratio, stratified by seed and pooled with
Mantel-Haenszel and the Robins-Breslow-Greenland variance, as in the main text:

                        promoter-proximal   distal
  observed link set             a             b
  trans link set                c             d

The trans pairing preserves the whole distribution of window sizes and of
proximal-to-distal ratios; the only thing it breaks is the pairing of a gene
with its own window. So this odds ratio reads directly as: how much of the
promoter enrichment in the detected links is POSITIONAL?

An odds ratio near 1 means that link set carries no positional information,
however many links it contains.

Usage:
  python3 c1_summary.py --grid results_methods/primary_ladder/c1_grid.csv
"""
from __future__ import annotations

import argparse
from ast import literal_eval
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy import stats as st


def g_test(obs: np.ndarray, ref: np.ndarray) -> Tuple[float, float]:
    """G test of the detected links' distance-band histogram against a reference.

    Uses all the links rather than only the promoter-proximal few, so it has far
    more power than the binary proximal/distal split -- a link set can shift its
    whole distance profile without crossing the 3 kb promoter boundary.

    Returns (G, p). NaN when either histogram is empty.
    """
    obs = np.asarray(obs, dtype=float)
    ref = np.asarray(ref, dtype=float)
    if obs.sum() < 1 or ref.sum() < 1:
        return float("nan"), float("nan")
    exp = ref / ref.sum() * obs.sum()
    m = exp > 0
    g = 2 * np.nansum(obs[m] * np.log(np.where(obs[m] > 0, obs[m] / exp[m], 1.0)))
    return float(g), float(st.chi2.sf(g, m.sum() - 1))


def band_median(bands: List[int], edges: List[int]) -> float:
    """Median distance estimated from the band histogram.

    Linear interpolation within the band containing the median; the raw
    distances are not retained, only the binned counts.
    """
    c = np.asarray(bands, dtype=float)
    if c.sum() == 0:
        return float("nan")
    cum = np.cumsum(c)
    half = c.sum() / 2
    i = int(np.searchsorted(cum, half))
    i = min(i, len(c) - 1)
    lo, hi = edges[i], edges[i + 1]
    before = cum[i - 1] if i > 0 else 0.0
    frac = (half - before) / max(c[i], 1e-9)
    return float(lo + frac * (hi - lo))


def mh(strata: List[Tuple[int, int, int, int]]) -> Tuple[float, float, float]:
    """Mantel–Haenszel 合并 OR，Robins–Breslow–Greenland 方差。照正文口径。"""
    Rs = Ss = 0.0
    PR = PS = QR = QS = 0.0
    for a, b, c, dd in strata:
        n = a + b + c + dd
        if n == 0:
            continue
        R = a * dd / n
        S = b * c / n
        Rs += R
        Ss += S
        P = (a + dd) / n
        Q = (b + c) / n
        PR += P * R
        PS += P * S
        QR += Q * R
        QS += Q * S
    if Ss == 0 or Rs == 0:
        return float("nan"), float("nan"), float("nan")
    OR = Rs / Ss
    se = np.sqrt(PR / (2 * Rs ** 2) + (PS + QR) / (2 * Rs * Ss) + QS / (2 * Ss ** 2))
    return OR, OR * np.exp(-1.96 * se), OR * np.exp(1.96 * se)


KEY = ["celltype", "n", "linker", "k", "match_tested", "archr_defaults", "window"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", required=True)
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    df = pd.read_csv(a.grid)
    obs = df[df.null == "none"].set_index(KEY + ["seed"])
    trans = df[df.null == "trans"].set_index(KEY + ["seed"])

    rows = []
    for key, g in obs.groupby(level=list(range(len(KEY)))):
        strata: List[Tuple[int, int, int, int]] = []
        links_o: List[int] = []
        links_t: List[int] = []
        ors: List[float] = []
        for idx, r in g.iterrows():
            if idx not in trans.index:
                continue
            t = trans.loc[idx]
            ao = int(r.prox_links)
            bo = int(r.n_links) - ao
            at = int(t.prox_links)
            bt = int(t.n_links) - at
            strata.append((ao, bo, at, bt))
            links_o.append(int(r.n_links))
            links_t.append(int(t.n_links))
            ors.append(float(r.OR))
        if not strata:
            continue
        # With an empty cell the Mantel-Haenszel RBG variance degenerates to
        # NaN. Apply the Haldane-Anscombe 0.5 to EVERY stratum, not just the
        # offending one, so the strata stay comparable, and record that it was
        # applied.
        ha = any(min(s) == 0 for s in strata)
        use = [tuple(v + 0.5 for v in s) for s in strata] if ha else strata
        pos, lo, hi = mh(use)
        ct, n, linker, k, matched, archr, window = key

        # Distance bands: sum the three seeds' histograms, then compare
        def stack(frame: pd.DataFrame, col: str) -> np.ndarray:
            vals = []
            for v in frame[col]:
                if isinstance(v, str):
                    vals.append(np.asarray(literal_eval(v), dtype=float))
                elif isinstance(v, (list, tuple, np.ndarray)):
                    vals.append(np.asarray(v, dtype=float))
            return np.sum(vals, axis=0) if vals else np.zeros(6)

        idxs = [i for i in g.index if i in trans.index]
        bl_o = stack(g.loc[idxs], "band_links")
        bt = stack(g.loc[idxs], "band_tested")
        bl_t = stack(trans.loc[idxs], "band_links")
        edges = [0, 3_000, 10_000, 50_000, 100_000, 250_000, 500_001]
        g_bg, p_bg = g_test(bl_o, bt)
        g_tr, p_tr = g_test(bl_o, bl_t)
        rows.append(dict(
            celltype=ct, n=int(n), linker=linker, k=int(k), match_tested=bool(matched),
            archr_defaults=bool(archr), window=int(window), seeds=len(strata),
            haldane_anscombe=bool(ha),
            links_obs=int(np.mean(links_o)), links_trans=int(np.mean(links_t)),
            inflation=float(np.mean(links_o)) / max(np.mean(links_t), 1e-9),
            prox_pct_obs=100 * sum(s[0] for s in strata) / max(sum(s[0] + s[1] for s in strata), 1),
            prox_pct_trans=100 * sum(s[2] for s in strata) / max(sum(s[2] + s[3] for s in strata), 1),
            promoter_OR=float(np.mean(ors)),
            positional_OR=pos, positional_lo=lo, positional_hi=hi,
            med_obs=band_median(bl_o, edges), med_trans=band_median(bl_t, edges),
            med_bg=band_median(bt, edges),
            G_vs_background=g_bg, P_vs_background=p_bg,
            G_vs_trans=g_tr, P_vs_trans=p_tr))

    out = pd.DataFrame(rows).sort_values(
        ["celltype", "window", "archr_defaults", "linker", "k", "n"])
    if a.out:
        out.to_csv(a.out, index=False)

    for (ct, window), g in out.groupby(["celltype", "window"]):
        print(f"\n== {ct}  ±{window // 1000} kb 窗口 ==")
        print(f"{'linker':14s} {'阈值':>6s} {'k':>4s} {'n':>6s} {'links':>9s} {'trans':>8s} "
              f"{'prox%':>6s} {'promOR':>7s} {'位置 OR (95% CI)':>23s} "
              f"{'中位kb':>7s} {'transkb':>8s}")
        for _, r in g.iterrows():
            ptr = "—" if not np.isfinite(r.P_vs_trans) else (
                "<1e-300" if r.P_vs_trans < 1e-300 else f"{r.P_vs_trans:.1e}")
            thr = "ArchR" if r.archr_defaults else "正文"
            print(f"{r.linker:14s} {thr:>6s} {r.k:>4d} {r.n:>6d} {r.links_obs:>9,} "
                  f"{r.links_trans:>8,} {r.prox_pct_obs:>6.2f} {r.promoter_OR:>7.3f} "
                  f"{r.positional_OR:>7.3f} ({r.positional_lo:.2f}–{r.positional_hi:.2f})"
                  f"{'*' if r.haldane_anscombe else ' '} "
                  f"{r.med_obs/1000:>7.0f} {r.med_trans/1000:>8.0f}")


if __name__ == "__main__":
    main()
