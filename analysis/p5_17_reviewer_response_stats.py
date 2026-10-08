#!/usr/bin/env python3
"""
Recomputations triggered by review, from deposited tables only
=============================================================
Nothing here re-runs an upstream scan; every block reads a deposited derived
table and writes back into data/derived_results/.

B2  sensitivity of the expression-layer exponent to the minimum-effective-n
    support floor, plus a paired comparison of the Seattle atlas observations
    against their own matched null ON THE SAME SUPPORT.
B3  the Layer 3 ladder. A seed is a re-draw of the same nuclei, not an
    independent stratum, so this pools cell types WITHIN a seed by
    Mantel-Haenszel and then combines across seeds.
B5  confidence interval for the composition-layer kappa, and its per-cell-type
    values.
B6  donor-clustered bootstrap test of the difference in composition exponent
    between the two cohorts.

Usage:
    python3 code/analysis/p5_17_reviewer_response_stats.py
"""

from __future__ import annotations

import csv
import math
import os
import re
from typing import Dict, List, Sequence, Tuple

import numpy as np

DR = "data/derived_results"
ST = "supplementary_tables"
NBOOT = 4000
SEED = 0


# ----------------------------------------------------------------------------- helpers
def read_csv(path: str) -> List[Dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def loglog_fit(n: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    """Fit log y = b log n + log a, returning (exponent b, coefficient a)."""
    b, la = np.polyfit(np.log(n), np.log(y), 1)
    return float(b), float(math.exp(la))


def boot_by_group(
    n: np.ndarray, y: np.ndarray, groups: np.ndarray, stat, nboot: int = NBOOT
) -> Tuple[float, float]:
    """Bootstrap by group (donor), returning the 95% interval of stat."""
    rng = np.random.default_rng(SEED)
    uniq = np.unique(groups)
    out = []
    for _ in range(nboot):
        pick = rng.choice(uniq, len(uniq), replace=True)
        idx = np.concatenate([np.flatnonzero(groups == g) for g in pick])
        if idx.size < 5:
            continue
        out.append(stat(n[idx], y[idx]))
    lo, hi = np.percentile(out, [2.5, 97.5])
    return float(lo), float(hi)


# ----------------------------------------------------------------------------- B2
def b2_exponent_support_sensitivity() -> List[Dict[str, object]]:
    """Exponent against the support floor, and observed against null.

    The motor cortex exponent is recomputed as the minimum n_eff is raised, to
    show the fit is not carried by the low-count tail. The Seattle atlas
    observations are then compared with their own matched null, paired row by
    row on the same support -- comparing two separately fitted intercepts
    instead would confound the support difference with the effect."""
    rows: List[Dict[str, object]] = []

    # The 50-nucleus pairing criterion governs the expression layer too: keep
    # only included_in_fit == yes, which is 26 donors
    mc = [x for x in read_csv(f"{ST}/Supplementary_Table_3a_expression_floors_motor_cortex.csv")
          if x.get("included_in_fit", "yes").startswith("yes")]
    n = np.array([float(x["n_eff"]) for x in mc])
    f = np.array([float(x["median_abs_log2FC"]) for x in mc])
    d = np.array([x["donor"] for x in mc])

    for thr in (0, 10, 25, 50):
        m = n >= thr if thr else n > 0
        b, a = loglog_fit(n[m], f[m])
        lo, hi = boot_by_group(n[m], f[m], d[m], lambda nn, yy: loglog_fit(nn, yy)[0])
        rows.append(dict(cohort="motor_cortex", series="observed", min_n_eff=thr,
                         n_obs=int(m.sum()), n_donors=int(len(set(d[m]))),
                         exponent=round(b, 4), exponent_lo=round(lo, 4),
                         exponent_hi=round(hi, 4), coefficient=round(a, 3)))

    # Seattle atlas: observation and null are paired row by row, and must be
    # compared on the same support
    ob = read_csv(f"{DR}/seaad_L2_obs.csv")
    nu = read_csv(f"{DR}/seaad_L2_null.csv")
    key = lambda x: (x["donor"], x["ct"])
    assert [key(x) for x in ob] == [key(x) for x in nu], "obs/null 行序不一致"
    ns = np.array([float(x["n_eff"]) for x in ob])
    fo = np.array([float(x["floor"]) for x in ob])
    fn = np.array([float(x["floor"]) for x in nu])
    ds = np.array([x["donor"] for x in ob])

    rng = np.random.default_rng(SEED)
    for thr in (0, 10, 25, 50):
        m = ns >= thr if thr else ns > 0
        bo, ao = loglog_fit(ns[m], fo[m])
        bn, an = loglog_fit(ns[m], fn[m])
        uniq, diffs = np.unique(ds[m]), []
        x, yo, yn, dl = np.log(ns[m]), np.log(fo[m]), np.log(fn[m]), ds[m]
        for _ in range(NBOOT):
            pick = rng.choice(uniq, len(uniq), replace=True)
            idx = np.concatenate([np.flatnonzero(dl == g) for g in pick])
            diffs.append(np.polyfit(x[idx], yo[idx], 1)[0]
                         - np.polyfit(x[idx], yn[idx], 1)[0])
        dlo, dhi = np.percentile(diffs, [2.5, 97.5])
        for series, bb, aa in (("observed", bo, ao), ("matched_null", bn, an)):
            rows.append(dict(cohort="seaad", series=series, min_n_eff=thr,
                             n_obs=int(m.sum()), n_donors=int(len(uniq)),
                             exponent=round(bb, 4), exponent_lo="", exponent_hi="",
                             coefficient=round(aa, 3)))
        rows.append(dict(cohort="seaad", series="observed_minus_null", min_n_eff=thr,
                         n_obs=int(m.sum()), n_donors=int(len(uniq)),
                         exponent=round(bo - bn, 4), exponent_lo=round(dlo, 4),
                         exponent_hi=round(dhi, 4), coefficient=""))
    return rows


# ----------------------------------------------------------------------------- B3
def _mh(tabs: Sequence[Tuple[float, float, float, float]]) -> Tuple[float, float]:
    """Mantel-Haenszel pooled odds ratio with the Robins-Breslow-Greenland
    standard error of the log odds ratio."""
    R = S = 0.0
    PR = QS = PSQR = 0.0
    for a, b, c, d in tabs:
        n = a + b + c + d
        r, s = a * d / n, b * c / n
        p, q = (a + d) / n, (b + c) / n
        R += r
        S += s
        PR += p * r
        QS += q * s
        PSQR += p * s + q * r
    var = PR / (2 * R**2) + PSQR / (2 * R * S) + QS / (2 * S**2)
    return R / S, math.sqrt(var)


def _tabs_from(rows: Sequence[Dict[str, str]],
               ha: bool = False,
               exclude_detected: bool = False
               ) -> List[Tuple[float, float, float, float]]:
    """Build the 2x2 tables from link and tested counts.

    One convention is used throughout, matching the body of Table 3:

      * ha=False. No cell in the ladder or the depth series is empty, so the
        Haldane-Anscombe correction is unnecessary. Applying it anyway would
        pull the estimates towards 1 and make them incomparable with every
        other odds ratio in the paper.
      * exclude_detected=False. The background is not reduced by the links
        already detected.

    Both switches are kept so the alternative convention can be reproduced
    explicitly and the size of the difference quantified, rather than asserted.
    """
    out = []
    for x in rows:
        a = int(x["prox_links"])
        b = int(x["n_links"]) - a
        c = int(x["prox_tested_all"])
        d = int(x["n_tested_all"]) - int(x["prox_tested_all"])
        if exclude_detected:
            c, d = c - a, d - b
        vals = (a, b, c, d)
        out.append(tuple(v + 0.5 for v in vals) if ha else tuple(float(v) for v in vals))
    return out


def b3_ladder_seed_handling() -> List[Dict[str, object]]:
    """A seed is a re-draw of the same nuclei: pool cell types within a seed,
    then combine across seeds."""
    rows: List[Dict[str, object]] = []
    for cohort, path in (("primary", f"{DR}/primary_L3_nucleus_ladder.csv"),
                         ("external", f"{DR}/nabec_L3_nucleus_ladder.csv")):
        if not os.path.exists(path):
            print(f"  [skip] {path} 不存在")
            continue
        data = read_csv(path)
        for nlev in sorted({int(x["n"]) for x in data}):
            at = [x for x in data if int(x["n"]) == nlev]

            # Old convention: every (cell type x seed) treated as an
            # independent stratum
            or_all, se_all = _mh(_tabs_from(at))
            # The variant with Haldane-Anscombe and the detected links removed
            # from the background, kept only to quantify what the convention
            # itself contributes
            or_ha, _ = _mh(_tabs_from(at, ha=True, exclude_detected=True))

            # New convention: pool cell types within each seed
            per: List[Tuple[float, float]] = []
            for sd in sorted({x["seed"] for x in at}):
                sub = [x for x in at if x["seed"] == sd]
                o, se = _mh(_tabs_from(sub))
                per.append((math.log(o), se))
            lv = np.array([p[0] for p in per])
            within = float(np.mean([p[1] ** 2 for p in per]))
            between = float(np.var(lv, ddof=1)) if len(lv) > 1 else 0.0
            # Conservative combination: within-seed sampling variance plus
            # between-seed downsampling variance
            se_comb = math.sqrt(within + between)
            mu = float(lv.mean())
            rows.append(dict(
                cohort=cohort, n_per_celltype=nlev,
                n_celltypes=len({x["celltype"] for x in at}),
                n_seeds=len(per),
                prox_links_total=sum(int(x["prox_links"]) for x in at),
                OR_seeds_as_strata=round(or_all, 4),
                OR_seeds_as_strata_with_HA=round(or_ha, 4),
                lo_seeds_as_strata=round(or_all * math.exp(-1.96 * se_all), 4),
                hi_seeds_as_strata=round(or_all * math.exp(1.96 * se_all), 4),
                OR_seeds_as_repeats=round(math.exp(mu), 4),
                lo_seeds_as_repeats=round(math.exp(mu - 1.96 * se_comb), 4),
                hi_seeds_as_repeats=round(math.exp(mu + 1.96 * se_comb), 4),
                per_seed_OR="; ".join(f"{math.exp(v):.3f}" for v in lv),
            ))
    return rows


# ----------------------------------------------------------------------------- B5
def b5_kappa() -> List[Dict[str, object]]:
    """Per-cell-type kappa, the RMS summary, and its bootstrap interval."""
    # Same source as b9: recomputed from the deposited per-pair tables rather
    # than read from a separate upstream summary. The two disagreed for 4 of 12
    # cell types; this is the one the text quotes.
    zz, nn, _ = _zn(COHORT_PAIR_TABLES[0][1], None)
    ct = np.array([x["celltype"] for x in read_csv(COHORT_PAIR_TABLES[0][1])])
    ratio = zz / np.sqrt(2.0 / nn)
    src = [{"celltype": t,
            "z_rms": math.sqrt(float((ratio[ct == t] ** 2).mean())),
            "mean_frac": float(np.mean([float(x["p"]) for x in read_csv(COHORT_PAIR_TABLES[0][1])
                                        if x["celltype"] == t]))}
           for t in sorted(set(ct))]
    # The published kappa = 4.28 is the RMS of z_rms, NOT the RMS of the
    # overdispersion_ratio column -- a distinction that is easy to lose
    k = np.array([float(x["z_rms"]) for x in src])
    ct = [x["celltype"] for x in src]
    frac = np.array([float(x["mean_frac"]) for x in src])

    rng = np.random.default_rng(SEED)
    boots = [math.sqrt(np.mean(rng.choice(k, k.size, replace=True) ** 2))
             for _ in range(NBOOT)]
    lo, hi = np.percentile(boots, [2.5, 97.5])
    rms = math.sqrt(float(np.mean(k**2)))

    rows = [dict(celltype=c, mean_fraction=round(float(fr), 5),
                 kappa=round(float(v), 3))
            for c, fr, v in sorted(zip(ct, frac, k), key=lambda t: -t[2])]
    rows.append(dict(celltype="__RMS_across_cell_types__", mean_fraction="",
                     kappa=round(rms, 3)))
    rows.append(dict(celltype="__RMS_bootstrap_95CI__", mean_fraction="",
                     kappa=f"{lo:.3f}-{hi:.3f}"))
    rows.append(dict(celltype="__median_across_cell_types__", mean_fraction="",
                     kappa=round(float(np.median(k)), 3)))

    # Design cost: nuclei per group needed to resolve one percentage point for
    # a cell type at 10% abundance
    def need(kappa: float, target_pp: float = 1.0, p: float = 0.10) -> float:
        return 2 * (1.96 * kappa * math.sqrt(p * (1 - p)) * 100 / target_pp) ** 2

    # Design cost uses the paper-wide kappa (b7's per-pair pooled value), not
    # the per-cell-type RMS
    pooled_kappa = _kappa(*_zn(COHORT_PAIR_TABLES[0][1], None)[:2])
    for lab, kv in (("design_nuclei_at_pooled_kappa", pooled_kappa),
                    ("design_nuclei_at_RMS", rms),
                    ("design_nuclei_at_median", float(np.median(k))),
                    ("design_nuclei_at_max", float(k.max())),
                    ("design_nuclei_multinomial", 1.0)):
        rows.append(dict(celltype=f"__{lab}__", mean_fraction="",
                         kappa=int(round(need(kv)))))
    # Correlation of kappa with abundance
    r = float(np.corrcoef(np.log(frac), k)[0, 1])
    rows.append(dict(celltype="__pearson_r_kappa_vs_log_abundance__",
                     mean_fraction="", kappa=round(r, 3)))
    return rows


# ----------------------------------------------------------------------------- B6
def b6_exponent_difference() -> List[Dict[str, object]]:
    """Donor-clustered bootstrap of the between-cohort difference in composition
    exponent."""
    path = f"{DR}/p5_two_layer_boot.csv"
    if not os.path.exists(path):
        print(f"  [skip] {path} 不存在，B6 跳过")
        return []
    rows = read_csv(path)
    cols = list(rows[0].keys())
    print(f"  p5_two_layer_boot.csv 列: {cols}")
    return [dict(note="inspect_columns", columns=";".join(cols), n_rows=len(rows))]


# ----------------------------------------------------------------------------- B6 / B7
COHORT_PAIR_TABLES = [
    ("Motor cortex", f"{ST}/Supplementary_Table_2a_composition_pairs_motor_cortex.csv", None),
    ("Seattle atlas", f"{ST}/Supplementary_Table_2b_composition_pairs_seattle.csv", None),
    ("PsychAD MSSM", f"{ST}/Supplementary_Table_9_psychad_composition_pairs.csv", "PsychAD MSSM"),
    ("PsychAD RADC", f"{ST}/Supplementary_Table_9_psychad_composition_pairs.csv", "PsychAD RADC"),
]


def _zn(path: str, cohort: str | None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = read_csv(path)
    if cohort is not None:
        rows = [x for x in rows if x["cohort"] == cohort]
    z = np.array([float(x["z"]) for x in rows])
    n = np.array([float(x["n_eff"]) for x in rows])
    d = np.array([x["donor"] for x in rows])
    return z, n, d


def _kappa(z: np.ndarray, n: np.ndarray) -> float:
    return float(math.sqrt(((z / np.sqrt(2.0 / n)) ** 2).mean()))


def b7_cohort_overdispersion() -> List[Dict[str, object]]:
    """Recompute kappa for all four cohorts from the deposited per-pair tables.

    The figure scripts read this rather than carrying literals, so a figure
    cannot drift away from the tables it claims to show."""
    rng = np.random.default_rng(SEED)
    out = []
    for label, path, cohort in COHORT_PAIR_TABLES:
        z, n, d = _zn(path, cohort)
        k = _kappa(z, n)
        uniq = np.unique(d)
        bs = []
        for _ in range(1000):
            pick = rng.choice(uniq, len(uniq), replace=True)
            idx = np.concatenate([np.flatnonzero(d == g) for g in pick])
            bs.append(_kappa(z[idx], n[idx]))
        lo, hi = np.percentile(bs, [2.5, 97.5])
        # Composition scaling exponent: log-log slope of |z| on n_eff
        b, _ = loglog_fit(n, np.abs(z) + 1e-12)
        out.append(dict(cohort=label, n_pairs=len(z), n_donors=len(uniq),
                        kappa=round(k, 4), kappa_lo=round(float(lo), 4),
                        kappa_hi=round(float(hi), 4), exponent=round(b, 4)))
    return out


def b6_cohort_exponent_difference() -> List[Dict[str, object]]:
    """Motor cortex against Seattle atlas: difference in composition exponent,
    donor-clustered bootstrap."""
    z1, n1, d1 = _zn(COHORT_PAIR_TABLES[0][1], None)
    z2, n2, d2 = _zn(COHORT_PAIR_TABLES[1][1], None)
    f = lambda z, n: loglog_fit(n, np.abs(z) + 1e-12)[0]
    obs = f(z1, n1) - f(z2, n2)
    rng = np.random.default_rng(SEED)
    u1, u2, diffs = np.unique(d1), np.unique(d2), []
    for _ in range(NBOOT):
        i1 = np.concatenate([np.flatnonzero(d1 == g) for g in rng.choice(u1, len(u1), True)])
        i2 = np.concatenate([np.flatnonzero(d2 == g) for g in rng.choice(u2, len(u2), True)])
        diffs.append(f(z1[i1], n1[i1]) - f(z2[i2], n2[i2]))
    diffs = np.array(diffs)
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    p = 2 * min((diffs <= 0).mean(), (diffs >= 0).mean())
    return [dict(comparison="motor_cortex_minus_seattle", difference=round(obs, 4),
                 lo=round(float(lo), 4), hi=round(float(hi), 4),
                 p_two_sided=round(float(p), 5), n_boot=NBOOT)]


def b8_depth_series() -> List[Dict[str, object]]:
    """The ATAC depth series, with seeds again treated as repeated measures
    rather than independent strata."""
    rows = read_csv(f"{DR}/nabec_L3_external_replication.csv")
    at150 = [x for x in rows if x["n"] == "150"]
    dep = [x for x in at150 if re.search(r"da\d+", x["file"])]
    base = [x for x in at150
            if not re.search(r"da\d+|fine", x["file"]) and x["celltype"] in ("Oligo", "ExN")]
    for x in base:
        x = x
    groups: Dict[int, List[Dict[str, str]]] = {5564: base}
    for x in dep:
        d = int(re.search(r"da(\d+)", x["file"]).group(1))
        groups.setdefault(d, []).append(x)

    out = []
    for depth in sorted(groups):
        at = groups[depth]
        or_all, se_all = _mh(_tabs_from(at))
        # The seed is encoded as seed<N> at the end of the filename
        per = []
        seeds = sorted({(re.search(r"seed(\d+)", x["file"]).group(1)) for x in at})
        for sd in seeds:
            sub = [x for x in at if re.search(r"seed(\d+)", x["file"]).group(1) == sd]
            o, se = _mh(_tabs_from(sub))
            per.append((math.log(o), se))
        lv = np.array([q[0] for q in per])
        within = float(np.mean([q[1] ** 2 for q in per]))
        between = float(np.var(lv, ddof=1)) if len(lv) > 1 else 0.0
        se_c = math.sqrt(within + between)
        mu = float(lv.mean())
        out.append(dict(depth_atac=depth, fold=round(depth / 5564, 2),
                        n_rows=len(at), n_seeds=len(per),
                        OR_seeds_as_strata=round(or_all, 4),
                        lo_seeds_as_strata=round(or_all * math.exp(-1.96 * se_all), 4),
                        hi_seeds_as_strata=round(or_all * math.exp(1.96 * se_all), 4),
                        OR_seeds_as_repeats=round(math.exp(mu), 4),
                        lo_seeds_as_repeats=round(math.exp(mu - 1.96 * se_c), 4),
                        hi_seeds_as_repeats=round(math.exp(mu + 1.96 * se_c), 4)))
    return out


def b9_kappa_by_celltype_all_cohorts() -> List[Dict[str, object]]:
    """Per-cell-type kappa in each of the four replicate sets.

    The paper originally reported one pooled kappa per cohort. But kappa varies
    far more between cell types WITHIN a dataset than it does BETWEEN datasets,
    and since design cost scales as kappa squared, that is the difference that
    actually determines a design.
    """
    out: List[Dict[str, object]] = []
    for label, path, cohort in COHORT_PAIR_TABLES:
        rows = read_csv(path)
        if cohort is not None:
            rows = [x for x in rows if x["cohort"] == cohort]
        ctk = "celltype" if "celltype" in rows[0] else "ct"
        groups: Dict[str, List[Tuple[float, float]]] = {}
        for x in rows:
            ratio = float(x["z"]) / math.sqrt(2.0 / float(x["n_eff"]))
            frac = float(x["p"]) if "p" in x and x["p"] != "" else float("nan")
            groups.setdefault(x[ctk], []).append((ratio, frac))
        ks, fr = [], []
        for ct, vals in sorted(groups.items()):
            if len(vals) < 5:
                continue
            r = np.array([a for a, _ in vals])
            k = math.sqrt(float((r ** 2).mean()))
            f = float(np.nanmean([b for _, b in vals]))
            ks.append(k); fr.append(f)
            out.append(dict(cohort=label, celltype=ct, n_pairs=len(vals),
                            mean_fraction=round(f, 6) if f == f else "",
                            kappa=round(k, 4)))
        ka = np.array(ks)
        corr = ""
        good = [(a, b) for a, b in zip(ks, fr) if b == b and b > 0]
        if len(good) > 2:
            corr = round(float(np.corrcoef(np.log([b for _, b in good]),
                                           [a for a, _ in good])[0, 1]), 3)
        out.append(dict(cohort=label, celltype="__summary__", n_pairs=len(ks),
                        mean_fraction="", kappa="",
                        kappa_min=round(float(ka.min()), 3),
                        kappa_max=round(float(ka.max()), 3),
                        spread=round(float(ka.max() / ka.min()), 2),
                        cost_spread=round(float((ka.max() / ka.min()) ** 2), 1),
                        r_kappa_vs_log_abundance=corr))
    return out


# ----------------------------------------------------------------------------- main
def write(path: str, rows: Sequence[Dict[str, object]]) -> None:
    if not rows:
        return
    keys: List[str] = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"  写出 {path}  ({len(rows)} 行)")


def main() -> None:
    print("B2 · 表达层指数支撑集敏感性 + SEA-AD 观测/null 配对比较")
    write(f"{DR}/p5_17_B2_exponent_support_sensitivity.csv",
          b2_exponent_support_sensitivity())
    print("B3 · Layer 3 阶梯的 seed 处理")
    write(f"{DR}/p5_17_B3_ladder_seed_handling.csv", b3_ladder_seed_handling())
    print("B5 · κ 的区间与逐细胞类型取值")
    write(f"{DR}/p5_17_B5_kappa_by_celltype.csv", b5_kappa())
    print("B7 · 四队列 κ 由逐对表重算（图脚本据此读取）")
    write(f"{DR}/p5_17_B7_cohort_overdispersion.csv", b7_cohort_overdispersion())
    print("B8 · ATAC 深度序列的 seed 处理")
    write(f"{DR}/p5_17_B8_depth_series.csv", b8_depth_series())
    print("B9 · 四个重复集的逐细胞类型 κ")
    write(f"{DR}/p5_17_B9_kappa_by_celltype_all_cohorts.csv", b9_kappa_by_celltype_all_cohorts())
    print("B6 · 组成层指数差异检验")
    write(f"{DR}/p5_17_B6_exponent_difference.csv", b6_cohort_exponent_difference())


if __name__ == "__main__":
    main()
