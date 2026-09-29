#!/usr/bin/env python3
"""
P5 步骤 17 · 审稿意见触发的三项重算（B2/B3/B5/B6）
==================================================
全部只用已沉积的派生表，不重跑上游扫描。输出写入 data/derived_results/。

B2  表达层标度指数对「最小有效核数」支撑集的敏感性，
    以及 SEA-AD 观测 vs 其自身 matched null 在**同一支撑集**上的配对比较。
B3  Layer 3 阶梯：seed 是同一批核的重抽样，不是独立 strata。
    改为「每个 seed 内按细胞类型做 MH 合并，再跨 seed 合成」。
B5  组成层 κ 的置信区间与逐细胞类型取值。
B6  两队列组成层指数差异的 donor-clustered bootstrap 检验。

用法：
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


# ----------------------------------------------------------------------------- 工具
def read_csv(path: str) -> List[Dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def loglog_fit(n: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    """返回 (指数 b, 系数 a)，拟合 log y ~ b log n + log a。"""
    b, la = np.polyfit(np.log(n), np.log(y), 1)
    return float(b), float(math.exp(la))


def boot_by_group(
    n: np.ndarray, y: np.ndarray, groups: np.ndarray, stat, nboot: int = NBOOT
) -> Tuple[float, float]:
    """按 group（供体）重抽样，返回 stat 的 95% 区间。"""
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
    """运动皮层观测指数随最小 n_eff 变化；SEA-AD 观测 vs 自身 null 的配对比较。"""
    rows: List[Dict[str, object]] = []

    # 50 核配对准则同样管表达层：只用 included_in_fit == yes 的行（26 个供体）
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

    # SEA-AD：观测与 null 逐行配对，必须在同一支撑集上比较
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
    """Mantel–Haenszel 合并 OR 与 Robins–Breslow–Greenland 的 log-OR 标准误。"""
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
    """由 link/tested 计数构造 2x2。

    全文统一口径（与 Table 3 表体一致）：
      · ha=False —— 阶梯与深度序列中没有任何一格为零，Haldane–Anscombe 校正不需要，
        加了反而把估计往 1 拉，并且与正文其余 OR 不可比。
      · exclude_detected=False —— 背景集不从检验集中扣除已检出的 link。
    两个开关保留，是为了能显式重现另一口径并量化差别。
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
    """seed 是同一批核的重抽样：先在 seed 内合并细胞类型，再跨 seed 合成。"""
    rows: List[Dict[str, object]] = []
    for cohort, path in (("primary", f"{DR}/primary_L3_nucleus_ladder.csv"),
                         ("external", f"{DR}/nabec_L3_nucleus_ladder.csv")):
        if not os.path.exists(path):
            print(f"  [skip] {path} 不存在")
            continue
        data = read_csv(path)
        for nlev in sorted({int(x["n"]) for x in data}):
            at = [x for x in data if int(x["n"]) == nlev]

            # 旧口径：所有 (细胞类型 x seed) 当作独立 strata
            or_all, se_all = _mh(_tabs_from(at))
            # 带 HA + 扣除检出的那一版，仅用于量化口径本身的贡献
            or_ha, _ = _mh(_tabs_from(at, ha=True, exclude_detected=True))

            # 新口径：每个 seed 内合并细胞类型
            per: List[Tuple[float, float]] = []
            for sd in sorted({x["seed"] for x in at}):
                sub = [x for x in at if x["seed"] == sd]
                o, se = _mh(_tabs_from(sub))
                per.append((math.log(o), se))
            lv = np.array([p[0] for p in per])
            within = float(np.mean([p[1] ** 2 for p in per]))
            between = float(np.var(lv, ddof=1)) if len(lv) > 1 else 0.0
            # 保守合成：种子内抽样方差 + 种子间下采样方差
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
    """κ 的逐细胞类型取值、RMS 汇总及其自助区间。"""
    src = read_csv(f"{DR}/p5_techrep_overdispersion.csv")
    # 已发表的 κ = 4.28 是 z_rms 的 RMS，不是 overdispersion_ratio 那一列的 RMS
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

    # 设计成本：10% 丰度的细胞类型，解析 1 个百分点所需每组核数
    def need(kappa: float, target_pp: float = 1.0, p: float = 0.10) -> float:
        return 2 * (1.96 * kappa * math.sqrt(p * (1 - p)) * 100 / target_pp) ** 2

    # 设计成本用全文统一的 κ（B7 的逐对合并值），不用逐细胞类型 RMS
    pooled_kappa = _kappa(*_zn(COHORT_PAIR_TABLES[0][1], None)[:2])
    for lab, kv in (("design_nuclei_at_pooled_kappa", pooled_kappa),
                    ("design_nuclei_at_RMS", rms),
                    ("design_nuclei_at_median", float(np.median(k))),
                    ("design_nuclei_at_max", float(k.max())),
                    ("design_nuclei_multinomial", 1.0)):
        rows.append(dict(celltype=f"__{lab}__", mean_fraction="",
                         kappa=int(round(need(kv)))))
    # κ 与丰度的相关
    r = float(np.corrcoef(np.log(frac), k)[0, 1])
    rows.append(dict(celltype="__pearson_r_kappa_vs_log_abundance__",
                     mean_fraction="", kappa=round(r, 3)))
    return rows


# ----------------------------------------------------------------------------- B6
def b6_exponent_difference() -> List[Dict[str, object]]:
    """两队列组成层指数差异的 donor-clustered bootstrap。"""
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
    """四个队列的 κ 全部由沉积的逐对表重算，供图脚本读取（取代写死的字面量）。"""
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
        # 组成层标度指数：|z| 对 n_eff 的对数-对数斜率
        b, _ = loglog_fit(n, np.abs(z) + 1e-12)
        out.append(dict(cohort=label, n_pairs=len(z), n_donors=len(uniq),
                        kappa=round(k, 4), kappa_lo=round(float(lo), 4),
                        kappa_hi=round(float(hi), 4), exponent=round(b, 4)))
    return out


def b6_cohort_exponent_difference() -> List[Dict[str, object]]:
    """运动皮层 vs Seattle 的组成层指数差异，供体聚类自助。"""
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
    """ATAC 深度序列：同样把 seed 当重复测量，而不是独立 strata。"""
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
        # seed 由文件名末尾的 seed<N> 标识
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
    """四个重复集各自的逐细胞类型 κ。

    论文原来只报每个队列一个合并 κ。但 κ 在数据集**内部**的跨细胞类型变异，
    远大于它在数据集**之间**的变异；由于成本按 κ² 走，这个差别决定了设计。
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
