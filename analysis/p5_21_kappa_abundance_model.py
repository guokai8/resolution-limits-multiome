#!/usr/bin/env python3
"""审稿意见 V1：κ 为何在 √(p(1−p)) 标准化之后仍随丰度上升？

z = |f₁−f₂| / √(p(1−p)) 已经按多项抽样的标准差除过一次，所以 κ = z/√(2/n_eff)
本该与丰度无关。它却与 log 丰度强相关（正文报 r = 0.874）。这个脚本把"报一个
相关系数"换成在候选机制之间做判别。

把 κ 对丰度的标度写成 κ ∝ p^γ，三个候选机制给出不同的 γ：

  γ ≈ 0     多项抽样加一个与丰度无关的乘性膨胀 —— 正文隐含的模型
  γ ≈ −1/2  |f₁−f₂| 含一个与丰度无关的**绝对**误差项（例如注释噪声按核数计）
  γ ≈ +1/2  每个细胞类型的比例带一个恒定的**相对**误差（恒定变异系数）

γ 显著为正且接近 1/2，说明组成层噪声的行为像比例上的恒定 CV，而不像计数上的
多项噪声；那么真正可能在数据集之间迁移的常数是 CV 而不是 κ。脚本因此同时比较
κ 与 CV 在细胞类型之间的离散程度，并检验丰度效应是否只是神经元/非神经元类别
的混淆。

区间全部用供体聚类 bootstrap。

用法：
    python3 code/analysis/p5_21_kappa_abundance_model.py
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

NBOOT = 2000
SEED = 0
logger = logging.getLogger("v1")

SOURCES = {
    "Motor cortex": "p5_L1_scaling_pairs.csv",
    "PsychAD MSSM": "psychad_L1_MSSM_pairs.csv",
    "PsychAD RADC": "psychad_L1_RADC_pairs.csv",
}


def cell_class(name: str) -> str:
    """把细胞类型名归到大类，用来检验丰度效应是否只是类别混淆。"""
    n = name.lower()
    if n.startswith(("exc", "en_", "l2", "l3", "l4", "l5", "l6")):
        return "excitatory"
    if n.startswith(("inh", "in_", "lamp5", "pvalb", "sst", "vip", "sncg",
                     "pax6", "chandelier")):
        return "inhibitory"
    return "non-neuronal"


def ols(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    return np.linalg.lstsq(X, y, rcond=None)[0]


def boot_by_donor(df: pd.DataFrame, fn, nboot: int = NBOOT) -> Tuple[float, float, float]:
    """按供体重抽样，返回 (点估计, 2.5%, 97.5%)。"""
    point = fn(df)
    rng = np.random.default_rng(SEED)
    donors = df.donor.to_numpy()
    uniq = np.unique(donors)
    idx_by = {d: np.flatnonzero(donors == d) for d in uniq}
    draws: List[float] = []
    for _ in range(nboot):
        pick = rng.choice(uniq, size=len(uniq), replace=True)
        idx = np.concatenate([idx_by[d] for d in pick])
        try:
            v = fn(df.iloc[idx])
        except Exception:
            continue
        if np.isfinite(v):
            draws.append(v)
    if len(draws) < 50:
        return point, float("nan"), float("nan")
    return point, float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))


def gamma_simple(d: pd.DataFrame) -> float:
    X = np.column_stack([np.ones(len(d)), np.log(d.p.to_numpy())])
    return float(ols(X, np.log(d.kappa.to_numpy()))[1])


def gamma_with_class(d: pd.DataFrame) -> float:
    """加入细胞大类的固定效应后，丰度的斜率是否还在。"""
    cls = pd.get_dummies(d.cls, drop_first=True).to_numpy(dtype=float)
    X = np.column_stack([np.ones(len(d)), np.log(d.p.to_numpy()), cls])
    return float(ols(X, np.log(d.kappa.to_numpy()))[1])


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    root = Path(__file__).resolve().parents[1]
    dr = root / "data" / "derived_results"
    if not dr.exists():
        dr = Path("data/derived_results")

    out_rows: List[Dict[str, object]] = []
    print("κ ∝ p^γ 的标度指数（供体聚类 bootstrap 95% 区间）")
    print(f"{'队列':<14s} {'对数':>6s} {'类型':>5s} {'γ（仅丰度）':>22s} "
          f"{'γ（+细胞大类）':>22s}")
    for label, fname in SOURCES.items():
        path = dr / fname
        if not path.exists():
            logger.warning("  [skip] %s 不存在", path)
            continue
        d = pd.read_csv(path)
        d = d.rename(columns={"ct": "celltype"})
        d = d[(d.p > 0) & (d.z > 0) & (d.n_eff > 0)].copy()
        d["kappa"] = d.z / np.sqrt(2.0 / d.n_eff)
        d = d[d.kappa > 0]
        d["cls"] = d.celltype.map(cell_class)

        g0, g0lo, g0hi = boot_by_donor(d, gamma_simple)
        g1, g1lo, g1hi = boot_by_donor(d, gamma_with_class)
        print(f"{label:<14s} {len(d):>6d} {d.celltype.nunique():>5d} "
              f"{g0:>8.3f} ({g0lo:.3f}–{g0hi:.3f}) "
              f"{g1:>8.3f} ({g1lo:.3f}–{g1hi:.3f})")

        # 恒定相对误差模型的完整检验。
        # n_eff 是供体级的库总量（在供体内对所有细胞类型相同，且与 p 不相关），
        # 所以 κ² − 1 = c²·N·p/(1−p) 预言 log(κ²−1) 对 log(N·p/(1−p)) 的斜率为 1。
        per = d.groupby("celltype").apply(
            lambda s: pd.Series(dict(
                kap2=float(np.mean(s.z ** 2 / (2.0 / s.n_eff))),
                N=float(np.mean(s.n_eff)), p=float(np.mean(s.p)), npair=len(s),
            )), include_groups=False)
        use = per[per.kap2 > 1.02].copy()
        slope = float("nan")
        c_implied = float("nan")
        if len(use) >= 4:
            x = np.log(use.N * use.p / (1 - use.p))
            y = np.log(use.kap2 - 1)
            slope, inter = np.polyfit(x, y, 1)
            c_implied = float(np.exp(inter / 2))
        use["kappa"] = np.sqrt(use.kap2)
        use["c"] = np.sqrt((use.kap2 - 1) * (1 - use.p) / (use.N * use.p))

        out_rows.append(dict(
            cohort=label, n_pairs=len(d), n_celltypes=int(d.celltype.nunique()),
            gamma=g0, gamma_lo=g0lo, gamma_hi=g0hi,
            gamma_adj_class=g1, gamma_adj_class_lo=g1lo, gamma_adj_class_hi=g1hi,
            model_slope=float(slope), model_slope_theory=1.0,
            c_implied=c_implied, c_median=float(use.c.median()),
            n_celltypes_used=int(len(use)),
            kappa_spread=float(use.kappa.max() / use.kappa.min()),
            c_spread=float(use.c.max() / use.c.min()),
            sd_log_kappa=float(np.std(np.log(use.kappa), ddof=1)),
            sd_log_c=float(np.std(np.log(use.c), ddof=1)),
        ))

    print("\n恒定相对误差模型 κ² − 1 = c²·N·p/(1−p)：斜率应为 1，且 c 应比 κ 更齐")
    print(f"{'队列':<14s} {'斜率':>7s} {'类型':>5s} {'κ 极差':>8s} {'c 极差':>8s} "
          f"{'sd(log κ)':>10s} {'sd(log c)':>10s}")
    for r in out_rows:
        print(f"{r['cohort']:<14s} {r['model_slope']:>7.3f} {r['n_celltypes_used']:>5d} "
              f"{r['kappa_spread']:>7.2f}× {r['c_spread']:>7.2f}× "
              f"{r['sd_log_kappa']:>10.3f} {r['sd_log_c']:>10.3f}")

    print("\nγ 的参照值：0 = 多项 + 与丰度无关的乘性膨胀（正文隐含）；")
    print("             −0.5 = 绝对误差项；+0.5 = 恒定相对误差（恒定 CV）")
    print("结论读法：γ 显著为正只出现在 process 重复里；但恒定 CV 的完整预言被拒绝，")
    print("          且在 split-aliquot 队列里 c 比 κ 更离散，所以没有哪个重参数化")
    print("          能把下限变成与丰度无关的常数。")
    df = pd.DataFrame(out_rows)
    dest = dr / "p5_21_kappa_abundance_model.csv"
    df.to_csv(dest, index=False)
    print(f"\n写出 {dest}")


if __name__ == "__main__":
    main()
