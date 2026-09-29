"""组成层下限：过度离散系数 κ。

论文的量：对每个 (供体, 细胞类型) 配对，
    z       = |f1 - f2| / sqrt(p(1-p))
    期望    = sqrt(2 / n_eff),  n_eff = 2 / (1/n1 + 1/n2)
    κ       = RMS(z / 期望)

用 RMS 而不是中位数，是因为它聚合的是标准差之比。κ 随细胞类型变化很大，所以
本模块既给合并值也给逐类型值——论文只报了前者，而后者才是设计时真正要用的。
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from .config import FloorConfig

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class CompositionFloor:
    """组成层下限的估计结果。

    Attributes:
        kappa: 合并的过度离散倍数（全部配对的 RMS）。
        kappa_lo: 自助下界。
        kappa_hi: 自助上界。
        per_celltype: 逐细胞类型的 κ 与平均丰度，按 κ 降序。
        n_pairs: 参与估计的配对数。
        n_donors: 参与估计的供体数。
        abundance_r: κ 与 log 丰度的皮尔逊相关；正值表示丰富类型更离散。
    """

    kappa: float
    kappa_lo: float
    kappa_hi: float
    per_celltype: pd.DataFrame
    n_pairs: int
    n_donors: int
    abundance_r: Optional[float] = None
    _notes: List[str] = field(default_factory=list)

    def summary(self) -> str:
        lines = [
            f"组成层过度离散 κ = {self.kappa:.3f} "
            f"(95% CI {self.kappa_lo:.3f}–{self.kappa_hi:.3f})",
            f"  基于 {self.n_pairs} 个配对 / {self.n_donors} 个供体",
        ]
        if len(self.per_celltype):
            lo = self.per_celltype.kappa.min()
            hi = self.per_celltype.kappa.max()
            lines.append(f"  逐细胞类型 κ 从 {lo:.2f} 到 {hi:.2f}（{hi/lo:.1f} 倍）")
        if self.abundance_r is not None:
            lines.append(f"  κ 与 log 丰度的相关 r = {self.abundance_r:+.2f}")
        return "\n".join(lines)


def _zscores(df: pd.DataFrame) -> pd.DataFrame:
    """由细胞计数算出 z、n_eff 与 z / 期望之比。"""
    out = df.copy()
    tot1 = out.groupby("donor").n1.transform("sum")
    tot2 = out.groupby("donor").n2.transform("sum")
    if (tot1 <= 0).any() or (tot2 <= 0).any():
        raise ValueError("有供体的某个文库细胞总数为零")
    f1 = out.n1 / tot1
    f2 = out.n2 / tot2
    p = (f1 + f2) / 2.0
    out["n_eff"] = 2.0 / (1.0 / out.n1.clip(lower=1) + 1.0 / out.n2.clip(lower=1))
    denom = np.sqrt(p * (1.0 - p))
    # p 为 0 或 1 时该细胞类型在两库都不存在（或占满），比值无定义
    out["z"] = np.where(denom > 0, (f1 - f2).abs() / denom, np.nan)
    out["expected"] = np.sqrt(2.0 / out.n_eff)
    out["ratio"] = out.z / out.expected
    out["p_mean"] = p
    return out


def _rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(values))))


def estimate_composition_floor(
    df: pd.DataFrame, config: FloorConfig = FloorConfig()
) -> CompositionFloor:
    """由重复配对的细胞计数估计组成层下限。

    Args:
        df: 长表，列为 donor, celltype, n1, n2。
        config: 口径配置。

    Returns:
        CompositionFloor。

    Raises:
        ValueError: 执行准则后没有可用配对。
    """
    if {"n1", "n2"}.issubset(df.columns):
        work = _zscores(df).dropna(subset=["ratio"])
    elif {"z", "n_eff"}.issubset(df.columns):
        # 预算好的形式：直接用沉积的 z 与 n_eff，不重建占比
        work = df.copy()
        work["expected"] = np.sqrt(2.0 / work.n_eff)
        work["ratio"] = work.z / work.expected
        work["p_mean"] = work["p"] if "p" in work.columns else np.nan
        work = work.dropna(subset=["ratio"])
    else:
        raise ValueError("需要 n1/n2 或 z/n_eff")
    if not len(work):
        raise ValueError("没有可用配对：所有细胞类型的合并占比为 0 或 1")

    kappa = _rms(work.ratio.to_numpy())

    rng = np.random.default_rng(config.seed)
    donors = work.donor.unique()
    draws: List[float] = []
    by_donor: Dict[str, np.ndarray] = {
        d: g.ratio.to_numpy() for d, g in work.groupby("donor")
    }
    for _ in range(config.n_bootstrap):
        pick = rng.choice(donors, size=len(donors), replace=True)
        draws.append(_rms(np.concatenate([by_donor[d] for d in pick])))
    lo, hi = (float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))) \
        if draws else (float("nan"), float("nan"))

    per = (
        work.groupby("celltype")
        .apply(lambda g: pd.Series({
            "n_pairs": len(g),
            "mean_fraction": float(g.p_mean.mean()),
            "kappa": _rms(g.ratio.to_numpy()),
        }), include_groups=False)
        .reset_index()
        .sort_values("kappa", ascending=False)
        .reset_index(drop=True)
    )

    r: Optional[float] = None
    if len(per) > 2 and (per.mean_fraction > 0).all():
        r = float(np.corrcoef(np.log(per.mean_fraction), per.kappa)[0, 1])

    logger.info("κ = %.3f（%d 配对 / %d 供体）", kappa, len(work), work.donor.nunique())
    return CompositionFloor(
        kappa=kappa, kappa_lo=lo, kappa_hi=hi, per_celltype=per,
        n_pairs=len(work), n_donors=int(work.donor.nunique()), abundance_r=r,
    )
