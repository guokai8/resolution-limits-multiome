"""把测得的下限换算成设计量。

组成层：丰度为 p 的细胞类型，可分辨的最小组间差异是
    F = z * kappa * sqrt(p(1-p)) * sqrt(2 / (D*N))
反过来解 N，代价随 kappa 的平方增长。

表达层：阈值直接由拟合的幂律给出。

两张表的作用是让下限在采数之前可用，而不是事后的告诫。
"""

from __future__ import annotations

import logging
import math
from typing import Iterable, Optional, Sequence

import pandas as pd

from .composition import CompositionFloor
from .config import FloorConfig
from .expression import ExpressionFloor

logger = logging.getLogger(__name__)

DEFAULT_DESIGNS: Sequence[tuple[int, int]] = ((10, 2000), (20, 5000), (50, 10000), (200, 20000))
DEFAULT_NUCLEI_LADDER: Sequence[int] = (25, 50, 100, 200, 500, 1000)


def resolvable_difference(
    kappa: float, donors: int, nuclei_per_donor: int, config: FloorConfig = FloorConfig()
) -> float:
    """返回可分辨的最小组间差异，单位百分点。"""
    if donors <= 0 or nuclei_per_donor <= 0:
        raise ValueError("供体数与每供体细胞数都必须为正")
    p = config.abundance
    total = donors * nuclei_per_donor
    return float(100.0 * config.z * kappa * math.sqrt(p * (1 - p)) * math.sqrt(2.0 / total))


def nuclei_required(
    kappa: float, config: FloorConfig = FloorConfig(), target_pp: Optional[float] = None
) -> int:
    """返回每组所需细胞数，以分辨 target_pp 个百分点的差异。"""
    p = config.abundance
    target = config.target_pp if target_pp is None else target_pp
    if target <= 0:
        raise ValueError("目标差异必须为正")
    return int(round(2.0 * (config.z * kappa * math.sqrt(p * (1 - p)) * 100.0 / target) ** 2))


def composition_design_table(
    floor: CompositionFloor,
    config: FloorConfig = FloorConfig(),
    designs: Iterable[tuple[int, int]] = DEFAULT_DESIGNS,
) -> pd.DataFrame:
    """按供体数 x 每供体细胞数，给出可分辨差异。

    同时给出多项式界（kappa = 1）、合并 kappa，以及逐细胞类型的中位与最差
    kappa——最后两列是论文只报合并值时看不到的东西。
    """
    per = floor.per_celltype
    k_median = float(per.kappa.median()) if len(per) else floor.kappa
    k_worst = float(per.kappa.max()) if len(per) else floor.kappa
    rows = []
    for donors, nuclei in designs:
        rows.append({
            "donors": donors,
            "nuclei_per_donor": nuclei,
            "multinomial_pp": round(resolvable_difference(1.0, donors, nuclei, config), 3),
            "pooled_kappa_pp": round(resolvable_difference(floor.kappa, donors, nuclei, config), 3),
            "median_celltype_pp": round(resolvable_difference(k_median, donors, nuclei, config), 3),
            "worst_celltype_pp": round(resolvable_difference(k_worst, donors, nuclei, config), 3),
        })
    return pd.DataFrame(rows)


def composition_requirement(
    floor: CompositionFloor, config: FloorConfig = FloorConfig()
) -> pd.DataFrame:
    """反向读法：分辨 target_pp 所需的每组细胞数。"""
    per = floor.per_celltype
    entries = [("multinomial bound", 1.0), ("pooled kappa", floor.kappa)]
    if len(per):
        entries += [
            ("median cell type", float(per.kappa.median())),
            ("most overdispersed cell type", float(per.kappa.max())),
        ]
    return pd.DataFrame([
        {"basis": label, "kappa": round(k, 3), "nuclei_per_group": nuclei_required(k, config)}
        for label, k in entries
    ])


def expression_design_table(
    floor: ExpressionFloor, ladder: Sequence[int] = DEFAULT_NUCLEI_LADDER
) -> pd.DataFrame:
    """按每细胞类型有效细胞数，给出 |log2 fold change| 的技术噪声阈值。"""
    return pd.DataFrame([
        {
            "n_eff": n,
            "log2fc_threshold": round(floor.threshold(n), 3),
            "extrapolated": bool(n > floor.median_n_eff * 10),
        }
        for n in ladder
    ])
