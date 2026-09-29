"""表达层下限及其标度。

论文的量：对每个 (供体, 细胞类型) 配对，
    floor = median_g | log2(CPM_1g + 1) - log2(CPM_2g + 1) |
并集取 CPM > 0 的基因，要求至少 min_genes 个。随后在对数-对数空间拟合
    floor = a * n_eff ** b
b ≈ -1/2 是纯多项式抽样的签名，所以 b 本身不是发现；有信息的是 a，以及
观测相对匹配零模型的超出量。
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .config import FloorConfig

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class ExpressionFloor:
    """表达层下限拟合的结果。

    Attributes:
        exponent: 拟合指数 b。
        exponent_lo: 按供体聚类自助的下界。
        exponent_hi: 上界。
        coefficient: 拟合系数 a。
        n_obs: 参与拟合的 (供体 x 细胞类型) 观测数。
        n_donors: 参与拟合的供体数。
        median_n_eff: 观测的有效细胞数中位数；远超它的外推没有依据。
        min_n_eff: 本次拟合所用的支撑集下限。
    """

    exponent: float
    exponent_lo: float
    exponent_hi: float
    coefficient: float
    n_obs: int
    n_donors: int
    median_n_eff: float
    min_n_eff: float

    def threshold(self, n_eff: float) -> float:
        """给定每细胞类型有效细胞数，返回 |log2 fold change| 的技术噪声阈值。"""
        if n_eff <= 0:
            raise ValueError("n_eff 必须为正")
        return float(self.coefficient * n_eff ** self.exponent)

    def summary(self) -> str:
        return "\n".join([
            f"表达层下限 = {self.coefficient:.2f} × n^{self.exponent:.3f} "
            f"(指数 95% CI {self.exponent_lo:.3f}–{self.exponent_hi:.3f})",
            f"  基于 {self.n_obs} 个观测 / {self.n_donors} 个供体，"
            f"有效细胞数中位 {self.median_n_eff:.0f}",
            "  纯多项式抽样的指数是 -0.500；接近该值说明下限与抽样不可区分",
        ])


def floor_from_counts(
    counts_a: Sequence[float], counts_b: Sequence[float], min_genes: int = 200
) -> Optional[float]:
    """由一对伪批量计数向量算出表达下限。

    两个向量必须按同一基因顺序对齐。并集规则：任一成员 CPM > 0 即纳入。

    Args:
        counts_a: 第一个文库的每基因计数。
        counts_b: 第二个文库的每基因计数，与 counts_a 等长同序。
        min_genes: 返回结果所需的最少并集基因数。

    Returns:
        下限，或并集基因不足 / 某库总计数为零时返回 None。

    Raises:
        ValueError: 两个向量长度不等。
    """
    a = np.asarray(counts_a, dtype=float)
    b = np.asarray(counts_b, dtype=float)
    if a.shape != b.shape:
        raise ValueError(f"计数向量长度不等：{a.shape} 与 {b.shape}")
    ta, tb = a.sum(), b.sum()
    if ta <= 0 or tb <= 0:
        return None
    cpm_a = a / ta * 1e6
    cpm_b = b / tb * 1e6
    keep = (cpm_a > 0) | (cpm_b > 0)
    if int(keep.sum()) < min_genes:
        return None
    diff = np.abs(np.log2(cpm_a[keep] + 1.0) - np.log2(cpm_b[keep] + 1.0))
    return float(np.median(diff))


def _loglog_fit(n: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    b, log_a = np.polyfit(np.log(n), np.log(y), 1)
    return float(b), float(math.exp(log_a))


def estimate_expression_floor(
    df: pd.DataFrame, config: FloorConfig = FloorConfig()
) -> ExpressionFloor:
    """拟合表达层下限对有效细胞数的标度关系。

    Args:
        df: 长表，列为 donor, celltype, n_eff, floor。
        config: 口径配置；`min_n_eff` 设定支撑集下限。

    Returns:
        ExpressionFloor。

    Raises:
        ValueError: 应用支撑集下限后不足 5 个观测。
    """
    work = df[df.n_eff >= config.min_n_eff] if config.min_n_eff > 0 else df
    work = work[(work.n_eff > 0) & (work.floor > 0)]
    if len(work) < 5:
        raise ValueError(
            f"min_n_eff={config.min_n_eff} 之后只剩 {len(work)} 个观测，不足以拟合"
        )

    n = work.n_eff.to_numpy()
    y = work.floor.to_numpy()
    donors = work.donor.to_numpy()
    exponent, coefficient = _loglog_fit(n, y)

    rng = np.random.default_rng(config.seed)
    uniq = np.unique(donors)
    draws: List[float] = []
    index_by_donor = {d: np.flatnonzero(donors == d) for d in uniq}
    for _ in range(config.n_bootstrap):
        pick = rng.choice(uniq, size=len(uniq), replace=True)
        idx = np.concatenate([index_by_donor[d] for d in pick])
        if idx.size < 5:
            continue
        draws.append(_loglog_fit(n[idx], y[idx])[0])
    lo, hi = (float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))) \
        if draws else (float("nan"), float("nan"))

    logger.info("表达层指数 %.4f，系数 %.3f（%d 观测）", exponent, coefficient, len(work))
    return ExpressionFloor(
        exponent=exponent, exponent_lo=lo, exponent_hi=hi, coefficient=coefficient,
        n_obs=len(work), n_donors=int(len(uniq)),
        median_n_eff=float(np.median(n)), min_n_eff=config.min_n_eff,
    )
