"""The expression floor and how it scales.

For each (donor, cell type) pair,

    floor = median over genes of |log2(CPM_1g + 1) - log2(CPM_2g + 1)|

taken over the union of genes with CPM > 0 in either member, requiring at
least min_genes of them. The floors are then fitted in log-log space as

    floor = a * n_eff ** b

An exponent near -1/2 is the signature of pure multinomial sampling, so b is
not itself a finding. What carries information is a, the magnitude, and how far
the observations sit above a matched null.
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
    """The fitted expression floor.

    Attributes:
        exponent: the fitted b.
        exponent_lo: lower bound, bootstrap clustered by donor.
        exponent_hi: upper bound.
        coefficient: the fitted a.
        n_obs: (donor x cell type) observations in the fit.
        n_donors: donors in the fit.
        median_n_eff: median effective nucleus count observed. Extrapolating
            far past it is not supported by anything in the data, which is why
            the design table flags it.
        min_n_eff: the support floor this fit was run under.
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
        """The |log2 fold change| noise threshold at a given effective n."""
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
    """Compute one expression floor from a pair of pseudobulk count vectors.

    The two vectors must be aligned to the same gene order. The union rule
    keeps a gene when CPM > 0 in *either* member; an intersection rule would
    quietly drop exactly the genes that differ most between the libraries.

    Args:
        counts_a: per-gene counts for the first library.
        counts_b: per-gene counts for the second, same length and order.
        min_genes: minimum union size before a floor is returned.

    Returns:
        The floor, or None when the union is too small or a library has zero
        total counts. None rather than an exception: an unusable pair is an
        ordinary outcome when scanning many of them.

    Raises:
        ValueError: the two vectors differ in length.
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
    """Fit how the expression floor scales with effective nucleus count.

    As in the composition layer, the bootstrap resamples donors rather than
    observations. A resampled draw that happens to contain fewer than five
    observations is skipped rather than fitted on too little data.

    Args:
        df: long-form table with donor, celltype, n_eff, floor.
        config: the conventions to use; `min_n_eff` sets the support floor.

    Returns:
        An ExpressionFloor.

    Raises:
        ValueError: fewer than five observations survive the support floor.
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
