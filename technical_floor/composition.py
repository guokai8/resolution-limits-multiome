"""The composition floor: the overdispersion factor kappa.

For each (donor, cell type) pair,

    z        = |f1 - f2| / sqrt(p(1-p))
    expected = sqrt(2 / n_eff),  n_eff = 2 / (1/n1 + 1/n2)
    kappa    = RMS(z / expected)

RMS rather than a median, because what is being aggregated is a ratio of
standard deviations and the squares are what add.

kappa varies a great deal between cell types, so this module reports both the
pooled value and the per-type values. The paper reports the pooled one; the
per-type ones are what a design calculation actually needs, since cost scales
as kappa squared and the spread within a dataset exceeds the spread between
cohorts.
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
    """The estimated composition floor.

    Attributes:
        kappa: the pooled overdispersion factor, the RMS over all pairs.
        kappa_lo: lower bootstrap bound.
        kappa_hi: upper bootstrap bound.
        per_celltype: kappa and mean abundance per cell type, kappa descending.
        n_pairs: pairs behind the estimate.
        n_donors: donors behind the estimate.
        abundance_r: Pearson correlation of kappa with log abundance. Positive
            means commoner cell types are more overdispersed, which is the
            pattern the paper reports and does not fully explain.
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
    """Turn nucleus counts into z, n_eff and the ratio z / expected."""
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
    # p of 0 or 1 means the cell type is absent from both libraries, or fills
    # them, and the ratio is undefined. NaN here, dropped by the caller.
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
    """Estimate the composition floor from replicate-pair nucleus counts.

    The bootstrap resamples donors, not pairs: pairs from one donor share a
    dissociation and a library prep, so resampling pairs would understate the
    interval.

    Args:
        df: long-form table with donor, celltype, n1, n2 -- or with donor,
            celltype, z, n_eff if the statistics are already computed.
        config: the conventions to use.

    Returns:
        A CompositionFloor.

    Raises:
        ValueError: neither column set is present, or no usable pairs remain.
    """
    if {"n1", "n2"}.issubset(df.columns):
        work = _zscores(df).dropna(subset=["ratio"])
    elif {"z", "n_eff"}.issubset(df.columns):
        # Precomputed form: use the deposited z and n_eff as they stand rather
        # than rebuilding proportions, whose denominator depends on which cell
        # types the original analysis included.
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
