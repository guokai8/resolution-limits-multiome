"""Turn the measured floors into quantities an experimenter can plan on.

Composition: for a cell type of abundance p, the smallest resolvable
between-group difference is

    F = z * kappa * sqrt(p(1-p)) * sqrt(2 / (D*N))

Solving for N instead shows the cost, which grows as kappa squared -- so a
twofold error in kappa is a fourfold error in the nuclei you need to collect.

Expression: the threshold comes straight from the fitted power law.

The point of both tables is to make a floor usable before the data are
collected, rather than a caution issued afterwards.
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
    """The smallest resolvable between-group difference, in percentage points."""
    if donors <= 0 or nuclei_per_donor <= 0:
        raise ValueError("供体数与每供体细胞数都必须为正")
    p = config.abundance
    total = donors * nuclei_per_donor
    return float(100.0 * config.z * kappa * math.sqrt(p * (1 - p)) * math.sqrt(2.0 / total))


def nuclei_required(
    kappa: float, config: FloorConfig = FloorConfig(), target_pp: Optional[float] = None
) -> int:
    """Nuclei needed per group to resolve a difference of target_pp points."""
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
    """Resolvable difference over a grid of donors x nuclei per donor.

    Four columns, in increasing order of realism: the multinomial bound
    (kappa = 1), the pooled kappa, and the median and worst cell type. The last
    two are the ones a pooled kappa hides -- within one dataset they can differ
    by more than the two cohorts differ from each other.
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
    """Read the other way: nuclei per group needed to resolve target_pp."""
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
    """|log2 fold change| noise thresholds over a ladder of effective n.

    `extrapolated` flags any rung more than tenfold past the median effective n
    actually observed, where the power law is no longer backed by data.
    """
    return pd.DataFrame([
        {
            "n_eff": n,
            "log2fc_threshold": round(floor.threshold(n), 3),
            "extrapolated": bool(n > floor.median_n_eff * 10),
        }
        for n in ladder
    ])
