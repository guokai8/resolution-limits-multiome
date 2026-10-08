"""Reading and validating the input tables.

The tool accepts tidy long-form CSV only. That is the format the paper
deposited, and it is the format easiest to export from any single-cell
pipeline. Each of the two input kinds has its own required column set.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Sequence

import pandas as pd

logger = logging.getLogger(__name__)

COMPOSITION_COUNT_COLUMNS: Sequence[str] = ("donor", "celltype", "n1", "n2")
COMPOSITION_ZSCORE_COLUMNS: Sequence[str] = ("donor", "celltype", "z", "n_eff")
EXPRESSION_COLUMNS: Sequence[str] = ("donor", "celltype", "n_eff", "floor")
# n1/n2 are optional in the expression input. Supplying them lets the same
# pairing criterion run on the expression layer too.


def _require_columns(df: pd.DataFrame, required: Sequence[str], path: Path) -> None:
    """Raise with the full column list, so a typo is obvious from the message."""
    missing: List[str] = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"{path} 缺少必需列 {missing}；需要 {list(required)}，实际有 {list(df.columns)}"
        )


def read_composition(path: Path | str) -> pd.DataFrame:
    """Read the composition-layer input.

    Two forms are accepted; either one is enough.

      * Raw counts -- donor, celltype, n1, n2. The natural export from your own
        pipeline.
      * Precomputed statistics -- donor, celltype, z, n_eff. The right choice
        for reproducing a published analysis, because the denominator of a
        proportion depends on which cell types entered that analysis, and
        rebuilding it from counts will not always reproduce it.

    When both are present the raw counts win.

    Args:
        path: path to the long-form CSV.

    Returns:
        The validated DataFrame.

    Raises:
        ValueError: a required column is missing, or a count is negative.
    """
    path = Path(path)
    df = pd.read_csv(path)
    has_counts = all(c in df.columns for c in COMPOSITION_COUNT_COLUMNS)
    has_z = all(c in df.columns for c in COMPOSITION_ZSCORE_COLUMNS)
    if not (has_counts or has_z):
        raise ValueError(
            f"{path} 必须提供 {list(COMPOSITION_COUNT_COLUMNS)} 或 "
            f"{list(COMPOSITION_ZSCORE_COLUMNS)}；实际有 {list(df.columns)}"
        )
    if has_counts and (df[["n1", "n2"]] < 0).any().any():
        raise ValueError(f"{path} 含负细胞计数")
    if has_z and (df.n_eff <= 0).any():
        raise ValueError(f"{path} 含非正 n_eff")
    logger.info("组成层输入：%d 行，%d 个供体，%d 个细胞类型",
                len(df), df.donor.nunique(), df.celltype.nunique())
    return df


def read_expression(path: Path | str) -> pd.DataFrame:
    """Read the expression-layer input.

    One row per (donor, cell type) pair, carrying a floor that has already been
    computed and the effective nucleus count behind it.
    `technical_floor.expression.floor_from_counts` builds this table from
    pseudobulk counts.

    Args:
        path: path to the long-form CSV; required columns are donor, celltype,
            n_eff and floor.

    Returns:
        The validated DataFrame.

    Raises:
        ValueError: a required column is missing, or n_eff or floor is not
            positive -- both are logged, so neither can be zero or negative.
    """
    path = Path(path)
    df = pd.read_csv(path)
    _require_columns(df, EXPRESSION_COLUMNS, path)
    bad = (df.n_eff <= 0) | (df.floor <= 0)
    if bad.any():
        raise ValueError(
            f"{path} 有 {int(bad.sum())} 行的 n_eff 或 floor 非正；对数-对数拟合无法使用"
        )
    logger.info("表达层输入：%d 行，%d 个供体", len(df), df.donor.nunique())
    return df


def apply_pair_criterion(df: pd.DataFrame, min_nuclei: int) -> pd.DataFrame:
    """Apply the donor-level minimum-nucleus criterion.

    The criterion is applied to the donor total across cell types, not per cell
    type: a donor qualifies when both of its libraries reach the threshold
    overall. A single small cell type should not disqualify a donor.

    It governs the expression layer as well as the composition layer. In the
    paper's own analysis it was at one point applied only to composition, which
    is why it lives in one function used by both.

    Args:
        df: long-form table with donor, n1 and n2.
        min_nuclei: minimum nucleus count per library.

    Returns:
        A copy carrying an `included_in_fit` column. Rows are flagged rather
        than dropped, so a caller can report what was excluded.
    """
    if "n1" not in df.columns or "n2" not in df.columns:
        raise ValueError("执行配对准则需要 n1 与 n2 两列")
    totals = df.groupby("donor")[["n1", "n2"]].sum()
    keep = totals[(totals.n1 >= min_nuclei) & (totals.n2 >= min_nuclei)].index
    dropped = sorted(set(totals.index) - set(keep))
    if dropped:
        logger.warning("按 >=%d 核准则排除 %d 个供体：%s", min_nuclei, len(dropped), dropped)
    out = df.copy()
    out["included_in_fit"] = out.donor.isin(keep)
    return out
