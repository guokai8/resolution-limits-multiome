"""输入表的读取与校验。

工具只接受长表（tidy CSV），因为那是论文自己沉积的格式，也是使用者最容易
从任意单细胞流程里导出的格式。两类输入各有一个必需列集合。
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
# n1/n2 可选；给了就能在表达层上执行同一条配对准则


def _require_columns(df: pd.DataFrame, required: Sequence[str], path: Path) -> None:
    missing: List[str] = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"{path} 缺少必需列 {missing}；需要 {list(required)}，实际有 {list(df.columns)}"
        )


def read_composition(path: Path | str) -> pd.DataFrame:
    """读取组成层输入。

    接受两种形式，二选一即可：

      · 原始计数 —— donor, celltype, n1, n2。适合从自己的流程直接导出。
      · 预算好的统计量 —— donor, celltype, z, n_eff。适合复现已发表的分析，
        因为占比的分母取决于哪些细胞类型进入了分析，由计数重建未必一致。

    两者都有时以原始计数为准。

    Args:
        path: 长表 CSV 路径。

    Returns:
        校验后的 DataFrame。

    Raises:
        ValueError: 缺列，或计数为负。
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
    """读取表达层输入。

    每行是一个 (供体, 细胞类型) 配对已算好的下限与有效细胞数。用
    `technical_floor.expression.floor_from_counts` 可以从伪批量计数生成这张表。

    Args:
        path: 长表 CSV 路径，必需列 donor, celltype, n_eff, floor。

    Returns:
        校验后的 DataFrame。

    Raises:
        ValueError: 缺列，或 n_eff/floor 非正。
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
    """执行供体级的最小细胞数准则。

    论文的准则是按供体汇总后两个文库都要达到门槛，而不是按细胞类型。这条准则
    同时管组成层和表达层——在论文自己的分析里它一度只落到了组成层上。

    Args:
        df: 含 donor, n1, n2 的长表。
        min_nuclei: 每个文库的最小细胞数。

    Returns:
        只含合格供体的副本，并带一列 `included_in_fit`。
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
