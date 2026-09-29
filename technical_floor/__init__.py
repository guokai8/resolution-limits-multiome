"""technical_floor —— 从重复文库测量单核多组学分析的分辨率下限。

三层各有一个下限，本包测量其中两层，并把它们换算成设计量：

    组成层   κ，过度离散倍数，逐细胞类型给出
    表达层   floor = a * n^b，以及由它得到的 |log2FC| 阈值
    设计     两张表，把上面两者换算成所需的供体数与细胞数

第三层（调控推断）需要原始的 peak x cell 矩阵，不在本包范围内；论文用的
promoter-enrichment 诊断另行提供。

用法：

    from technical_floor import (
        FloorConfig, read_composition, estimate_composition_floor,
    )

    df = read_composition("pairs.csv")
    floor = estimate_composition_floor(df, FloorConfig())
    print(floor.summary())

命令行：

    python3 -m technical_floor --composition pairs.csv --expression floors.csv --out results/
"""

from .composition import CompositionFloor, estimate_composition_floor
from .config import FloorConfig
from .design import (
    composition_design_table,
    composition_requirement,
    expression_design_table,
    nuclei_required,
    resolvable_difference,
)
from .expression import ExpressionFloor, estimate_expression_floor, floor_from_counts
from .io import apply_pair_criterion, read_composition, read_expression

__version__ = "1.0.0"

__all__ = [
    "FloorConfig",
    "CompositionFloor",
    "ExpressionFloor",
    "estimate_composition_floor",
    "estimate_expression_floor",
    "floor_from_counts",
    "read_composition",
    "read_expression",
    "apply_pair_criterion",
    "composition_design_table",
    "composition_requirement",
    "expression_design_table",
    "nuclei_required",
    "resolvable_difference",
    "__version__",
]
