"""估计量的不可变配置。

论文里每一个口径选择在这里都有一个显式字段，默认值就是论文所用的那一套。
改动任何一个字段都会改变结果，所以它们是配置而不是常量。
"""

from __future__ import annotations

from dataclasses import dataclass

# 论文的组成层统计量：z = |f1 - f2| / sqrt(p(1-p))，与多项式期望 sqrt(2/n_eff) 相比。
# n_eff 取两库细胞数的调和平均。
DEFAULT_MIN_NUCLEI_PER_PAIR = 50
DEFAULT_MIN_GENES = 200
# 1000 是论文沉积结果所用的次数；改动它会在第三位小数上移动区间，
# 从而让工具与论文报告的区间不再逐位一致。
DEFAULT_BOOTSTRAP = 1000
DEFAULT_SEED = 0


@dataclass(frozen=True)
class FloorConfig:
    """一次下限估计的全部口径选择。

    Attributes:
        min_nuclei_per_pair: 供体级入选门槛。两个文库都必须达到该细胞数，
            论文用 50。这条准则同时管组成层与表达层。
        min_genes: 计算表达下限所需的最少并集基因数，论文用 200。
        min_n_eff: 进入标度拟合的最小有效细胞数。0 表示不设下限（论文的已发表
            口径）；设为 10 或 25 可检验结论对低计数尾部的稳健性。
        n_bootstrap: 自助重抽样次数，按供体聚类。
        seed: 自助的随机种子。
        abundance: 设计表所针对的细胞类型丰度，论文用 0.10。
        target_pp: 设计表所要分辨的组间差异，单位百分点，论文用 1.0。
        z: 双侧显著性对应的正态分位数，论文用 1.96。
    """

    min_nuclei_per_pair: int = DEFAULT_MIN_NUCLEI_PER_PAIR
    min_genes: int = DEFAULT_MIN_GENES
    min_n_eff: float = 0.0
    n_bootstrap: int = DEFAULT_BOOTSTRAP
    seed: int = DEFAULT_SEED
    abundance: float = 0.10
    target_pp: float = 1.0
    z: float = 1.96

    def __post_init__(self) -> None:
        if not 0.0 < self.abundance < 1.0:
            raise ValueError(f"abundance 必须在 (0,1) 内，收到 {self.abundance}")
        if self.target_pp <= 0:
            raise ValueError(f"target_pp 必须为正，收到 {self.target_pp}")
        if self.n_bootstrap < 0:
            raise ValueError(f"n_bootstrap 不能为负，收到 {self.n_bootstrap}")
