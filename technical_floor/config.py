"""Immutable configuration for the floor estimators.

Every convention the paper had to choose appears here as an explicit field, and
the defaults are the ones the paper used. Changing any field changes the result,
which is why these are configuration rather than constants buried in the code.
"""

from __future__ import annotations

from dataclasses import dataclass

# The composition statistic is z = |f1 - f2| / sqrt(p(1-p)), compared against
# the multinomial expectation sqrt(2/n_eff), where n_eff is the harmonic mean of
# the two libraries' nucleus counts.
DEFAULT_MIN_NUCLEI_PER_PAIR = 50
DEFAULT_MIN_GENES = 200
# 1000 is the number of resamples behind the deposited results. Changing it
# moves the interval in the third decimal, so the tool would stop agreeing
# digit for digit with the intervals the paper reports.
DEFAULT_BOOTSTRAP = 1000
DEFAULT_SEED = 0


@dataclass(frozen=True)
class FloorConfig:
    """Every convention behind one floor estimate.

    Frozen on purpose: a config is passed down through the estimators, and a
    field that could be mutated halfway would make a result impossible to
    attribute to a stated set of choices.

    Attributes:
        min_nuclei_per_pair: donor-level inclusion threshold. Both libraries
            must reach this nucleus count; the paper used 50. The criterion
            governs the expression layer as well as the composition layer.
        min_genes: minimum size of the union gene set needed to compute an
            expression floor; the paper used 200.
        min_n_eff: minimum effective nucleus count for a pair to enter the
            scaling fit. 0 means no floor, which is the published convention;
            10 or 25 tests whether the conclusion survives dropping the
            low-count tail.
        n_bootstrap: resamples, clustered by donor rather than by pair, since
            pairs from one donor are not independent.
        seed: random seed for the bootstrap.
        abundance: the cell-type abundance the design tables are computed at;
            the paper used 0.10.
        target_pp: the between-group difference the design tables resolve, in
            percentage points; the paper used 1.0.
        z: normal quantile for the two-sided interval; the paper used 1.96.
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
        # Validate here rather than at the point of use: a bad abundance would
        # otherwise surface as a silently wrong design table.
        if not 0.0 < self.abundance < 1.0:
            raise ValueError(f"abundance 必须在 (0,1) 内，收到 {self.abundance}")
        if self.target_pp <= 0:
            raise ValueError(f"target_pp 必须为正，收到 {self.target_pp}")
        if self.n_bootstrap < 0:
            raise ValueError(f"n_bootstrap 不能为负，收到 {self.n_bootstrap}")
