"""technical_floor -- measure resolution floors from replicate libraries.

The paper describes three inference layers, each with its own floor. This
package measures two of them and converts both into design quantities:

    composition   kappa, the overdispersion factor, per cell type
    expression    floor = a * n^b, and the |log2FC| thresholds it implies
    design        two tables turning those into donors and nuclei required

The third layer, regulatory inference, needs the raw peak-by-cell matrices and
is out of scope here; the promoter-enrichment diagnostic the paper uses for it
lives in the analysis scripts instead.

Two of the three constants the paper reports do not transfer between datasets.
That is the reason this package exists: what travels is the procedure, not a
table of numbers, so the tool measures the floors in *your* data rather than
applying ours.

Library use:

    from technical_floor import (
        FloorConfig, read_composition, estimate_composition_floor,
    )

    df = read_composition("pairs.csv")
    floor = estimate_composition_floor(df, FloorConfig())
    print(floor.summary())

Command line:

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
