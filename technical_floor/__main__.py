"""Command-line entry point: python3 -m technical_floor

Reads one or both input tables, applies the pairing criterion, runs the
estimators and writes the design tables. Everything printed here is also
written to CSV, so the console output is a summary rather than the deliverable.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from .composition import estimate_composition_floor
from .config import DEFAULT_BOOTSTRAP, DEFAULT_SEED, FloorConfig
from .design import (
    composition_design_table,
    composition_requirement,
    expression_design_table,
)
from .expression import estimate_expression_floor
from .io import apply_pair_criterion, read_composition, read_expression

logger = logging.getLogger("technical_floor")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="technical_floor",
        description="Measure composition and expression resolution floors from replicate libraries.",
    )
    p.add_argument("--composition", type=Path,
                   help="long-form CSV: donor, celltype, n1, n2")
    p.add_argument("--expression", type=Path,
                   help="long-form CSV: donor, celltype, n_eff, floor")
    p.add_argument("--out", type=Path, default=Path("technical_floor_out"),
                   help="output directory (default technical_floor_out)")
    p.add_argument("--min-nuclei", type=int, default=50,
                   help="donor-level minimum nucleus criterion, which both libraries must meet (default 50)")
    p.add_argument("--min-n-eff", type=float, default=0.0,
                   help="minimum effective nucleus count to enter the scaling fit (default 0, no floor)")
    p.add_argument("--abundance", type=float, default=0.10,
                   help="cell-type abundance the design tables are computed at (default 0.10)")
    p.add_argument("--target-pp", type=float, default=1.0,
                   help="difference the design tables resolve, in percentage points (default 1.0)")
    p.add_argument("--bootstrap", type=int, default=DEFAULT_BOOTSTRAP,
                   help=f"bootstrap resamples (default {DEFAULT_BOOTSTRAP}, as used for the deposited results)")
    p.add_argument("--seed", type=int, default=DEFAULT_SEED,
                   help=f"random seed (default {DEFAULT_SEED})")
    p.add_argument("-v", "--verbose", action="store_true")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)s %(name)s: %(message)s",
    )
    if not args.composition and not args.expression:
        logger.error("give at least one of --composition or --expression")
        return 2

    config = FloorConfig(
        min_nuclei_per_pair=args.min_nuclei, min_n_eff=args.min_n_eff,
        n_bootstrap=args.bootstrap, seed=args.seed,
        abundance=args.abundance, target_pp=args.target_pp,
    )
    args.out.mkdir(parents=True, exist_ok=True)
    print(f"technical_floor · abundance {config.abundance:.0%} · target {config.target_pp} percentage points\n")

    comp = None
    if args.composition:
        # Two accepted input shapes. Raw counts can be filtered by the pairing
        # criterion here; precomputed z/n_eff cannot, because the criterion
        # needs the per-library nucleus counts, so the user is warned instead.
        df = read_composition(args.composition)
        if {"n1", "n2"}.issubset(df.columns):
            df = apply_pair_criterion(df, config.min_nuclei_per_pair)
            df = df[df.included_in_fit]
        else:
            logger.warning(
                "the composition input is the precomputed z/n_eff form and carries no "
                "n1/n2, so the >=%d nucleus criterion cannot be applied; check that the "
                "input is already filtered", config.min_nuclei_per_pair)
        comp = estimate_composition_floor(df, config)
        print(comp.summary())
        comp.per_celltype.to_csv(args.out / "composition_kappa_by_celltype.csv", index=False)
        composition_design_table(comp, config).to_csv(
            args.out / "composition_design_table.csv", index=False)
        req = composition_requirement(comp, config)
        req.to_csv(args.out / "composition_requirement.csv", index=False)
        print("\n  nuclei per group to resolve "
              f"{config.target_pp} percentage points:")
        for _, r in req.iterrows():
            print(f"    {r.basis:<32} κ={r.kappa:<6} {int(r.nuclei_per_group):>9,}")
        print()

    if args.expression:
        edf = read_expression(args.expression)
        # Same two shapes, same reason for the warning branch.
        if {"n1", "n2"}.issubset(edf.columns):
            # The pairing criterion governs the expression layer too. In the
            # paper's own analysis it was at one point applied to composition
            # only, which is why both branches call the same function here.
            edf = apply_pair_criterion(edf, config.min_nuclei_per_pair)
            edf = edf[edf.included_in_fit]
        else:
            logger.warning(
                "the expression input carries no n1/n2, so the >=%d nucleus pair "
                "criterion cannot be applied; the fit will use every donor",
                config.min_nuclei_per_pair)
        expr = estimate_expression_floor(edf, config)
        print(expr.summary())
        tbl = expression_design_table(expr)
        tbl.to_csv(args.out / "expression_design_table.csv", index=False)
        print("\n  |log2 fold change| technical noise thresholds:")
        for _, r in tbl.iterrows():
            flag = "  (extrapolated)" if r.extrapolated else ""
            print(f"    n_eff={int(r.n_eff):>5}   {r.log2fc_threshold:.2f}{flag}")
        print()

    print(f"written to {args.out}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
