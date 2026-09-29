"""命令行入口：python3 -m technical_floor"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from .composition import estimate_composition_floor
from .config import FloorConfig
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
        description="从重复文库测量单核多组学的组成层与表达层分辨率下限。",
    )
    p.add_argument("--composition", type=Path,
                   help="长表 CSV：donor, celltype, n1, n2")
    p.add_argument("--expression", type=Path,
                   help="长表 CSV：donor, celltype, n_eff, floor")
    p.add_argument("--out", type=Path, default=Path("technical_floor_out"),
                   help="输出目录（默认 technical_floor_out）")
    p.add_argument("--min-nuclei", type=int, default=50,
                   help="供体级最小细胞数准则，两个文库都要达到（默认 50）")
    p.add_argument("--min-n-eff", type=float, default=0.0,
                   help="进入标度拟合的最小有效细胞数（默认 0，即不设限）")
    p.add_argument("--abundance", type=float, default=0.10,
                   help="设计表所针对的细胞类型丰度（默认 0.10）")
    p.add_argument("--target-pp", type=float, default=1.0,
                   help="设计表要分辨的差异，单位百分点（默认 1.0）")
    p.add_argument("--bootstrap", type=int, default=2000, help="自助次数（默认 2000）")
    p.add_argument("--seed", type=int, default=0, help="随机种子（默认 0）")
    p.add_argument("-v", "--verbose", action="store_true")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)s %(name)s: %(message)s",
    )
    if not args.composition and not args.expression:
        logger.error("至少要给 --composition 或 --expression 之一")
        return 2

    config = FloorConfig(
        min_nuclei_per_pair=args.min_nuclei, min_n_eff=args.min_n_eff,
        n_bootstrap=args.bootstrap, seed=args.seed,
        abundance=args.abundance, target_pp=args.target_pp,
    )
    args.out.mkdir(parents=True, exist_ok=True)
    print(f"technical_floor · 丰度 {config.abundance:.0%} · 目标 {config.target_pp} 个百分点\n")

    comp = None
    if args.composition:
        df = read_composition(args.composition)
        if {"n1", "n2"}.issubset(df.columns):
            df = apply_pair_criterion(df, config.min_nuclei_per_pair)
            df = df[df.included_in_fit]
        else:
            logger.warning(
                "组成层输入是预算好的 z/n_eff 形式，没有 n1/n2，无法执行 >=%d 核准则；"
                "请确认输入已经过滤", config.min_nuclei_per_pair)
        comp = estimate_composition_floor(df, config)
        print(comp.summary())
        comp.per_celltype.to_csv(args.out / "composition_kappa_by_celltype.csv", index=False)
        composition_design_table(comp, config).to_csv(
            args.out / "composition_design_table.csv", index=False)
        req = composition_requirement(comp, config)
        req.to_csv(args.out / "composition_requirement.csv", index=False)
        print("\n  每组所需细胞数，分辨 "
              f"{config.target_pp} 个百分点：")
        for _, r in req.iterrows():
            print(f"    {r.basis:<32} κ={r.kappa:<6} {int(r.nuclei_per_group):>9,}")
        print()

    if args.expression:
        edf = read_expression(args.expression)
        if {"n1", "n2"}.issubset(edf.columns):
            # 同一条配对准则必须同时管表达层；论文原来只把它落到了组成层上
            edf = apply_pair_criterion(edf, config.min_nuclei_per_pair)
            edf = edf[edf.included_in_fit]
        else:
            logger.warning(
                "表达层输入没有 n1/n2，无法执行 >=%d 核配对准则；"
                "拟合将使用全部供体", config.min_nuclei_per_pair)
        expr = estimate_expression_floor(edf, config)
        print(expr.summary())
        tbl = expression_design_table(expr)
        tbl.to_csv(args.out / "expression_design_table.csv", index=False)
        print("\n  |log2 fold change| 技术噪声阈值：")
        for _, r in tbl.iterrows():
            flag = "  (外推)" if r.extrapolated else ""
            print(f"    n_eff={int(r.n_eff):>5}   {r.log2fc_threshold:.2f}{flag}")
        print()

    print(f"写出到 {args.out}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
