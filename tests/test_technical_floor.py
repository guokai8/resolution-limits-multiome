"""technical_floor 的测试。

最重要的一组是 TestReproducesPublished：工具跑在论文自己沉积的表上，必须
还原论文报告的每一个数。工具与论文之间任何一处不合，要么是工具错了，要么
是论文错了——两种都必须先查清再往下走。
"""

from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
# 包在本仓库位于根目录，在论文工程里位于 code/；两种布局都要能跑
for candidate in (ROOT, ROOT / "code"):
    if (candidate / "technical_floor").is_dir():
        sys.path.insert(0, str(candidate))
        break

from technical_floor import (  # noqa: E402
    FloorConfig,
    apply_pair_criterion,
    composition_requirement,
    estimate_composition_floor,
    estimate_expression_floor,
    expression_design_table,
    floor_from_counts,
    nuclei_required,
    read_composition,
    read_expression,
    resolvable_difference,
)

EXAMPLES = ROOT / "data" / "example_inputs"


class TestFloorFromCounts(unittest.TestCase):
    def test_identical_libraries_give_zero_floor(self) -> None:
        counts = np.arange(1, 501, dtype=float)
        self.assertAlmostEqual(floor_from_counts(counts, counts, min_genes=200), 0.0, places=12)

    def test_returns_none_below_min_genes(self) -> None:
        self.assertIsNone(floor_from_counts([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], min_genes=200))

    def test_returns_none_on_empty_library(self) -> None:
        z = np.zeros(300)
        self.assertIsNone(floor_from_counts(np.ones(300), z, min_genes=200))

    def test_length_mismatch_raises(self) -> None:
        with self.assertRaises(ValueError):
            floor_from_counts([1.0, 2.0], [1.0, 2.0, 3.0])

    def test_union_rule_includes_genes_seen_in_one_member(self) -> None:
        a = np.concatenate([np.ones(250), np.zeros(50)])
        b = np.concatenate([np.ones(250), np.ones(50)])
        self.assertIsNotNone(floor_from_counts(a, b, min_genes=300))


class TestDesignArithmetic(unittest.TestCase):
    def test_cost_scales_as_kappa_squared(self) -> None:
        one = nuclei_required(1.0)
        two = nuclei_required(2.0)
        self.assertAlmostEqual(two / one, 4.0, places=6)

    def test_multinomial_bound_matches_published(self) -> None:
        self.assertEqual(nuclei_required(1.0), 6915)

    def test_resolvable_difference_falls_as_sqrt_n(self) -> None:
        a = resolvable_difference(1.0, 20, 5000)
        b = resolvable_difference(1.0, 20, 20000)
        self.assertAlmostEqual(a / b, 2.0, places=6)

    def test_rejects_nonpositive_design(self) -> None:
        with self.assertRaises(ValueError):
            resolvable_difference(1.0, 0, 5000)


class TestConfigValidation(unittest.TestCase):
    def test_rejects_impossible_abundance(self) -> None:
        for bad in (0.0, 1.0, -0.2, 1.5):
            with self.assertRaises(ValueError):
                FloorConfig(abundance=bad)

    def test_rejects_nonpositive_target(self) -> None:
        with self.assertRaises(ValueError):
            FloorConfig(target_pp=0.0)

    def test_is_frozen(self) -> None:
        cfg = FloorConfig()
        with self.assertRaises(Exception):
            cfg.abundance = 0.2  # type: ignore[misc]


class TestPairCriterion(unittest.TestCase):
    def test_excludes_donor_failing_either_library(self) -> None:
        df = pd.DataFrame({
            "donor": ["A", "A", "B", "B"],
            "celltype": ["x", "y", "x", "y"],
            "n1": [100, 100, 100, 100],
            "n2": [100, 100, 10, 10],
        })
        out = apply_pair_criterion(df, min_nuclei=50)
        self.assertTrue(out[out.donor == "A"].included_in_fit.all())
        self.assertFalse(out[out.donor == "B"].included_in_fit.any())

    def test_criterion_is_per_donor_not_per_celltype(self) -> None:
        # 单个细胞类型很小不该踢掉供体；供体汇总达标即可
        df = pd.DataFrame({
            "donor": ["A", "A"], "celltype": ["big", "tiny"],
            "n1": [95, 5], "n2": [95, 5],
        })
        self.assertTrue(apply_pair_criterion(df, 50).included_in_fit.all())


@unittest.skipUnless(EXAMPLES.exists(), "example inputs not built")
class TestReproducesPublished(unittest.TestCase):
    """工具跑论文的沉积数据，必须还原论文的数。"""

    @classmethod
    def setUpClass(cls) -> None:
        cls.cfg = FloorConfig(n_bootstrap=400, seed=0)
        cls.comp = estimate_composition_floor(
            read_composition(EXAMPLES / "motor_cortex_composition_pairs.csv"), cls.cfg)
        expr = read_expression(EXAMPLES / "motor_cortex_expression_floors.csv")
        expr = apply_pair_criterion(expr, cls.cfg.min_nuclei_per_pair)
        cls.expr = estimate_expression_floor(expr[expr.included_in_fit], cls.cfg)

    def test_kappa(self) -> None:
        self.assertAlmostEqual(self.comp.kappa, 4.27, places=2)

    def test_pair_and_donor_counts(self) -> None:
        self.assertEqual(self.comp.n_pairs, 300)
        self.assertEqual(self.comp.n_donors, 26)

    def test_per_celltype_kappa_range(self) -> None:
        self.assertAlmostEqual(self.comp.per_celltype.kappa.min(), 1.53, places=2)
        self.assertAlmostEqual(self.comp.per_celltype.kappa.max(), 8.50, places=2)

    def test_kappa_rises_with_abundance(self) -> None:
        self.assertGreater(self.comp.abundance_r, 0.8)

    def test_design_requirements(self) -> None:
        req = composition_requirement(self.comp, self.cfg).set_index("basis").nuclei_per_group
        self.assertEqual(int(req["multinomial bound"]), 6915)
        self.assertEqual(int(req["pooled kappa"]), 125949)
        self.assertEqual(int(req["median cell type"]), 44435)
        self.assertEqual(int(req["most overdispersed cell type"]), 499575)

    def test_expression_fit(self) -> None:
        self.assertAlmostEqual(self.expr.exponent, -0.497, places=3)
        self.assertAlmostEqual(self.expr.coefficient, 4.58, places=2)
        self.assertEqual(self.expr.n_obs, 300)
        self.assertEqual(self.expr.n_donors, 26)

    def test_expression_thresholds(self) -> None:
        tbl = expression_design_table(self.expr).set_index("n_eff").log2fc_threshold
        for n, expected in ((25, 0.93), (50, 0.66), (100, 0.47), (200, 0.33), (1000, 0.15)):
            self.assertAlmostEqual(tbl[n], expected, places=2, msg=f"n_eff={n}")

    def test_the_pair_criterion_changes_the_answer(self) -> None:
        """不执行准则会得到论文最初报告的、已被撤回的那组数。"""
        raw = read_expression(EXAMPLES / "motor_cortex_expression_floors.csv")
        unfiltered = estimate_expression_floor(raw, self.cfg)
        self.assertEqual(unfiltered.n_donors, 27)
        self.assertAlmostEqual(unfiltered.exponent, -0.507, places=3)
        self.assertNotAlmostEqual(unfiltered.exponent, self.expr.exponent, places=3)


if __name__ == "__main__":
    unittest.main(verbosity=2)
