#!/usr/bin/env python3
"""N1 汇总：trans 零假设的稳定性与染色体敏感性。

输出三块：
  A  10 次独立 trans 配对下位置 OR 的分布（看它是否取决于某一次配对）
  B  留一染色体 jackknife（看它是否被个别染色体主导）
  C  高信息量档上的表现（观测应明显带位置信息，trans 不应当）

用法：
  python3 c1_trans_null_summary.py --grid n1_trans_validation.csv
"""
from __future__ import annotations

import argparse
from ast import literal_eval
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


def ha_or(a: float, b: float, c: float, d: float) -> float:
    """Haldane–Anscombe 校正的 OR。"""
    return ((a + .5) * (d + .5)) / ((b + .5) * (c + .5))


def arm_of(r: pd.Series) -> str:
    return "single-cell" if r.linker == "single_cell" else "ArchR 默认"


def parse_chrom(v) -> Dict[str, List[int]]:
    return literal_eval(v) if isinstance(v, str) else (v or {})


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", required=True)
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    df = pd.read_csv(a.grid)
    obs = df[df.null == "none"]
    tr = df[df.null == "trans"]

    rows = []
    print("A  10 次独立 trans 配对下的位置 OR")
    print(f"{'臂':<13s} {'n':>6s} {'链接':>7s} {'观测 OR':>8s} {'trans 中位':>10s} "
          f"{'trans 范围':>16s} {'位置 OR':>8s} {'位置范围':>16s} {'重复':>3s}")
    for (n, linker, archr), g in obs.groupby(["n", "linker", "archr_defaults"]):
        o = g.iloc[0]
        t = tr[(tr.n == n) & (tr.linker == linker) & (tr.archr_defaults == archr)]
        if t.empty:
            continue
        o_or = ha_or(o.prox_links, o.n_links - o.prox_links,
                     o.prox_tested_all, o.n_tested_all - o.prox_tested_all)
        t_ors = np.array([ha_or(r.prox_links, r.n_links - r.prox_links,
                                r.prox_tested_all, r.n_tested_all - r.prox_tested_all)
                          for _, r in t.iterrows()])
        pos = o_or / t_ors
        flag = " ⚠ 链接过少" if o.n_links < 200 else ""
        print(f"{arm_of(o):<13s} {int(n):>6d} {int(o.n_links):>7,} {o_or:>8.3f} "
              f"{np.median(t_ors):>10.3f} {t_ors.min():>7.3f}–{t_ors.max():<8.3f} "
              f"{np.median(pos):>8.3f} {pos.min():>7.3f}–{pos.max():<8.3f} "
              f"{len(t):>3d}{flag}")
        rows.append(dict(arm=arm_of(o), n=int(n), observed_OR=o_or,
                         trans_OR_median=float(np.median(t_ors)),
                         trans_OR_min=float(t_ors.min()), trans_OR_max=float(t_ors.max()),
                         trans_OR_cv=float(t_ors.std(ddof=1) / t_ors.mean()),
                         positional_OR_median=float(np.median(pos)),
                         positional_OR_min=float(pos.min()),
                         positional_OR_max=float(pos.max()), n_reps=len(t)))

    print("\nB  留一染色体 jackknife（位置 OR 用一次 trans 配对；列出偏离最大的三条）")
    print(f"{'臂':<13s} {'n':>6s} {'全量':>8s} {'jackknife 范围':>18s} {'影响最大的染色体':>24s}")
    for (n, linker, archr), g in obs.groupby(["n", "linker", "archr_defaults"]):
        o = g.iloc[0]
        t = tr[(tr.n == n) & (tr.linker == linker) & (tr.archr_defaults == archr)
               & (tr.trans_rep == 1)]
        if t.empty:
            continue
        oc = parse_chrom(o.per_chrom)
        tc = parse_chrom(t.iloc[0].per_chrom)
        if not oc or not tc:
            continue
        chroms = sorted(set(oc) | set(tc))
        def pos_drop(drop: str) -> float:
            oa = ob = ta = tb = 0
            for ch in chroms:
                if ch == drop:
                    continue
                x = oc.get(ch, [0, 0, 0, 0]); y = tc.get(ch, [0, 0, 0, 0])
                oa += x[0]; ob += x[1]; ta += y[0]; tb += y[1]
            if ob == 0 or tb == 0:
                return float("nan")
            return ha_or(oa, ob, 0, 0) / ha_or(ta, tb, 0, 0) if False else \
                   ((oa + .5) / (ob + .5)) / ((ta + .5) / (tb + .5))
        full = pos_drop("__none__")
        vals = {ch: pos_drop(ch) for ch in chroms}
        vals = {k: v for k, v in vals.items() if np.isfinite(v)}
        if not vals or not np.isfinite(full):
            print(f"{arm_of(o):<13s} {int(n):>6d} "
                  f"{'链接太少，jackknife 无意义':>44s}")
            continue
        worst = sorted(vals.items(), key=lambda kv: -abs(kv[1] - full))[:3]
        print(f"{arm_of(o):<13s} {int(n):>6d} {full:>8.3f} "
              f"{min(vals.values()):>8.3f}–{max(vals.values()):<8.3f} "
              f"{', '.join(f'{c}:{v:.2f}' for c, v in worst):>24s}")
        rows.append(dict(arm=arm_of(o), n=int(n), jackknife_full=full,
                         jackknife_min=min(vals.values()),
                         jackknife_max=max(vals.values()),
                         jackknife_worst=worst[0][0]))

    print("\n读法：A 中 trans OR 的范围窄、且位置 OR 的范围不跨过 1，说明结论不依赖")
    print("      单次配对；B 中 jackknife 范围窄，说明不被个别染色体主导；高核数档上")
    print("      观测 OR 远高于 trans OR，说明零假设在该有位置结构的地方确实识别得出。")
    if a.out:
        pd.DataFrame(rows).to_csv(a.out, index=False)
        print(f"\n写出 {a.out}")


if __name__ == "__main__":
    main()
