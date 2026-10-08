#!/usr/bin/env python3
"""Summarise the Layer 3 sweep in the external NABEC/HBCC cohort.

Reports the promoter-enrichment odds ratio against two things at once: the
nucleus count n, and the ATAC depth.

Both are needed because in an equalised design the two are INVERSELY COUPLED.
Equalising forces every cell type to the smallest one's nucleus count, and the
common depth threshold is then whatever all of those nuclei clear -- so raising
n lowers the depth every nucleus is held to. Reporting OR against n alone would
attribute to nucleus number an effect that is partly a loss of depth.
"""
import json, glob, numpy as np, pandas as pd
rows=[json.load(open(f)) for f in sorted(glob.glob("results_methods/nabec_L3/*.json"))]
d=pd.DataFrame(rows).sort_values(["celltype","n"])
# Significance is read off the interval rather than a p-value: the question is
# whether the enrichment is distinguishable from unity in either direction.
d["enriched"]=d.lo>1; d["depleted"]=d.hi<1
print(f"配置数 {len(d)}   细胞类型 {d.celltype.nunique()}\n")
print(d[["celltype","n","depth_atac","n_tested_all","n_links","OR","lo","hi"]]
      .to_string(index=False,float_format=lambda x:f"{x:.3f}"))
print(f"\n显著富集 (CI>1): {int(d.enriched.sum())}/{len(d)}   显著耗竭 (CI<1): {int(d.depleted.sum())}/{len(d)}")
print("\n=== OR 与 n / 深度 的关系（log-log 斜率）===")
# Fit log OR against each of the two coupled resources separately. Neither
# slope is a causal estimate; together they show the measure responds to both.
for v in ("n","depth_atac"):
    ok=(d[v]>0)&(d.OR>0)
    b=np.polyfit(np.log(d.loc[ok,v]),np.log(d.loc[ok,"OR"]),1)[0]
    r=np.corrcoef(np.log(d.loc[ok,v]),np.log(d.loc[ok,"OR"]))[0,1]
    print(f"  log OR ~ log {v:11s}: 斜率 {b:+.3f}  r={r:+.3f}")
print("\n=== 与原队列对照 ===")
print(f"  原队列 (ALS, 25 亚型, n=150, 深度 5,564): OR 0.663 [0.568, 0.768]")
print(f"  本数据集深度范围: {d.depth_atac.min():,.0f}–{d.depth_atac.max():,.0f}"
      f"  （原队列的 {d.depth_atac.min()/5564:.0f}–{d.depth_atac.max()/5564:.0f} 倍）")
d.to_csv("results_methods/nabec_L3/summary.csv",index=False)
