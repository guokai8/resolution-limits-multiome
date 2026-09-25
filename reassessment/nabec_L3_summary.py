#!/usr/bin/env python3
"""汇总 NABEC/HBCC 的 Layer 3 扫描：OR 对 细胞数 n 与 ATAC 深度 的关系。

⚠️ 等化设计里 n 与深度是【反向耦合】的——n 越大，共同深度阈值越低。
   所以必须把两者分开看，不能只报 OR vs n。
"""
import json, glob, numpy as np, pandas as pd
rows=[json.load(open(f)) for f in sorted(glob.glob("results_methods/nabec_L3/*.json"))]
d=pd.DataFrame(rows).sort_values(["celltype","n"])
d["enriched"]=d.lo>1; d["depleted"]=d.hi<1
print(f"配置数 {len(d)}   细胞类型 {d.celltype.nunique()}\n")
print(d[["celltype","n","depth_atac","n_tested_all","n_links","OR","lo","hi"]]
      .to_string(index=False,float_format=lambda x:f"{x:.3f}"))
print(f"\n显著富集 (CI>1): {int(d.enriched.sum())}/{len(d)}   显著耗竭 (CI<1): {int(d.depleted.sum())}/{len(d)}")
print("\n=== OR 与 n / 深度 的关系（log-log 斜率）===")
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
