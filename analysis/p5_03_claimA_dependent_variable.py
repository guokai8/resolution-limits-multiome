#!/usr/bin/env python3
"""
P5 步骤 0.6 · 命题 A 的因变量：逐细胞类型的 TDP-43 病理富集度
================================================================
新颖性重估后（decisions/DECISION_P5_02），命题 A 是唯一存活的主轴，
本脚本产出的就是它要预测的**因变量**。这一步不成立，整个项目不成立。

要回答的问题：**哪些神经元类型在核内 TDP-43 丧失的核里被富集？**

两个独立的统计量（必须都算，并比较）：

  OR_internal  = (n_Low_c / n_High_c) / (n_Low_rest / n_High_rest)
       在 FANS 内部，以同一组织、同一供体、同一次分选得到的 TDP-43 High 核
       作为基线。这是**供体配对的病例-对照设计**，不需要外部分母。
       ⚠️ 风险：High 分选本身可能有组成偏倚（分选阈值不是随机抽样）。

  OR_external = (n_Low_c / n_Low_rest) / (n_Multi_c / n_Multi_rest)
       以 Multiome 中**疾病供体的神经元组成**作为基线。
       ⚠️ 风险：FANS 是 NeuN+ 分选，Multiome 未分选，两者流程不同。
       ⚠️ 分母必须用我们自算的组成，**不能用 Supplementary Data 22**——
          该表的 n_Exc_THEMIS_L6_ITC 整列为 0（见 FINDING_P5_01）。

  两个 OR 一致 → 结论稳健；不一致 → 必须报告并解释，不能挑一个。

置信区间：对 FANS 的**供体池**做 bootstrap 重抽（不是对细胞），
因为细胞不是独立观测单位。Multiome 端对供体重抽。

只用元数据，不碰计数矩阵，秒级完成。

用法：
  python3 p5_03_claimA_dependent_variable.py --data DIR --out results/ [--level WNN_L2.5]
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

FANS_LEVEL = {"WNN_L2.5": "ID_WNN_L25_Predicted",
              "WNN_L3":   "ID_WNN_L3_Predicted",
              "WNN_L4":   "ID_WNN_L4_Predicted"}


def load_fans(data_dir, level, drop_inconsistent=True, multi=None):
    """
    ⚠️ 必须用 *_Predicted 列。FANS 元数据的 ID_WNN_L1..L4 六列是坏的
    （每行的值都是前六列的**列名**字符串）。见 docs/04_DATA_NOTES.md。
    """
    fa = pd.read_csv(os.path.join(data_dir, "FANS_Dataset_Metadata.txt"),
                     sep="\t", dtype=str)
    col = FANS_LEVEL[level]
    if col not in fa.columns:
        sys.exit(f"[FATAL] FANS 元数据缺少 {col}")
    bad = fa[col].dropna().unique()
    if len(bad) == 1 and str(bad[0]).startswith("ID_WNN"):
        sys.exit(f"[FATAL] {col} 是常量列，读错了列。用 *_Predicted。")

    n0 = len(fa)
    if drop_inconsistent and multi is not None and level != "WNN_L4":
        # 用 Multiome 建立 L4 -> level 的真映射，过滤 FANS 中自相矛盾的核
        m = (multi[["WNN_L4", level]].drop_duplicates().dropna()
             .set_index("WNN_L4")[level].to_dict())
        l4 = fa["ID_WNN_L4_Predicted"]
        known = l4.isin(m)
        implied = l4.map(m)
        drop = known & (implied != fa[col])
        fa = fa[~drop]
        print(f"  FANS: 过滤 L4 与 {level} 自相矛盾的核 {int(drop.sum())} "
              f"({100*drop.sum()/n0:.1f}%)，剩 {len(fa)}")
    fa = fa[fa["TDP43"].isin(["High", "Low"])]
    return fa[[col, "TDP43", "Sample_donor"]].rename(columns={col: "ct"})


def or_with_ci(tab_fn, groups, n_boot=2000, seed=0):
    """tab_fn(group_subset) -> DataFrame[ct, a, b]；对 groups 做 bootstrap。"""
    rng = np.random.default_rng(seed)
    point = tab_fn(groups)
    boots = []
    for _ in range(n_boot):
        g = rng.choice(groups, size=len(groups), replace=True)
        try:
            boots.append(tab_fn(list(g))["OR"])
        except Exception:
            continue
    B = pd.concat(boots, axis=1)
    point["lo"] = B.quantile(0.025, axis=1)
    point["hi"] = B.quantile(0.975, axis=1)
    point["boot_p"] = (B <= 1).mean(axis=1).clip(lower=1 / len(boots))
    return point


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", default="results")
    ap.add_argument("--level", default="WNN_L2.5",
                    choices=list(FANS_LEVEL))
    ap.add_argument("--min-cells", type=int, default=100,
                    help="FANS 中细胞数低于此值的类型不报（NeuN+ 分选使胶质不可用）")
    ap.add_argument("--n-boot", type=int, default=2000)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    print(f"[1/4] 读 Multiome 元数据 (level={args.level}) …")
    mm = pd.read_csv(os.path.join(args.data, "Multiome_Dataset_Metadata.txt"),
                     sep="\t", usecols=["ID", "Case", "WNN_L1", "WNN_L4",
                                        args.level], dtype=str)

    print("[2/4] 读并清洗 FANS …")
    fa = load_fans(args.data, args.level, multi=mm)
    print(f"  FANS: {len(fa)} 核，High={int((fa.TDP43=='High').sum())}, "
          f"Low={int((fa.TDP43=='Low').sum())}，"
          f"{fa.Sample_donor.nunique()} 个供体池")

    # ---- 基线：Multiome 的疾病供体 + 仅神经元（对齐 NeuN+ 分选）------------
    dis = mm[mm["Case"].isin(["ALS", "ALS_FTD"]) & (mm["WNN_L1"] == "Neuronal")]
    print(f"  Multiome 基线: {len(dis)} 个疾病供体的神经元核，"
          f"{dis.ID.nunique()} 个供体")

    # 只保留两边都有、且 FANS 中细胞数够的神经元类型
    fa_n = fa[fa["ct"].isin(set(dis[args.level]))]
    keep = [c for c, n in fa_n["ct"].value_counts().items()
            if n >= args.min_cells]
    fa_n = fa_n[fa_n["ct"].isin(keep)]
    dis = dis[dis[args.level].isin(keep)]
    print(f"  可用神经元类型: {len(keep)} 个 -> {sorted(keep)}")
    if len(keep) < 4:
        print("  ⚠️ 类型数过少，命题 A 的回归几乎没有功效。")

    # ---- 统计量 1：FANS 内部 Low vs High --------------------------------
    def internal(pools):
        d = fa_n[fa_n.Sample_donor.isin(pools)]
        t = pd.crosstab(d["ct"], d["TDP43"]).reindex(keep).fillna(0)
        for c in ("High", "Low"):
            if c not in t:
                t[c] = 0
        a, b = t["Low"].astype(float), t["High"].astype(float)
        # Haldane-Anscombe 校正，避免 0 导致 OR 发散
        orr = ((a + .5) / (b + .5)) / (((a.sum() - a) + .5) /
                                       ((b.sum() - b) + .5))
        return pd.DataFrame({"n_High": t["High"], "n_Low": t["Low"], "OR": orr})

    print("[3/4] 统计量 1：FANS 内部 (Low vs High)，供体池 bootstrap …")
    r1 = or_with_ci(internal, sorted(fa_n.Sample_donor.unique()),
                    n_boot=args.n_boot)

    # ---- 统计量 2：FANS-Low vs Multiome 疾病神经元组成 -------------------
    low = fa_n[fa_n.TDP43 == "Low"]
    base = dis[args.level].value_counts().reindex(keep).fillna(0).astype(float)

    def external(donors):
        b = dis[dis.ID.isin(donors)][args.level].value_counts() \
               .reindex(keep).fillna(0).astype(float)
        a = low["ct"].value_counts().reindex(keep).fillna(0).astype(float)
        orr = ((a + .5) / (b + .5)) / (((a.sum() - a) + .5) /
                                       ((b.sum() - b) + .5))
        return pd.DataFrame({"OR": orr})

    print("[4/4] 统计量 2：FANS-Low vs Multiome 基线，供体 bootstrap …")
    r2 = or_with_ci(external, sorted(dis.ID.unique()), n_boot=args.n_boot)

    out = pd.DataFrame({
        "celltype": keep,
        "n_High": r1["n_High"].reindex(keep).values,
        "n_Low": r1["n_Low"].reindex(keep).values,
        "n_Multiome_baseline": base.reindex(keep).values,
        "OR_internal": r1["OR"].reindex(keep).values,
        "OR_int_lo": r1["lo"].reindex(keep).values,
        "OR_int_hi": r1["hi"].reindex(keep).values,
        "OR_external": r2["OR"].reindex(keep).values,
        "OR_ext_lo": r2["lo"].reindex(keep).values,
        "OR_ext_hi": r2["hi"].reindex(keep).values,
    }).sort_values("OR_internal", ascending=False)
    out["rank_internal"] = out["OR_internal"].rank(ascending=False).astype(int)
    out["rank_external"] = out["OR_external"].rank(ascending=False).astype(int)

    f = os.path.join(args.out, f"p5_claimA_pathology_OR_{args.level.replace('.','')}.csv")
    out.to_csv(f, index=False)

    print("\n" + "=" * 92)
    print(f"命题 A 的因变量 · {args.level}")
    print("=" * 92)
    print(f"{'细胞类型':<20}{'High':>7}{'Low':>7}{'OR_int':>9}"
          f"{'95%CI':>18}{'OR_ext':>9}{'一致?':>8}")
    for _, r in out.iterrows():
        agree = "✅" if (r.OR_internal - 1) * (r.OR_external - 1) > 0 else "⚠️不同向"
        sig = "*" if (r.OR_int_lo > 1 or r.OR_int_hi < 1) else " "
        print(f"{r.celltype:<20}{int(r.n_High):>7}{int(r.n_Low):>7}"
              f"{r.OR_internal:>8.2f}{sig}"
              f"{f'[{r.OR_int_lo:.2f},{r.OR_int_hi:.2f}]':>18}"
              f"{r.OR_external:>9.2f}{agree:>8}")

    rho = out[["rank_internal", "rank_external"]].corr(method="spearman").iloc[0, 1]
    print(f"\n两个统计量的排序一致性 Spearman ρ = {rho:.3f}")
    if rho < 0.5:
        print("  ⚠️ 两个基线给出的排序差异大。**不能挑一个报**，必须并列并讨论。")
    else:
        print("  ✅ 两个独立基线给出一致排序 —— 因变量稳健。")
    print(f"\n[done] {f}")
    print("\n下一步：把这个排序作为命题 A 的预测目标（A3 的留一交叉验证）。")


if __name__ == "__main__":
    main()
