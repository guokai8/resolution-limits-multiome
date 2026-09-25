#!/usr/bin/env python3
"""
P5 Fig 1 · 混杂审计
====================
在做任何差异分析之前，把 79 个供体的 年龄 × 性别 × 队列 × 分组 结构摊开。

论文用 SVA 吸收了 PMI / 起病部位等隐藏变量。SVA 同时会吸收与分组共线的年龄效应，
使得"疾病效应"与"年龄效应"不可分。本脚本量化这个共线性，并输出一个
年龄匹配子集，供后续所有分析做敏感性复现。

输入（只读，全部已在盘）：
  Multiome_Dataset_Samples_Summary.txt   79 行，Case/Cohort/Sex/Age
  Multiome_Dataset_Metadata.txt          180,016 行逐核元数据（供体 ID 在 col 13）

输出：
  out/p5_fig1_donor_table.csv            供体级汇总表
  out/p5_fig1_agematched_donors.csv      年龄匹配子集的供体 ID 列表
  out/p5_fig1_confounding_stats.txt      文字报告

用法：
  python3 p5_00_confounding_audit.py --data /path/to/ResearchD --out out/
"""

import argparse
import os
import sys
from collections import Counter, defaultdict

import numpy as np
import pandas as pd


# --------------------------------------------------------------------------
# 供体表：Samples_Summary 没有供体 ID，只有 Case/Cohort/Sex/Age，
# 顺序与 Metadata 中的 ID 出现顺序不保证一致。因此以 Metadata 为准重建供体表，
# 再用 Samples_Summary 做行数与分布的一致性校验。
# --------------------------------------------------------------------------
def build_donor_table(data_dir):
    meta_path = os.path.join(data_dir, "Multiome_Dataset_Metadata.txt")
    usecols = ["ID", "Sex", "Case", "C9ORF72", "Case_Type",
               "WNN_L2", "WNN_L2.5", "WNN_L4", "Hist_Layer_CellType",
               "nCount_RNA", "nCount_ATAC"]
    meta = pd.read_csv(meta_path, sep="\t", usecols=usecols,
                       dtype={"ID": str, "Sex": str, "Case": str,
                              "Case_Type": str, "C9ORF72": "Int8"})
    donors = (meta.groupby("ID")
                  .agg(Sex=("Sex", "first"),
                       Case=("Case", "first"),
                       Case_Type=("Case_Type", "first"),
                       C9=("C9ORF72", "first"),
                       n_nuclei=("Case", "size"),
                       median_nCount_RNA=("nCount_RNA", "median"),
                       median_nCount_ATAC=("nCount_ATAC", "median"))
                  .reset_index())

    # Age / Cohort 只在 Samples_Summary 里。按 (Case, Sex) 分层做顺序配对，
    # 这是唯一可用的连接方式 —— 若作者提供了带 ID 的临床表（Supp Data 1/2），
    # 应改用那份表，本函数会在检测到 clinical.csv 时优先使用。
    clin = os.path.join(data_dir, "p5_clinical_with_ID.csv")
    if os.path.exists(clin):
        c = pd.read_csv(clin)
        donors = donors.merge(c[["ID", "Age", "Cohort"]], on="ID", how="left")
        donors.attrs["age_source"] = "clinical_table"
        return donors, meta

    summ = pd.read_csv(os.path.join(data_dir,
                                    "Multiome_Dataset_Samples_Summary.txt"),
                       sep="\t")
    summ["Case"] = summ["Case"].replace({"Control": "HC"})
    summ["Age"] = pd.to_numeric(summ["Age"], errors="coerce")

    # 分层顺序配对（近似；仅用于分布层面的审计，不用于逐供体建模）
    donors = donors.sort_values(["Case", "Sex", "ID"]).reset_index(drop=True)
    summ = summ.sort_values(["Case", "Sex"]).reset_index(drop=True)
    matched = []
    for (case, sex), grp in donors.groupby(["Case", "Sex"], sort=False):
        pool = summ[(summ.Case == case) & (summ.Sex == sex)]
        for i, (_, row) in enumerate(grp.iterrows()):
            if i < len(pool):
                matched.append((row.ID, pool.iloc[i].Age, pool.iloc[i].Cohort))
            else:
                matched.append((row.ID, np.nan, "UNKNOWN"))
    m = pd.DataFrame(matched, columns=["ID", "Age", "Cohort"])
    donors = donors.merge(m, on="ID", how="left")
    donors.attrs["age_source"] = "stratified_orderwise_match_APPROXIMATE"
    return donors, meta


def confounding_report(donors, fh):
    def p(*a):
        print(*a); print(*a, file=fh)

    p("=" * 74)
    p("P5 Fig 1 · 混杂审计")
    p("=" * 74)
    p(f"\n年龄来源: {donors.attrs.get('age_source')}")
    if "APPROXIMATE" in str(donors.attrs.get("age_source")):
        p("  ⚠️  当前用分层顺序配对近似。仅可用于分布层面的结论。")
        p("      逐供体建模前必须拿到带 ID 的临床表（Supplementary Data 1/2），")
        p("      放为 p5_clinical_with_ID.csv 后重跑本脚本。")

    p("\n--- 分组 × 年龄 ---")
    for case, g in donors.groupby("Case"):
        a = g.Age.dropna()
        p(f"  {case:10s} n={len(g):3d}  mean={a.mean():5.1f}  median={a.median():5.1f}"
          f"  IQR=[{a.quantile(.25):.0f},{a.quantile(.75):.0f}]  range={a.min():.0f}-{a.max():.0f}")

    # 年龄 ~ 分组 的方差解释比例
    grand = donors.Age.dropna().mean()
    ss_tot = ((donors.Age.dropna() - grand) ** 2).sum()
    ss_bet = sum(len(g.Age.dropna()) * (g.Age.dropna().mean() - grand) ** 2
                 for _, g in donors.groupby("Case"))
    p(f"\n  年龄方差被【分组】解释的比例  eta^2 = {ss_bet/ss_tot:.3f}")

    ss_bet_c = sum(len(g.Age.dropna()) * (g.Age.dropna().mean() - grand) ** 2
                   for _, g in donors.groupby("Cohort"))
    p(f"  年龄方差被【队列】解释的比例  eta^2 = {ss_bet_c/ss_tot:.3f}")
    p("  → 两者都高 = 年龄、队列、分组三者共线。任何只调整其一的模型都不够。")

    p("\n--- 队列 × 分组 列联表 ---")
    ct = pd.crosstab(donors.Cohort, donors.Case)
    p(ct.to_string())
    zero_ctrl = [c for c in ct.index if ct.loc[c].get("HC", 0) == 0]
    if zero_ctrl:
        p(f"\n  ⚠️  无对照的队列: {zero_ctrl}")
        p("      这些队列的样本对 疾病-vs-对照 的贡献完全依赖跨队列比较，")
        p("      任何队列效应都会直接伪装成疾病效应。")
        p("      对策：报告 (a) 全样本 (b) 剔除无对照队列 两个版本。")

    p("\n--- 对照组年龄的双峰结构（最危险的一处） ---")
    for coh, g in donors[donors.Case == "HC"].groupby("Cohort"):
        p(f"  {coh:10s} n={len(g):2d}  {sorted(g.Age.dropna().astype(int).tolist())}")
    p("  → 若对照的年龄分布是双峰的，'平均年龄' 是无意义的匹配指标。")
    p("    必须做分布匹配（如倾向性得分或精确分箱），不能只匹配均值。")
    return


def age_matched_subset(donors, caliper=6.0, seed=0):
    """
    最近邻年龄匹配（无放回，caliper 内），HC 匹配到 ALS。
    返回匹配上的供体 ID。这是所有主结论的敏感性复现子集。
    """
    rng = np.random.default_rng(seed)
    cases = donors[(donors.Case == "ALS") & donors.Age.notna()].copy()
    ctrls = donors[(donors.Case == "HC") & donors.Age.notna()].copy()
    used, pairs = set(), []
    for _, cs in cases.sample(frac=1, random_state=seed).iterrows():
        cand = ctrls[~ctrls.ID.isin(used)]
        if cand.empty:
            break
        d = (cand.Age - cs.Age).abs()
        j = d.idxmin()
        if d.loc[j] <= caliper:
            used.add(cand.loc[j, "ID"])
            pairs.append((cs.ID, cand.loc[j, "ID"], cs.Age, cand.loc[j, "Age"]))
    return pd.DataFrame(pairs, columns=["case_ID", "ctrl_ID", "case_Age", "ctrl_Age"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="ResearchD 目录")
    ap.add_argument("--out", default="out")
    ap.add_argument("--caliper", type=float, default=6.0,
                    help="年龄匹配容差（岁）")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    donors, meta = build_donor_table(args.data)
    donors.to_csv(os.path.join(args.out, "p5_fig1_donor_table.csv"), index=False)

    with open(os.path.join(args.out, "p5_fig1_confounding_stats.txt"), "w") as fh:
        confounding_report(donors, fh)

        pairs = age_matched_subset(donors, caliper=args.caliper)
        msg = (f"\n--- 年龄匹配子集 (caliper={args.caliper} 岁) ---\n"
               f"  匹配上 {len(pairs)} 对 ALS-HC\n"
               f"  病例年龄 mean={pairs.case_Age.mean():.1f}  "
               f"对照年龄 mean={pairs.ctrl_Age.mean():.1f}  "
               f"绝对差 mean={(pairs.case_Age-pairs.ctrl_Age).abs().mean():.2f} 岁\n"
               f"  → 全部主结论必须在这个子集里存活，否则改写为'部分由年龄介导'。")
        print(msg); print(msg, file=fh)

        # 细胞类型 × 分组 的核数结构（功效均衡的输入）
        ct = pd.crosstab(meta["WNN_L2.5"], meta["Case"])
        s = "\n--- 细胞类型 × 分组 核数（DEG 数比较前必须做的功效均衡输入）---\n"
        print(s + ct.to_string()); print(s + ct.to_string(), file=fh)

    pairs.to_csv(os.path.join(args.out, "p5_fig1_agematched_donors.csv"), index=False)
    print(f"\n[done] 输出写入 {args.out}/")


if __name__ == "__main__":
    main()
