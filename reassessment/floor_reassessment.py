#!/usr/bin/env python3
"""把已发表的细胞类型组成差异放到技术下限上重新评估。

预注册：plan/2026-09-22_floor_reassessment_prespecification.md
（SHA-256 b5efb71d878b11c1c531f7408d6d9430c1d545cb0554649fdca41dfd03c75e94）

⚠️ 手稿自身的结论是过度离散不可跨数据集迁移，故**主分析用多项式下限 κ=1**
（数学下界，任何真实数据集都达不到），脑队列实测的 3.90/4.28 仅作标注过的敏感性。

每个 contrast = (数据集 × 注释层级 × 细胞类型 × 一对分组)：
  Δ            = 两组 donor 级比例均值之差的绝对值
  Var_κ(Δ)     = κ² p(1−p) [ Σ_A(1/n_d)/D_A² + Σ_B(1/n_d)/D_B² ]
  下限 F_κ     = 1.96 √Var_κ(Δ)
  常规检验     = donor 级比例的 Welch t 检验
主要结局：在常规检验 P<0.05 的 contrast 中，Δ ≤ F_κ 的比例。
"""

from __future__ import annotations

import argparse, itertools, json, logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats

logger = logging.getLogger("floor")
MIN_CELLS_PER_DONOR = 50
MIN_DONORS_PER_GROUP = 5
KAPPAS = {"multinomial": 1.0, "brain_seaad": 3.90, "brain_als": 4.28}


# ---------------------------------------------------------------- 数据加载
def _h5ad_obs(path: Path, cols: List[str]) -> pd.DataFrame:
    """只读 obs 的指定列（含 categorical 解码），不载入矩阵。"""
    import h5py
    out = {}
    with h5py.File(path, "r") as f:
        o = f["obs"]
        for c in cols:
            if c not in o:
                continue
            g = o[c]
            if isinstance(g, h5py.Group):           # categorical
                cats = np.asarray(g["categories"][:]).astype(str)
                out[c] = cats[g["codes"][:]]
            else:
                v = np.asarray(g[:])
                out[c] = v.astype(str) if v.dtype.kind in "SOU" else v
    return pd.DataFrame(out)


def load_D1(root: Path) -> Optional[Dict]:
    p = root / "Data/Other_Datasets/Multiome_Dataset/files/Multiome_Dataset_Metadata.txt"
    if not p.exists():
        return None
    m = pd.read_csv(p, sep="\t", quoting=3, dtype=str,
                    usecols=["ID", "Case", "WNN_L1", "WNN_L2.5"])
    return dict(name="D1_ALS_motor_cortex", tissue="brain", donor="ID", group="Case",
                levels={"coarse": "WNN_L1", "fine": "WNN_L2.5"}, obs=m,
                contrasts=[("ALS", "HC"), ("ALS_FTD", "HC")])


def load_D2(root: Path) -> Optional[Dict]:
    p = root / ("Data/Other_Datasets/SEAAD_MTG/files/"
                "SEAAD_MTG_RNAseq_final-nuclei_metadata.2026-06-22.csv")
    if not p.exists():
        return None
    cand = pd.read_csv(p, nrows=0).columns.tolist()
    grp = next((c for c in ("Overall AD neuropathological Change", "Cognitive Status",
                            "ADNC", "Braak stage") if c in cand), None)
    if grp is None:
        return None
    m = pd.read_csv(p, usecols=["Donor ID", "Class", "Subclass", grp],
                    dtype=str, low_memory=False).rename(columns={"Donor ID": "donor"})
    vals = [v for v in m[grp].dropna().unique()]
    # 预注册 §9：序数变量取两端
    order = ["Not AD", "Low", "Intermediate", "High"]
    present = [v for v in order if v in vals]
    con = [(present[-1], present[0])] if len(present) >= 2 else []
    return dict(name="D2_SEAAD_MTG", tissue="brain", donor="donor", group=grp,
                levels={"coarse": "Class", "fine": "Subclass"}, obs=m, contrasts=con)


def load_D3(root: Path) -> Optional[Dict]:
    p = root / "GSE290359_BA44_BA46_COUNTS.h5ad"
    if not p.exists():
        return None
    m = _h5ad_obs(p, ["Subject ID", "Edinburgh ID", "Disease?", "Classifcation",
                      "Cell_Type", "Cell_Subtype"])
    donor = "Subject ID" if "Subject ID" in m else "Edinburgh ID"
    if donor not in m or "Disease?" not in m:
        return None
    vals = sorted(m["Disease?"].dropna().unique())
    ctrl = next((v for v in vals if v.lower() in ("control", "no", "ctrl", "hc", "false")), None)
    con = [(v, ctrl) for v in vals if ctrl and v != ctrl]
    return dict(name="D3_GSE290359_BA44_BA46", tissue="brain", donor=donor, group="Disease?",
                levels={"coarse": "Cell_Type", "fine": "Cell_Subtype"}, obs=m, contrasts=con)


def load_D4(root: Path) -> Optional[Dict]:
    out = []
    for tag, fn in (("MTC", "GSE330130_MTC_SNRNASEQ_COUNTS.h5ad"),
                    ("LSC", "GSE330130_LSC_SNRNASEQ_COUNTS.h5ad")):
        p = root / fn
        if not p.exists():
            continue
        m = _h5ad_obs(p, ["Subject ID", "Diagnosis", "Cell_Type", "Cell_Subtype"])
        if "Subject ID" not in m or "Diagnosis" not in m:
            continue
        vals = sorted(m["Diagnosis"].dropna().unique())
        # 子串匹配：本数据集的对照写作 "Non-neurological control"
        ctrl = next((v for v in vals if "control" in v.lower() or v.lower() in ("ctrl", "hc")), None)
        con = [(v, ctrl) for v in vals if ctrl and v != ctrl]
        out.append(dict(name=f"D4_GSE330130_{tag}", tissue="brain", donor="Subject ID",
                        group="Diagnosis",
                        levels={"coarse": "Cell_Type", "fine": "Cell_Subtype"},
                        obs=m, contrasts=con))
    return out or None


def load_D5(root: Path) -> Optional[Dict]:
    cell = root / "Data/COVID_feasibility/raw_metadata/GSE158055_cell_annotation.csv.gz"
    meta = root / "Data/COVID_feasibility/raw_metadata/GSE158055_sample_metadata.xlsx"
    if not (cell.exists() and meta.exists()):
        return None
    c = pd.read_csv(cell, usecols=["sampleID", "celltype", "majorType"], dtype=str)
    s = pd.read_excel(meta, skiprows=20)
    s.columns = [str(x).strip() for x in s.columns]
    col_sample = next((x for x in s.columns if x.lower().startswith("sample name")), None)
    col_pat = next((x for x in s.columns if "Patients" in x), None)
    col_sev = next((x for x in s.columns if "severity" in x.lower()), None)
    col_type = next((x for x in s.columns if "Sample type" in x), None)
    if not all([col_sample, col_pat, col_sev]):
        return None
    s = s[[col_sample, col_pat, col_sev] + ([col_type] if col_type else [])].dropna(
        subset=[col_sample])
    s.columns = ["sampleID", "donor", "severity"] + (["sample_type"] if col_type else [])
    if "sample_type" in s:                       # 只保留 PBMC（血），避免混入其他组织
        s = s[s.sample_type.astype(str).str.contains("PBMC", case=False, na=False)]
    m = c.merge(s, on="sampleID", how="inner")
    vals = [v for v in m.severity.dropna().unique()]
    order = ["control", "mild/moderate", "severe/critical"]
    low = next((v for v in vals if v.lower().startswith("control")), None)
    hi = next((v for v in vals if "severe" in v.lower()), None)
    con = [(hi, low)] if (low and hi) else []
    return dict(name="D5_GSE158055_COVID_PBMC", tissue="blood", donor="donor",
                group="severity", levels={"coarse": "majorType", "fine": "celltype"},
                obs=m, contrasts=con)


# ---------------------------------------------------------------- 核心计算
def analyse(ds: Dict) -> pd.DataFrame:
    rows = []
    obs = ds["obs"]
    dcol, gcol = ds["donor"], ds["group"]
    obs = obs.dropna(subset=[dcol, gcol])
    # donor → 组（取该 donor 的唯一值；不唯一则丢弃并记录）
    dg = obs.groupby(dcol)[gcol].nunique()
    bad = dg[dg > 1].index.tolist()
    if bad:
        logger.warning("%s: %d 个 donor 分组不唯一，已排除", ds["name"], len(bad))
        obs = obs[~obs[dcol].isin(bad)]
    donor_group = obs.drop_duplicates(dcol).set_index(dcol)[gcol]
    ncell = obs.groupby(dcol).size()
    keep_donor = ncell[ncell >= MIN_CELLS_PER_DONOR].index
    obs = obs[obs[dcol].isin(keep_donor)]
    ncell = ncell.loc[keep_donor]

    for lvl_name, lvl_col in ds["levels"].items():
        if lvl_col not in obs.columns:
            continue
        sub = obs.dropna(subset=[lvl_col])
        tab = sub.groupby([dcol, lvl_col]).size().unstack(fill_value=0)
        tab = tab.reindex(index=ncell.index, fill_value=0)
        frac = tab.div(ncell, axis=0)
        pooled = tab.sum(0) / tab.sum().sum()
        for ga, gb in ds["contrasts"]:
            A = [d for d in frac.index if donor_group.get(d) == ga]
            B = [d for d in frac.index if donor_group.get(d) == gb]
            if len(A) < MIN_DONORS_PER_GROUP or len(B) < MIN_DONORS_PER_GROUP:
                continue
            invA = float(np.sum(1.0 / ncell.loc[A])) / len(A) ** 2
            invB = float(np.sum(1.0 / ncell.loc[B])) / len(B) ** 2
            for ct in frac.columns:
                fa, fb = frac.loc[A, ct].to_numpy(), frac.loc[B, ct].to_numpy()
                p = float(pooled[ct])
                if p <= 0:
                    continue
                delta = abs(fa.mean() - fb.mean())
                t, pv = stats.ttest_ind(fa, fb, equal_var=False)
                r = dict(dataset=ds["name"], tissue=ds["tissue"], level=lvl_name,
                         cell_type=str(ct), group_a=ga, group_b=gb,
                         n_donor_a=len(A), n_donor_b=len(B),
                         median_cells=float(ncell.loc[A + B].median()),
                         pooled_p=p, delta=delta, t=float(t), p_value=float(pv))
                for k, kv in KAPPAS.items():
                    var = kv ** 2 * p * (1 - p) * (invA + invB)
                    F = 1.96 * np.sqrt(var)
                    r[f"floor_{k}"] = float(F)
                    r[f"above_{k}"] = bool(delta > F)
                r["delta_over_floor_multinomial"] = delta / r["floor_multinomial"] \
                    if r["floor_multinomial"] > 0 else np.nan
                # ---- 事后补充（非预注册）：把细胞当独立个体的合并检验 ----
                # 文献中常见的做法：把组内各 donor 的细胞合并，对
                # (该类型 vs 其余) × (组A vs 组B) 做 2x2 检验。该检验假设细胞独立，
                # 正是手稿所否定的假设。
                ca, cb = int(tab.loc[A, ct].sum()), int(tab.loc[B, ct].sum())
                ta, tb = int(ncell.loc[A].sum()), int(ncell.loc[B].sum())
                table = np.array([[ca, ta - ca], [cb, tb - cb]])
                if table.min() >= 0 and table.sum() > 0:
                    try:
                        chi2, pnaive, _, _ = stats.chi2_contingency(table)
                    except ValueError:
                        chi2, pnaive = np.nan, np.nan
                else:
                    chi2, pnaive = np.nan, np.nan
                r["p_naive_pooled"] = float(pnaive) if pnaive == pnaive else np.nan
                r["prop_pooled_a"] = ca / ta if ta else np.nan
                r["prop_pooled_b"] = cb / tb if tb else np.nan
                rows.append(r)
    d = pd.DataFrame(rows)
    if len(d):                                    # 每个 数据集×层级×contrast 内做 BH
        d["q_bh"] = np.nan
        for _, idx in d.groupby(["level", "group_a", "group_b"]).groups.items():
            pv = d.loc[idx, "p_value"].to_numpy()
            order = np.argsort(pv); m = len(pv); q = np.empty(m); prev = 1.0
            for rank, i in enumerate(order[::-1], 1):
                prev = min(prev, pv[i] * m / (m - rank + 1)); q[i] = prev
            d.loc[idx, "q_bh"] = q
    return d


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=Path("."),
                    help="数据根目录；见 data/DATA_SOURCES.md")
    ap.add_argument("--out", type=Path, default=Path("results_methods/floor_reassessment"))
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S")
    a.out.mkdir(parents=True, exist_ok=True)

    loaded, dropped = [], []
    for fn in (load_D1, load_D2, load_D3, load_D4, load_D5):
        try:
            r = fn(a.root)
        except Exception as e:                     # 记录，不静默跳过
            dropped.append(dict(loader=fn.__name__, reason=f"{type(e).__name__}: {e}"))
            logger.warning("%s 失败: %s", fn.__name__, e); continue
        if r is None:
            dropped.append(dict(loader=fn.__name__, reason="必需字段或文件缺失"))
            logger.warning("%s: 缺字段/文件，按预注册 §9 丢弃", fn.__name__); continue
        loaded.extend(r if isinstance(r, list) else [r])

    allres = []
    for ds in loaded:
        if not ds["contrasts"]:
            dropped.append(dict(loader=ds["name"], reason="无可用的分组对比"))
            logger.warning("%s: 无 contrast，丢弃", ds["name"]); continue
        logger.info("%s | %s | donor=%s group=%s | contrasts=%s",
                    ds["name"], ds["tissue"], ds["donor"], ds["group"], ds["contrasts"])
        d = analyse(ds)
        logger.info("   -> %d 个 contrast", len(d))
        if len(d):
            allres.append(d)

    res = pd.concat(allres, ignore_index=True) if allres else pd.DataFrame()
    res.to_csv(a.out / "contrasts.csv", index=False)
    json.dump(dropped, open(a.out / "dropped.json", "w"), indent=2, ensure_ascii=False)
    logger.info("总 contrast: %d，丢弃记录 %d 条 -> %s", len(res), len(dropped), a.out)


if __name__ == "__main__":
    main()
