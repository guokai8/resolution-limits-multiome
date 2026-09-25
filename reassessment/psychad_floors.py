#!/usr/bin/env python3
"""PsychAD（MSSM + RADC）的组成与表达技术下限。

口径与手稿完全一致（照 p5_10 / p5_14）：
  组成: z = |f1-f2| / sqrt(p(1-p))，多项式期望 sqrt(2/n_eff)，
        过度离散 = RMS of z/sqrt(2/n_eff)
  表达: median |log2(CPM1+1) - log2(CPM2+1)|，基因取并集 (CPM>0 任一)，≥200 基因
        零模型 = 合并 profile 按各自深度多项式重抽样
  标度: log(floor) ~ log(n_eff) OLS，按 donor 聚类 bootstrap

⚠️ PsychAD 的重复是【同一管悬液分装成两份并行上机】（Sci Data 2025 方法：
   "2 aliquots of 60,000 pooled nuclei ... processed in parallel"），
   因此它测的是【上机+建库+测序】那一层，不含核分离与组织取样。
   这正是与 Ruf et al. 队列（含组织分装+核分离+不同实验批次）对比的价值所在。
"""
from __future__ import annotations
import argparse, json, logging, time
from pathlib import Path
import numpy as np, pandas as pd, h5py, scipy.sparse as sp

logger = logging.getLogger("psychad")
MIN_CELLS = 50

def read_obs(path: Path, cols):
    out = {}
    with h5py.File(path, "r") as f:
        o = f["obs"]
        for c in cols:
            if c not in o: continue
            g = o[c]
            if isinstance(g, h5py.Group):
                out[c] = np.asarray(g["categories"][:]).astype(str)[g["codes"][:]]
            else:
                v = np.asarray(g[:]); out[c] = v.astype(str) if v.dtype.kind in "SOU" else v
    return pd.DataFrame(out)

def fit(x, y, g, rng, nb=2000):
    ok = np.isfinite(x)&np.isfinite(y)&(x>0)&(y>0)
    lx, ly, gg = np.log(x[ok]), np.log(y[ok]), g[ok]
    b = float(np.polyfit(lx, ly, 1)[0])
    u, inv = np.unique(gg, return_inverse=True)
    ix = [np.flatnonzero(inv==i) for i in range(u.size)]
    bs = [np.polyfit(lx[s], ly[s], 1)[0] for s in
          (np.concatenate([ix[i] for i in rng.integers(0,u.size,u.size)]) for _ in range(nb))]
    return [b, float(np.percentile(bs,2.5)), float(np.percentile(bs,97.5))]

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5ad", type=Path, required=True)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--level", default="subclass")
    ap.add_argument("--out", type=Path, default=Path("results_methods/psychad_floors"))
    ap.add_argument("--seed", type=int, default=2026)
    a = ap.parse_args(); a.out.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S")
    rng = np.random.default_rng(a.seed); t0=time.time()

    obs = read_obs(a.h5ad, ["barcodekey","donor_id",a.level,"AD_status"])
    obs["mid"] = obs.barcodekey.str.extract(r"^.+?-(\d+)-[ACGT]+-\d+$")[0]
    obs = obs.dropna(subset=["mid", a.level])
    logger.info("%s: %d 核, %d donor, %d 个 %s",
                a.cohort, len(obs), obs.donor_id.nunique(), obs[a.level].nunique(), a.level)

    # 每 donor 取细胞数最多的两个 aliquot，各需 >=50
    n = obs.groupby(["donor_id","mid"]).size().rename("n").reset_index()
    n = n[n.n >= MIN_CELLS]
    top2 = n.sort_values("n", ascending=False).groupby("donor_id").head(2)
    pairs = top2.groupby("donor_id").filter(lambda g: len(g)==2)
    sel = set(zip(pairs.donor_id, pairs.mid))
    logger.info("  重复对: %d 个 donor", pairs.donor_id.nunique())

    keep = np.array([(d,m) in sel for d,m in zip(obs.donor_id, obs.mid)])
    sub = obs[keep]
    types = sorted(sub[a.level].unique())
    cnt = sub.groupby(["donor_id","mid",a.level]).size()
    tot = sub.groupby(["donor_id","mid"]).size()

    rows=[]
    for donor, g in pairs.groupby("donor_id"):
        m1, m2 = sorted(g.mid)
        n1, n2 = int(tot[(donor,m1)]), int(tot[(donor,m2)])
        neff = 2/(1/n1+1/n2)
        for ct in types:
            c1 = int(cnt.get((donor,m1,ct),0)); c2 = int(cnt.get((donor,m2,ct),0))
            f1, f2 = c1/n1, c2/n2
            p = (c1+c2)/(n1+n2)
            if p <= 0 or p >= 1: continue
            z = abs(f1-f2)/np.sqrt(p*(1-p))
            if z <= 0: continue
            rows.append(dict(donor=donor, celltype=ct, n1=n1, n2=n2, n_eff=neff, p=p, z=z))
    d = pd.DataFrame(rows)
    d.to_csv(a.out/f"L1_{a.cohort}_pairs.csv", index=False)

    od = d.z/np.sqrt(2.0/d.n_eff)
    b = fit(d.n_eff.to_numpy(), d.z.to_numpy(), d.donor.to_numpy(), rng)
    res = dict(cohort=a.cohort, level=a.level, n_pairs=int(d.donor.nunique()),
               n_obs=int(len(d)), n_types=len(types),
               replicate_unit="two aliquots of one pooled nuclear suspension, processed in parallel",
               overdispersion_rms=float(np.sqrt((od**2).mean())),
               overdispersion_median=float(od.median()),
               exponent=b,
               n_eff_median=float(d.n_eff.median()))
    json.dump(res, open(a.out/f"L1_{a.cohort}_summary.json","w"), indent=2, ensure_ascii=False)
    logger.info("组成下限 | 过度离散 RMS %.2f× (中位 %.2f) | 指数 %.3f (%.3f, %.3f) | n_eff 中位 %.0f | %.0fs",
                res["overdispersion_rms"], res["overdispersion_median"], *b, res["n_eff_median"], time.time()-t0)

if __name__ == "__main__":
    main()
