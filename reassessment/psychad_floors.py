#!/usr/bin/env python3
"""Composition floor in the PsychAD replicates (MSSM and RADC).

Conventions are identical to p5_10 and p5_14 so the numbers are comparable:
  composition  z = |f1-f2| / sqrt(p(1-p)), multinomial expectation
               sqrt(2/n_eff), overdispersion = RMS of z / sqrt(2/n_eff)
  scaling      OLS of log(z) on log(n_eff), bootstrap clustered by donor

What a PsychAD replicate actually repeats matters for interpreting the result.
These are two aliquots of ONE pooled nuclear suspension run in parallel (Sci
Data 2025: "2 aliquots of 60,000 pooled nuclei ... processed in parallel"), so
they capture loading, library prep and sequencing -- and nothing before that.
Tissue sampling and nuclear dissociation are shared between the two members and
cannot contribute to the measured spread.

That is precisely why they are worth measuring. Contrasting them with the Ruf
et al. cohort, whose replicates do span tissue sub-sampling and dissociation,
is what separates the part of the floor that comes from the bench from the part
that comes from the instrument.
"""
from __future__ import annotations
import argparse, json, logging, time
from pathlib import Path
import numpy as np, pandas as pd, h5py, scipy.sparse as sp

logger = logging.getLogger("psychad")
MIN_CELLS = 50

def read_obs(path: Path, cols):
    """Read selected obs columns from an .h5ad without loading X.

    AnnData stores a categorical as a group of codes plus categories, and a
    plain array otherwise, so both shapes have to be handled. Reading obs alone
    keeps the memory cost independent of the matrix size.
    """
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
    """OLS slope of log(y) on log(x), with a bootstrap clustered by g.

    Resampling donors rather than observations: the pairs from one donor share
    a suspension, so treating them as independent would understate the
    interval. Returns [slope, 2.5th percentile, 97.5th percentile].
    """
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

    # Take the two largest aliquots per donor, each needing at least 50 nuclei.
    # Largest two rather than all pairs: a donor with three aliquots would
    # otherwise contribute three correlated pairs and be weighted threefold.
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
            # p of 0 or 1 means the cell type is absent from both aliquots or
            # fills them; sqrt(p(1-p)) is then 0 and z is undefined.
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
