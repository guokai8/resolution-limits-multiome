#!/usr/bin/env python3
"""Expression floor in the PsychAD replicates, matching p5_14 exactly.

    floor = median over genes of |log2(CPM1+1) - log2(CPM2+1)|

taken over the union of genes with CPM > 0 in either aliquot, requiring at
least 200 of them.

The matched null resamples both members multinomially from their pooled profile
at their own observed depths, then runs the identical floor computation. It
therefore differs from the observation in one respect only: it has no process
variation. The ratio of the two is the quantity of interest.

As in psychad_floors.py, these replicates span loading and library prep but not
dissociation, which bounds what the measured floor can include.
"""
from __future__ import annotations
import argparse, json, logging, time
from pathlib import Path
import numpy as np, pandas as pd, h5py, scipy.sparse as sp

logger = logging.getLogger("L2")
MIN_CELLS, MIN_TYPE = 50, 10

def read_obs(path, cols):
    out={}
    with h5py.File(path,"r") as f:
        o=f["obs"]
        for c in cols:
            if c not in o: continue
            g=o[c]
            out[c]=(np.asarray(g["categories"][:]).astype(str)[g["codes"][:]]
                    if isinstance(g,h5py.Group) else np.asarray(g[:]).astype(str))
    return pd.DataFrame(out)

def floor_from_counts(c1,c2,min_cpm=0.0):
    """The expression floor for one pair of pseudobulk count vectors.

    Union rule, not intersection: a gene counts when CPM > 0 in EITHER member.
    An intersection would silently drop the genes that differ most between the
    two libraries, which are exactly the ones the floor is meant to measure.

    Returns NaN rather than raising when the pair is unusable, since scanning
    many pairs will always turn up some that are.
    """
    t1,t2=c1.sum(),c2.sum()
    if t1<=0 or t2<=0: return np.nan
    p1,p2=c1/t1*1e6, c2/t2*1e6
    keep=(p1>min_cpm)|(p2>min_cpm)
    if keep.sum()<200: return np.nan
    return float(np.median(np.abs(np.log2(p1[keep]+1)-np.log2(p2[keep]+1))))

def main():
    ap=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5ad",type=Path,required=True); ap.add_argument("--cohort",required=True)
    ap.add_argument("--level",default="subclass"); ap.add_argument("--chunk",type=int,default=20000)
    ap.add_argument("--out",type=Path,default=Path("results_methods/psychad_floors"))
    ap.add_argument("--seed",type=int,default=2026)
    a=ap.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    logging.basicConfig(level=logging.INFO,format="%(asctime)s %(message)s",datefmt="%H:%M:%S")
    rng=np.random.default_rng(a.seed); t0=time.time()

    obs=read_obs(a.h5ad,["barcodekey","donor_id",a.level])
    obs["mid"]=obs.barcodekey.str.extract(r"^.+?-(\d+)-[ACGT]+-\d+$")[0]
    obs=obs.dropna(subset=["mid",a.level]).reset_index(drop=True)
    n=obs.groupby(["donor_id","mid"]).size().rename("n").reset_index()
    n=n[n.n>=MIN_CELLS]
    pr=n.sort_values("n",ascending=False).groupby("donor_id").head(2)
    pr=pr.groupby("donor_id").filter(lambda g:len(g)==2)
    sel={(d,m) for d,m in zip(pr.donor_id,pr.mid)}
    keep=np.array([(d,m) in sel for d,m in zip(obs.donor_id,obs.mid)])
    types=sorted(obs[a.level].unique()); t_map={t:i for i,t in enumerate(types)}
    grp=sorted(sel); g_map={g:i for i,g in enumerate(grp)}
    # Encode (aliquot, cell type) as one integer key per nucleus, -1 for the
    # nuclei not in any selected pair. A single key lets the pseudobulk sum be
    # done as one sparse matrix product per chunk instead of a Python loop.
    key=np.full(len(obs),-1,np.int64)
    idx=np.flatnonzero(keep)
    key[idx]=(np.array([g_map[(d,m)] for d,m in zip(obs.donor_id[idx],obs.mid[idx])])*len(types)
              + obs.loc[idx,a.level].map(t_map).to_numpy())
    NK=len(grp)*len(types)
    logger.info("%s: %d 对, %d 组合, %d 细胞参与",a.cohort,pr.donor_id.nunique(),NK,len(idx))

    with h5py.File(a.h5ad,"r") as f:
        NG=int(f["X"].attrs["shape"][1])
        dptr,ddat,didx=f["X/indptr"],f["X/data"],f["X/indices"]
        acc=np.zeros((NK,NG),np.float32); nc=np.zeros(NK,np.int64)
        # Stream X in row chunks and accumulate group sums. The indptr slice is
        # rebased to the chunk (ip-b) so each chunk builds a standalone CSR.
        for s0 in range(0,len(obs),a.chunk):
            s1=min(s0+a.chunk,len(obs)); k=key[s0:s1]; good=k>=0
            if not good.any(): continue
            ip=dptr[s0:s1+1]; b,e=int(ip[0]),int(ip[-1])
            M=sp.csr_matrix((ddat[b:e],didx[b:e],ip-b),shape=(s1-s0,NG))
            kk=k[good]
            # G is a group-membership indicator; G @ M sums each group's rows.
            G=sp.csr_matrix((np.ones(kk.size,np.float32),(kk,np.flatnonzero(good))),shape=(NK,s1-s0))
            acc+=np.asarray((G@M).todense(),dtype=np.float32); nc+=np.bincount(kk,minlength=NK)
            if (s0//a.chunk)%40==0: logger.info("  %d/%d 行 (%.0fs)",s1,len(obs),time.time()-t0)

    rows=[]
    for donor,g in pr.groupby("donor_id"):
        m1,m2=sorted(g.mid); i1,i2=g_map[(donor,m1)],g_map[(donor,m2)]
        for ti,ct in enumerate(types):
            r1,r2=i1*len(types)+ti, i2*len(types)+ti
            n1,n2=int(nc[r1]),int(nc[r2])
            if min(n1,n2)<MIN_TYPE: continue
            c1,c2=acc[r1].astype(np.float64),acc[r2].astype(np.float64)
            f_=floor_from_counts(c1,c2)
            if not np.isfinite(f_): continue
            # The matched null: draw both members from the pooled profile at
            # their own observed totals, so depth and gene set are preserved and
            # only the process variation is removed.
            pool=c1+c2; p=pool/pool.sum(); nz=p>0
            s1=np.zeros_like(c1); s2=np.zeros_like(c2)
            s1[nz]=rng.multinomial(int(c1.sum()),p[nz]); s2[nz]=rng.multinomial(int(c2.sum()),p[nz])
            f0=floor_from_counts(s1,s2)
            rows.append(dict(donor=donor,celltype=ct,n1=n1,n2=n2,n_eff=2/(1/n1+1/n2),
                             floor=f_,floor_null=f0,ratio=f_/f0 if f0 and np.isfinite(f0) else np.nan))
    d=pd.DataFrame(rows); d.to_csv(a.out/f"L2_{a.cohort}_pairs.csv",index=False)

    def fit(x,y,g,nb=1000):
        ok=np.isfinite(x)&np.isfinite(y)&(x>0)&(y>0)
        lx,ly,gg=np.log(x[ok]),np.log(y[ok]),g[ok]; b=float(np.polyfit(lx,ly,1)[0])
        u,inv=np.unique(gg,return_inverse=True); ix=[np.flatnonzero(inv==i) for i in range(u.size)]
        bs=[np.polyfit(lx[s],ly[s],1)[0] for s in
            (np.concatenate([ix[i] for i in rng.integers(0,u.size,u.size)]) for _ in range(nb))]
        return [b,float(np.percentile(bs,2.5)),float(np.percentile(bs,97.5))]
    g=d.donor.to_numpy(); rt=d.ratio.dropna()
    res=dict(cohort=a.cohort,n_pairs=int(d.donor.nunique()),n_obs=int(len(d)),
             exponent=fit(d.n_eff.to_numpy(),d.floor.to_numpy(),g),
             null_exponent=fit(d.n_eff.to_numpy(),d.floor_null.to_numpy(),g),
             ratio_median=float(rt.median()),
             ratio_iqr=[float(rt.quantile(.25)),float(rt.quantile(.75))],
             floor_median=float(d.floor.median()))
    json.dump(res,open(a.out/f"L2_{a.cohort}_summary.json","w"),indent=2,ensure_ascii=False)
    logger.info("表达下限 | 指数 %.3f (%.3f,%.3f) | 零模型 %.3f | 观测/零模型 中位 %.2f (IQR %.2f-%.2f) | %.1f 分",
                *res["exponent"],res["null_exponent"][0],res["ratio_median"],*res["ratio_iqr"],(time.time()-t0)/60)

if __name__=="__main__": main()
