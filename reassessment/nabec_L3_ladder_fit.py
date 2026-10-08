#!/usr/bin/env python3
"""The nucleus ladder: pool seeds, fit OR against nucleus number, find where
it crosses unity.

Each rung of the ladder was run under three downsampling seeds. The seeds are
pooled by Mantel-Haenszel here, and the fit is then used to read off the
nucleus number at which the diagnostic would cross 1.00 and the number at which
it would reach 2.45, the lower edge of the full-depth reference range.

Those crossing points are extrapolations from a fitted line, not measurements,
and the paper quotes them as such.
"""
import json, glob, re, numpy as np
from collections import defaultdict

def mh(ds, bg="all"):
    """Mantel-Haenszel pooling with the Robins-Breslow-Greenland variance.

    Each stratum contributes a 2x2 table: promoter-proximal or not, among
    detected links against among all tested pairs. `bg` selects the background
    convention -- "all" counts every tested pair, "var" excludes peaks that
    cannot be called at the depth analysed.
    """
    Rs=Ss=0.0; PR=PS=QR=QS=0.0
    for d in ds:
        a=d["prox_links"]; b=d["n_links"]-a
        c=d[f"prox_tested_{bg}"]; dd=d[f"n_tested_{bg}"]-c; n=a+b+c+dd
        R=a*dd/n; S=b*c/n; Rs+=R; Ss+=S
        P=(a+dd)/n; Q=(b+c)/n
        PR+=P*R; PS+=P*S; QR+=Q*R; QS+=Q*S
    OR=Rs/Ss; se=np.sqrt(PR/(2*Rs**2)+(PS+QR)/(2*Rs*Ss)+QS/(2*Ss**2))
    return OR, OR*np.exp(-1.96*se), OR*np.exp(1.96*se), se

# The run files encode the rung and the seed in the filename, so group them
# back into (cell type, nucleus number) with the three seeds together.
g=defaultdict(list)
for f in glob.glob('results_methods/nabec_L3/*ladder_n*_s*.json'):
    j=json.load(open(f)); m=re.search(r'_n(\d+)_s(\d)',f)
    g[(j['celltype'], int(m.group(1)))].append(j)

NS=sorted({n for _,n in g}); CTS=sorted({c for c,_ in g})
print("== 三种子 MH 合并 ==")
print(f"{'type':6s} {'n':>6s} {'links':>7s} {'prox':>5s} {'MH OR':>7s} {'95% CI':>16s}")
pts=defaultdict(list)
for ct in CTS:
    for n in NS:
        ds=g[(ct,n)]
        if not ds: continue
        o,lo,hi,se=mh(ds)
        lk=sum(d['n_links'] for d in ds); px=sum(d['prox_links'] for d in ds)
        print(f"{ct:6s} {n:6d} {lk:7d} {px:5d} {o:7.3f} [{lo:5.3f},{hi:6.3f}]")
        pts[ct].append((n,o,se))
    print()

print("== 对数线性拟合 log(OR) ~ log(n)，用 1/se^2 加权，仅取 n>=600 的上升段 ==")
for ct in CTS:
    # Fit only the rising arm. Below 600 nuclei the curve is flat near the
    # floor, and including it would bias the slope and the crossing point.
    P=[(n,o,se) for n,o,se in pts[ct] if n>=600]
    x=np.log([p[0] for p in P]); y=np.log([p[1] for p in P]); w=1/np.array([p[2] for p in P])**2
    # Weight by inverse variance: the low rungs carry far wider intervals and
    # should not pull the line as hard as the well-determined high rungs.
    b,a=np.polyfit(x,y,1,w=np.sqrt(w))
    n1=np.exp(-a/b); n245=np.exp((np.log(2.45)-a)/b)
    r=np.corrcoef(x,y)[0,1]
    print(f"{ct:6s} log(OR) = {a:.3f} + {b:.3f} log(n)   r = {r:.3f}")
    print(f"        OR = 1.00 处 n ≈ {n1:,.0f}      OR = 2.45（全深度参照下沿）处 n ≈ {n245:,.0f}")
