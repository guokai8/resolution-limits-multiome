#!/usr/bin/env python3
"""Layer 3 外部复现汇总：7 个细胞类 + n 梯队阳性对照 + Mantel-Haenszel 合并。"""
import json, glob, numpy as np

def load(p): return json.load(open(p))

def mh(ds, bg="all"):
    Rs = Ss = 0.0; PR = PS = QR = QS = 0.0
    for d in ds:
        a = d["prox_links"]; b = d["n_links"] - a
        c = d[f"prox_tested_{bg}"]; dd = d[f"n_tested_{bg}"] - c
        n = a + b + c + dd
        R = a*dd/n; S = b*c/n; Rs += R; Ss += S
        P = (a+dd)/n; Q = (b+c)/n
        PR += P*R; PS += P*S; QR += Q*R; QS += Q*S
    OR = Rs/Ss
    se = np.sqrt(PR/(2*Rs**2) + (PS+QR)/(2*Rs*Ss) + QS/(2*Ss**2))
    return OR, OR*np.exp(-1.96*se), OR*np.exp(1.96*se)

def row(d, lab):
    pl = 100*d["prox_links"]/max(d["n_links"], 1)
    return (f"{lab:12s} {d['n']:6d} {d['n_tested_all']:10d} {d['n_links']:7d} "
            f"{pl:6.2f} {d['OR']:7.3f} [{d['lo']:5.3f},{d['hi']:6.3f}]")

HDR = f"{'':12s} {'n':>6s} {'tested':>10s} {'links':>7s} {'prox%':>6s} {'OR':>7s} {'95% CI':>16s}"

classes = ["Oligo", "MG", "ExN", "InN", "Astro", "OPC", "VC"]
coarse = []
print("== 7 个细胞类，n=150，ATAC 5564 / RNA 5265，口径全对齐 ==")
print(HDR)
for c in classes:
    d = load(f"results_methods/nabec_L3/{c}_n150p508match_seed0.json"); coarse.append(d)
    print(row(d, c))

fine = [(l, f"results_methods/nabec_L3/{k}_n150p508match_fine_seed0.json")
        for l, k in [("Oligo cl18","18"),("MG cl21","21"),("ExN cl7","7"),("InN cl37","37")]]
print("\n== 4 个细亚型（粒度与原队列一致） ==")
print(HDR)
fds = []
for lab, p in fine:
    d = load(p); fds.append(d); print(row(d, lab))

print("\n== 阳性对照：深度锁死 5564/5265，只提细胞数 ==")
print(HDR)
for c in ["Oligo", "ExN"]:
    d0 = load(f"results_methods/nabec_L3/{c}_n150p508match_seed0.json")
    print(row(d0, c))
    for nn in [600, 2400, 5000]:
        p = f"results_methods/nabec_L3/{c}_n{nn}p508match_n{nn}_seed0.json"
        try: print(row(load(p), c))
        except FileNotFoundError: print(f"{c:12s} {nn:6d}  (pending)")

print("\n== Mantel-Haenszel 合并 ==")
for lab, ds in [("7 个细胞类 (n=150)", coarse), ("4 个细亚型 (n=150)", fds),
                ("全部 11 个 (n=150)", coarse+fds)]:
    o, lo, hi = mh(ds, "all"); ov, lov, hiv = mh(ds, "var")
    print(f"{lab:22s} 手稿口径 {o:6.3f} [{lo:5.3f},{hi:5.3f}]   "
          f"varpeaks {ov:6.3f} [{lov:5.3f},{hiv:5.3f}]")
print("\n原队列 (手稿 Fig.4c): tested 3029251, links 33227, prox% 0.63, OR 0.663 [0.568,0.768]")
print("全深度 cellranger-arc 参照区间: 2.45 - 4.01")
