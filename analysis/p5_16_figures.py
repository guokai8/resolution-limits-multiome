#!/usr/bin/env python3
"""
Superseded matplotlib figure script, kept for provenance
=======================================================
The published figures come from the R/ggplot2 scripts under figures/. This file
is retained because it is the version that produced the figures circulated
during the analysis, and because its panel-by-panel notes record why several
panels are laid out the way they are.

Figure 2 is the paper's central visual argument: one layout, two opposite
results.
  a  expression layer, two cohorts, fitted lines PARALLEL -- same exponent
     (-0.507 against -0.505), different intercept (4.79 against 3.65).
     Not COINCIDENT. "The exponent transfers, the coefficient does not" is the
     claim, and an annotation that says otherwise would be seen through at once.
  b  composition layer, the two cohorts' lines DIVERGE (-0.152 against -0.584)
  c  overdispersion agrees between cohorts (4.28x against 3.90x) -- the half
     that does replicate
  Panels a and b must share axes and scale, or the layout itself will absorb
  the parallel-against-divergent contrast the figure exists to show.

Figure 5 is the diagnostic's reference range.
  a  forest plot of odds ratios over six link sets, log axis, unity marked
  b  window sensitivity: the same data gives 8.84 at +/-1.7 Mb and 5.40 at 1 Mb
  c  28 external multiome libraries: linkage against (N x depth) on log-log,
     exponent +0.42, no saturation

matplotlib only, no seaborn: fewer dependencies, and it ran on a 3 GB machine.
A missing input skips its panel with a warning rather than aborting, because
the analyses finished at different times.

Usage:
  python3 p5_16_figures.py --results results --out figures
"""

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Fixed cohort colours, consistent across every figure
C_ULM, C_SEA, C_NULL = "#1f4e79", "#c0504d", "#7f7f7f"

# Constants already published or already settled upstream
ULM = dict(L1_b=-0.152, L1_ci=(-0.340, 0.111), L2_b=-0.507, L2_ci=(-0.569, -0.441),
           L2_a=4.79, od=4.28)
SEA = dict(L1_b=-0.584, L1_ci=(-0.838, -0.369), L2_b=-0.505, L2_ci=(-0.524, -0.483),
           L2_a=3.65, od=3.90)
OR_SETS = [  # (标签, OR, lo, hi)  —— 升序，本项目置于最下以突出其在 1 的另一侧
    ("pbmc_granulocyte_10k", 2.451, 2.399, 2.506),
    ("pbmc_granulocyte_3k", 3.074, 2.961, 3.196),
    ("lymph_node_lymphoma_14k", 3.575, 3.466, 3.693),
    ("human_brain_3k", 4.012, 3.906, 4.125),
    ("Published AD multiome", 5.40, 5.30, 5.51),
    ("This study (150 nuclei)", 0.663, 0.568, 0.768),
]


def col(df, *cands):
    """Column names differ between the analysis scripts (floor against
    median_abs_log2FC), so try each in turn."""
    if df is None:
        return None
    for c in cands:
        if c in df.columns:
            return df[c].values
    print(f"    ⚠️ 未找到列 {cands}，实有 {list(df.columns)}")
    return None


def loglog_panel(ax, datasets, ylab, title, ann):
    """datasets: [(label, color, x, y, b, ci)]. Shared axes and scale, which is
    what makes the two panels comparable."""
    for lab, cc, x, y, b, ci in datasets:
        if x is not None and y is not None:
            m = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
            x, y = x[m], y[m]
        col = cc
        if x is not None and len(x):
            ax.scatter(x, y, s=4, alpha=.12, color=col, edgecolors="none",
                       rasterized=True)
        xs = np.array([np.nanmin(x), np.nanmax(x)]) if x is not None and len(x) \
            else np.array([1e2, 1e4])
        # Median of the point cloud as the intercept; used only to place the
        # illustrative line, never as an estimate
        if x is not None and len(x):
            a = np.exp(np.median(np.log(y) - b * np.log(x)))
        else:
            a = 1.0
        ax.plot(xs, a * xs ** b, color=col, lw=2.2,
                label=f"{lab}  b = {b:+.3f} [{ci[0]:+.3f}, {ci[1]:+.3f}]")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Effective number of nuclei ($n_{eff}$)")
    ax.set_ylabel(ylab)
    ax.set_title(title, loc="left", fontsize=11, fontweight="bold")
    ax.legend(fontsize=7.5, frameon=False, loc="lower left")
    ax.text(.97, .95, ann, transform=ax.transAxes, ha="right", va="top",
            fontsize=9, color="#333")


def fig2(res, out):
    f = plt.figure(figsize=(12, 4.2))
    gs = f.add_gridspec(1, 3, width_ratios=[1, 1, .75], wspace=.42)

    # ---- a: expression layer, the two cohorts
    ax = f.add_subplot(gs[0])
    s2 = _first(res, "seaad_L2_obs.csv", "seaad_L2_pairs_genelevel.csv")
    u2 = _read(res, "p5_floor_scaling_pairs.csv")   # Ulm，名称若不同请改
    loglog_panel(ax, [
        ("Motor cortex", C_ULM,
         col(u2, "n_eff"), col(u2, "floor", "median_abs_log2FC"),
         ULM["L2_b"], ULM["L2_ci"]),
        ("Seattle atlas", C_SEA,
         col(s2, "n_eff"), col(s2, "floor", "median_abs_log2FC"),
         SEA["L2_b"], SEA["L2_ci"]),
    ], "Expression floor  (median |Δlog₂CPM|)",
        "a   Expression layer", "parallel:\nsame slope,\ndifferent intercept")

    # ---- b: composition layer, the two cohorts diverge (same scale as a)
    ax2 = f.add_subplot(gs[1])
    s1 = _read(res, "seaad_L1_pairs.csv")

    u1 = _first(res, "p5_L1_scaling_pairs.csv", "p5_two_layer_L1_pairs.csv")
    loglog_panel(ax2, [
        ("Motor cortex", C_ULM,
         col(u1, "n_eff"), col(u1, "z"), ULM["L1_b"], ULM["L1_ci"]),
        ("Seattle atlas", C_SEA,
         col(s1, "n_eff"), col(s1, "z"), SEA["L1_b"], SEA["L1_ci"]),
    ], "Composition discrepancy  (standardised)",
        "b   Composition layer", "divergent:\nslopes differ")

    # ---- c: overdispersion agrees -- the half that does replicate
    ax3 = f.add_subplot(gs[2])
    ax3.bar([0, 1], [ULM["od"], SEA["od"]], color=[C_ULM, C_SEA], width=.55)
    ax3.axhline(1, ls="--", lw=1, color="k")
    ax3.text(1.02, 1, "multinomial\nexpectation", va="center", fontsize=7.5,
             transform=ax3.get_yaxis_transform())
    ax3.set_xticks([0, 1]); ax3.set_xticklabels(["Motor\ncortex", "Seattle\natlas"])
    ax3.set_ylabel("Composition overdispersion (×)")
    ax3.set_title("c   Overdispersion", loc="left",
                  fontsize=11, fontweight="bold")
    for i, v in enumerate([ULM["od"], SEA["od"]]):
        ax3.text(i, v + .12, f"{v:.2f}×", ha="center", fontsize=9)

    _save(f, out, "Fig2_transferability")


def fig5(res, out):
    f = plt.figure(figsize=(12, 4.2))
    gs = f.add_gridspec(1, 3, width_ratios=[1.15, .7, .9], wspace=.38)

    # ---- a: odds-ratio forest plot
    ax = f.add_subplot(gs[0])
    ys = np.arange(len(OR_SETS))[::-1]
    for y, (lab, o, lo, hi) in zip(ys, OR_SETS):
        col = "#c0504d" if o < 1 else "#1f4e79"
        ax.plot([lo, hi], [y, y], color=col, lw=2.4)
        ax.plot(o, y, "o", color=col, ms=6)
        ax.text(hi * 1.06, y, f"{o:.2f}", va="center", fontsize=8, color=col)
    ax.axvline(1, ls="--", color="k", lw=1.1)
    ax.axvspan(2.4, 5.6, color="#1f4e79", alpha=.07)
    ax.set_xscale("log")
    ax.set_yticks(ys); ax.set_yticklabels([s[0] for s in OR_SETS], fontsize=8)
    ax.set_xlabel("Promoter-enrichment odds ratio  (log scale)")
    ax.set_title("a   Reference range across link sets", loc="left",
                 fontsize=11, fontweight="bold")

    # ---- b: window sensitivity
    ax2 = f.add_subplot(gs[1])
    ax2.bar([0, 1], [5.40, 8.84], color=["#1f4e79", "#9dc3e6"], width=.55)
    ax2.set_xticks([0, 1]); ax2.set_xticklabels(["±1 Mb\n(matched)", "±1.7 Mb\n(own)"])
    ax2.set_ylabel("Odds ratio, same link set")
    ax2.set_title("b   Window matters", loc="left", fontsize=11, fontweight="bold")
    for i, v in enumerate([5.40, 8.84]):
        ax2.text(i, v + .2, f"{v:.2f}", ha="center", fontsize=9)

    # ---- c: linkage saturation curve
    ax3 = f.add_subplot(gs[2])
    lib = _read(res, "seaad_L3_libraries.csv")
    if lib is not None and {"n_cells", "link", "atac_frag"} <= set(lib.columns):
        s = lib.dropna(subset=["n_cells", "link", "atac_frag"])
        x = s.n_cells * s.atac_frag
        ax3.scatter(x, s.link, s=26, color=C_SEA, alpha=.8, edgecolors="none")
        xs = np.array([x.min(), x.max()])
        b = 0.421
        a = np.exp(np.median(np.log(s.link) - b * np.log(x)))
        ax3.plot(xs, a * xs ** b, color="k", lw=1.8,
                 label=f"b = +{b:.2f}  (no plateau)")
        ax3.legend(fontsize=8, frameon=False)
    else:
        ax3.text(.5, .5, "seaad_L3_libraries.csv\nnot found", ha="center",
                 va="center", transform=ax3.transAxes, color="#999")
    ax3.set_xscale("log"); ax3.set_yscale("log")
    ax3.set_xlabel("Nuclei × fragments per nucleus")
    ax3.set_ylabel("Feature linkages detected")
    ax3.set_title("c   No saturation (independent pipeline)", loc="left",
                  fontsize=11, fontweight="bold")

    _save(f, out, "Fig5_diagnostic")


def fig3(res, out):
    """Null-model comparison and within-type stratification. The paper's most
    important self-check: b near -0.5 is what the null model predicts, so the
    exponent alone is not evidence of anything."""
    f = plt.figure(figsize=(12, 4.2))
    gs = f.add_gridspec(1, 3, width_ratios=[1, .8, 1.15], wspace=.42)
    o = _read(res, "seaad_L2_obs.csv"); n = _read(res, "seaad_L2_null.csv")

    # a: observed against the null model
    ax = f.add_subplot(gs[0])
    for d, lab, c, b in ((o, "Observed", C_SEA, -0.505),
                         (n, "Sampling null", C_NULL, -0.482)):
        if d is None:
            continue
        x, y = d.n_eff.values, d.floor.values
        m = (x > 0) & (y > 0)
        ax.scatter(x[m], y[m], s=4, alpha=.10, color=c, edgecolors="none",
                   rasterized=True)
        xs = np.array([x[m].min(), x[m].max()])
        a = np.exp(np.median(np.log(y[m]) - b * np.log(x[m])))
        ax.plot(xs, a * xs ** b, color=c, lw=2.2, label=f"{lab}  b = {b:+.3f}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Effective number of nuclei ($n_{eff}$)")
    ax.set_ylabel("Expression floor")
    ax.set_title("a   Observed vs sampling null", loc="left",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=8, frameon=False, loc="lower left")

    # b: the per-pair ratio, which is where all the information is
    ax2 = f.add_subplot(gs[1])
    if o is not None and n is not None:
        j = o.merge(n, on=["donor", "ct"], suffixes=("_o", "_m"))
        r = (j.floor_o / j.floor_m).replace([np.inf, -np.inf], np.nan).dropna()
        ax2.hist(r, bins=45, color=C_SEA, alpha=.85)
        ax2.axvline(1, ls="--", color="k", lw=1.2)
        ax2.axvline(r.median(), color="#1f4e79", lw=2)
        ax2.text(r.median() * 1.04, ax2.get_ylim()[1] * .9,
                 f"median\n{r.median():.2f}×", fontsize=9, color="#1f4e79")
    ax2.set_xlabel("Observed / null floor, per pair")
    ax2.set_ylabel("Pairs")
    ax2.set_title("b   Excess over pure sampling", loc="left",
                  fontsize=11, fontweight="bold")

    # c: forest plot of within-type slopes
    ax3 = f.add_subplot(gs[2])
    st = _read(res, "seaad_L2_by_celltype.csv")
    if st is not None:
        st = st.sort_values("b")
        ys = np.arange(len(st))
        # The widest-span types are the best determined, so they are labelled
        # separately -- including the unexplained steep neuronal slopes
        wide = st.span >= st.span.quantile(.85)
        ax3.scatter(st.b, ys, s=[34 if w else 16 for w in wide],
                    color=[C_ULM if w else C_SEA for w in wide], zorder=3)
        ax3.axvline(-0.5, ls="--", color="k", lw=1.2)
        ax3.set_yticks(ys); ax3.set_yticklabels(st.ct, fontsize=6.5)
        ax3.set_xlabel("Within-cell-type exponent")
        ax3.set_title("c   Not driven by abundance", loc="left",
                      fontsize=11, fontweight="bold")
        ax3.axvline(st.b.median(), color=C_ULM, lw=1.2, alpha=.6)
        ax3.text(.03, .97, f"median {st.b.median():+.3f}\nlarge = widest "
                 f"$n_{{eff}}$ range\n(best determined)", transform=ax3.transAxes,
                 fontsize=7.5, va="top")
    _save(f, out, "Fig3_null_and_stratified")


def _first(res, *names):
    """Cannot be written as _read(a) or _read(b): the truth value of a DataFrame
    raises rather than being falsy when empty."""
    for n in names:
        p = os.path.join(res, n)
        if os.path.exists(p):
            return pd.read_csv(p)
    print(f"  ⚠️ 缺 {' / '.join(names)}，该面板只画拟合线")
    return None


def _read(res, name):
    p = os.path.join(res, name)
    if not os.path.exists(p):
        print(f"  ⚠️ 缺 {name}，该面板只画拟合线")
        return None
    return pd.read_csv(p)


def _save(f, out, name):
    os.makedirs(out, exist_ok=True)
    for ext in ("pdf", "png"):
        f.savefig(os.path.join(out, f"{name}.{ext}"), dpi=300,
                  bbox_inches="tight")
    plt.close(f)
    print(f"  ✅ {out}/{name}.pdf / .png")


# ==================================================================== Fig 1
def fig1(res, out, meta_dir=None):
    """Resource and design. Panel a needs the well-to-chip mapping, which lives in
    the raw metadata rather than in results/."""
    f = plt.figure(figsize=(13, 4.4))
    gs = f.add_gridspec(1, 4, width_ratios=[1.25, 1, 1.45, .85], wspace=.55)
    don = _read(res, "p5_fig1_donor_table.csv")

    # a: donor-by-chip layout, marking the cross-chip replicates
    ax = f.add_subplot(gs[0])
    m = None
    if meta_dir:
        fp = os.path.join(meta_dir, "Multiome_Dataset_Metadata.txt")
        if os.path.exists(fp):
            m = pd.read_csv(fp, sep="\t", usecols=["ID", "Well"], dtype=str)
            m["Chip"] = m.Well.str.extract(r"^(Chip\d+)")
    if m is not None:
        ct = m.groupby(["ID", "Chip"]).size().unstack(fill_value=0)
        nchip = (ct > 0).sum(1)
        order = nchip.sort_values(ascending=False).index          # 重复供体排最上
        ct = ct.loc[order]
        chips = sorted(ct.columns, key=lambda c: int(c.replace("Chip", "")))
        ct = ct[chips]
        ax.imshow(np.log10(ct.values + 1), aspect="auto", cmap="Blues",
                  interpolation="nearest")
        nrep = int((nchip >= 2).sum())
        ax.axhline(nrep - .5, color="#c0504d", lw=1.6)
        ax.text(len(chips) * .5, nrep - 1.6,
                f"{nrep} donors on two chips\n(26 with $\\geq$50 nuclei on both — used)",
                color="#c0504d", fontsize=7.5, ha="center", va="bottom")
        ax.set_xticks(range(len(chips)))
        ax.set_xticklabels([c.replace("Chip", "") for c in chips], fontsize=7)
        ax.set_xlabel("Chip"); ax.set_ylabel("Donor (ordered)")
    else:
        ax.text(.5, .5, "metadata not reachable\n(pass --meta-dir)", ha="center",
                va="center", transform=ax.transAxes, color="#999")
    ax.set_title("a   Donor × chip layout", loc="left", fontsize=11,
                 fontweight="bold")

    # b: quality dimensions, replicated against non-replicated donors
    ax2 = f.add_subplot(gs[1])
    if don is not None and m is not None:
        rep = set((m.groupby("ID").Chip.nunique() >= 2).pipe(lambda s: s[s].index))
        d = don.copy(); d["rep"] = d.ID.isin(rep)
        dims = [("n_nuclei", "Nuclei"), ("median_nCount_RNA", "RNA"),
                ("median_nCount_ATAC", "ATAC")]
        for i, (c, lab) in enumerate(dims):
            if c not in d.columns:
                continue
            for j, (g, cc) in enumerate(((True, C_ULM), (False, "#bbbbbb"))):
                v = pd.to_numeric(d.loc[d.rep == g, c], errors="coerce").dropna()
                v = v / v.median() if len(v) else v          # 各维度归一后同轴
                ax2.scatter(np.full(len(v), i + (j - .5) * .26)
                            + np.random.default_rng(0).normal(0, .035, len(v)),
                            v, s=7, color=cc, alpha=.7, edgecolors="none")
        ax2.set_xticks(range(len(dims)))
        ax2.set_xticklabels([l for _, l in dims], fontsize=8)
        ax2.set_ylabel("Value / median")
        ax2.text(.02, .97, "blue = replicated\ngrey = not\nP = 0.13–0.62",
                 transform=ax2.transAxes, fontsize=7.5, va="top")
    ax2.set_title("b   Quality is unbiased", loc="left", fontsize=11,
                  fontweight="bold")

    # c: the systematic search table, which turns "only one dataset exists"
    # into an argument about scarcity rather than an admission
    ax3 = f.add_subplot(gs[2]); ax3.axis("off")
    rows = [("NIH-CARD PFC", "362 samples", "1 library/donor"),
            ("SEA-AD multiome", "28 libraries", "28 donors"),
            ("ROSMAP multi-region", "283 libraries", "1 per donor×region"),
            ("GSE212630", "14 donors", "below minimum"),
            ("This cohort", "79 donors", "26 replicate pairs")]
    ax3.text(0, 1.02, "c   Survey for replicates", fontsize=11,
             fontweight="bold", transform=ax3.transAxes)
    for i, (a, b, c) in enumerate(rows):
        y = .86 - i * .175
        ok = i == len(rows) - 1
        ax3.text(0, y, ("✓ " if ok else "✗ ") + a, fontsize=8.5,
                 color=("#1f4e79" if ok else "#999"),
                 fontweight="bold" if ok else "normal", transform=ax3.transAxes)
        ax3.text(.60, y, b, fontsize=7.5, color="#666", transform=ax3.transAxes)
        ax3.text(.60, y - .07, c, fontsize=7.5,
                 color=("#1f4e79" if ok else "#999"), transform=ax3.transAxes)

    # d: the three-layer schematic
    ax4 = f.add_subplot(gs[3]); ax4.axis("off")
    ax4.text(0, 1.06, "d   Three layers", fontsize=11, fontweight="bold",
             transform=ax4.transAxes)
    layers = [("L1  Composition", "cell assignment", "#9dc3e6"),
              ("L2  Expression", "aggregated profile", "#5b9bd5"),
              ("L3  Regulation", "per-element, per-cell", "#1f4e79")]
    for i, (t1, t2, c) in enumerate(layers):
        y = .74 - i * .27
        ax4.add_patch(plt.Rectangle((0, y), .95, .19, color=c,
                                    transform=ax4.transAxes, clip_on=False))
        ax4.text(.06, y + .115, t1, fontsize=8.5, fontweight="bold",
                 color="w" if i else "#1f4e79", transform=ax4.transAxes)
        ax4.text(.06, y + .04, t2, fontsize=7,
                 color="w" if i else "#1f4e79", transform=ax4.transAxes)
    ax4.annotate("", xy=(1.02, .70), xytext=(1.02, .12),
                 xycoords=ax4.transAxes, textcoords=ax4.transAxes,
                 arrowprops=dict(arrowstyle="<-", lw=1.4, color="#555"))
    ax4.text(1.06, .40, "increasing\ndemand", fontsize=7.5, rotation=90,
             va="center", transform=ax4.transAxes)

    _save(f, out, "Fig1_resource")


# ==================================================================== Fig 4
def fig4(res, out):
    """The regulatory floor: three lines of evidence plus the power curve."""
    f = plt.figure(figsize=(13, 4.2))
    gs = f.add_gridspec(1, 4, width_ratios=[1, 1, 1.1, 1], wspace=.52)
    fe = _read(res, "p5_claimA_features.csv")
    bg = _read(res, "p5_gene_peak_window_counts.csv")

    # a: the equalised design
    ax = f.add_subplot(gs[0])
    cells = _read(res, "A2_cells_seed0.csv")
    if cells is not None and "celltype" in cells.columns:
        n = cells.groupby("celltype").size()
        ax.bar(range(len(n)), n.values, color=C_ULM, width=.85)
        ax.axhline(150, ls="--", color="k", lw=1)
        ax.text(len(n) * .5, 158, "150 nuclei per subtype", ha="center", fontsize=8)
        ax.set_ylim(0, 200); ax.set_xticks([])
        ax.set_xlabel(f"{len(n)} neuronal subtypes"); ax.set_ylabel("Nuclei")
        ax.text(.03, .06, "depth ratio 1.0000\nno nucleus below target",
                transform=ax.transAxes, fontsize=7.5)
    ax.set_title("a   Equalised design", loc="left", fontsize=11,
                 fontweight="bold")

    # b: independent links per gene -- the primary variable is zero throughout
    ax2 = f.add_subplot(gs[1])
    if fe is not None and "redundancy" in fe.columns:
        # The histogram shows per-type MEANS and the red line the MEDIAN. Two
        # different statistics, so both must be labelled, or the panel looks
        # self-contradictory.
        ax2.hist(fe.redundancy_mean, bins=14, color=C_ULM, alpha=.85,
                 label="mean per subtype")
        ax2.axvline(0, color="#c0504d", lw=2.2, label="median (all subtypes)")
        ax2.legend(fontsize=7.5, frameon=False, loc="upper right")
        ax2.text(.04, .60, "median = 0 in\nall 25 subtypes;\nmeans span only\n"
                 "0.17–0.26", color="#c0504d", transform=ax2.transAxes,
                 fontsize=8, va="top")
    ax2.set_xlabel("Independent links per expressed gene")
    ax2.set_ylabel("Subtypes")
    ax2.set_title("b   No between-subtype variance", loc="left", fontsize=11,
                  fontweight="bold")

    # c: promoter-enrichment odds ratio and the power curve -- the decisive panel
    ax3 = f.add_subplot(gs[2])
    n_link = np.logspace(2, 5, 120)
    p0 = 0.0095                      # 检验集中启动子邻近比例
    for orv, c in ((3.0, "#9dc3e6"), (2.0, "#5b9bd5"), (1.5, "#1f4e79")):
        p1 = orv * p0 / (1 - p0 + orv * p0)
        se = np.sqrt(p1 * (1 - p1) / n_link + p0 * (1 - p0) / 3.03e6)
        ax3.plot(n_link, (p1 - p0) / se, color=c, lw=1.8, label=f"OR = {orv}")
    ax3.axhline(1.96, ls="--", color="k", lw=1)
    ax3.axvline(33227, color="#c0504d", lw=2)
    ax3.text(33227 * 1.15, 1.5, "33,227 links\nobserved", color="#c0504d",
             fontsize=8)
    ax3.set_xscale("log"); ax3.set_xlabel("Number of detected links")
    ax3.set_ylabel("Detection z-score")
    ax3.set_ylim(0, 12); ax3.legend(fontsize=7.5, frameon=False, loc="upper left")
    ax3.set_title("c   Power was ample", loc="left", fontsize=11,
                  fontweight="bold")

    # d: the measured odds ratio, significantly below 1
    ax4 = f.add_subplot(gs[3])
    ax4.barh([0], [0.663], xerr=[[0.663 - 0.568], [0.768 - 0.663]],
             color="#c0504d", height=.4, error_kw=dict(lw=1.5))
    ax4.axvline(1, ls="--", color="k", lw=1.2)
    ax4.set_xlim(0, 1.35); ax4.set_yticks([]); ax4.set_ylim(-.8, .8)
    ax4.set_xlabel("Promoter-enrichment odds ratio")
    ax4.text(.68, .30, "0.663\n[0.568, 0.768]", fontsize=9, color="#c0504d",
             transform=ax4.transAxes)
    ax4.text(.04, .78, "significant depletion,\nnot absence of signal",
             transform=ax4.transAxes, fontsize=7.5)
    ax4.set_title("d   Observed: depletion", loc="left", fontsize=11,
                  fontweight="bold")

    _save(f, out, "Fig4_regulatory_floor")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    ap.add_argument("--out", default="figures")
    ap.add_argument("--meta-dir", default="/sessions/wonderful-relaxed-clarke/mnt/ResearchD",
                    help="Multiome_Dataset_Metadata.txt 所在目录（Fig 1a 需要 Well→Chip）")
    args = ap.parse_args()
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "pdf.fonttype": 42})
    print("Fig 2 …"); fig2(args.results, args.out)
    print("Fig 1 …"); fig1(args.results, args.out, args.meta_dir)
    print("Fig 3 …"); fig3(args.results, args.out)
    print("Fig 4 …"); fig4(args.results, args.out)
    print("Fig 5 …"); fig5(args.results, args.out)
    print("\n⚠️ Ulm 侧的点云文件名若与脚本假设不符（p5_floor_scaling_pairs.csv / "
          "p5_two_layer_L1_pairs.csv），改 _read 的参数即可；\n"
          "   拟合线与统计量来自 decisions/ 中已锁定的值，不受影响。")


if __name__ == "__main__":
    main()


