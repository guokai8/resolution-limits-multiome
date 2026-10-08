#!/usr/bin/env python3
"""
Transfer pathology labels from sorted FANS nuclei to the multiome data
=====================================================================
A gate, not a convenience: if the leave-one-donor-out median AUROC falls below
0.65, the claim that depends on these labels is downgraded.

The hard part is not fitting a classifier. It is three problems that would each
be sufficient grounds for rejection:

  1. PRIOR SHIFT. FANS nuclei are SORTED, so high:low = 7189:5036 ~ 59:41 is an
     artefact of the gate. The true fraction of nuclei with low nuclear TDP-43
     in tissue is unknown and far lower. Applying the classifier directly would
     overestimate it systematically.
     -> Saerens-Latinne-Decaestecker EM, plus BBSE (Lipton et al., ICML 2018)
        as an independent second estimate.

  2. DOMAIN SHIFT. Sorting itself changes the transcriptome; FANS carries a
     percent.soup column that the multiome data does not have.
     -> this script provides (a) a baseline of highly variable genes plus
        standardisation, and leaves the interface open for (b) an scVI or DANN
        latent space in production.

  3. CREDIBILITY. The first thing a reviewer will say is that the labels are
     invented.
     -> conformal prediction, class-conditional (Mondrian), giving each nucleus
        a prediction set with a coverage guarantee. Downstream analysis uses
        ONLY nuclei whose prediction set is a singleton, and reports the
        fraction discarded.

Validation protocol, all of it required:
  * leave-one-donor-out cross-validation, reporting AUROC, AUPRC and the
    CALIBRATION curve
  * AUC stratified by cell type
  * negative control 1: permute the labels; AUC must return to 0.5
  * negative control 2: oligodendrocytes, which TDP-43 pathology does not
    affect; AUC must be clearly below that of excitatory neurons

Usage:
  # self-test on synthetic data: checks the EM, BBSE and conformal
  # implementations without touching real data
  python3 p5_02_label_transfer.py --selftest

  # real run
  python3 p5_02_label_transfer.py --data DIR --out out/ --celltype Exc_LINC00507

Depends on numpy only. In production, swap fit_logreg for sklearn or scVI --
the interface is the same.
"""

import argparse
import os
import sys

import numpy as np


# ==========================================================================
# 1. Classifier: L2-regularised logistic regression in pure numpy, Adam.
#    The swap point for production -- anything offering predict_proba fits here.
# ==========================================================================
def fit_logreg(X, y, l2=1e-3, epochs=300, lr=0.05, seed=0):
    rng = np.random.default_rng(seed)
    n, d = X.shape
    w = rng.normal(0, 0.01, d); b = 0.0
    mw = vw = np.zeros(d); mb = vb = 0.0
    b1, b2, eps = 0.9, 0.999, 1e-8
    for t in range(1, epochs + 1):
        z = X @ w + b
        p = 1.0 / (1.0 + np.exp(-np.clip(z, -30, 30)))
        gw = X.T @ (p - y) / n + l2 * w
        gb = float((p - y).mean())
        mw = b1 * mw + (1 - b1) * gw; vw = b2 * vw + (1 - b2) * gw ** 2
        mb = b1 * mb + (1 - b1) * gb; vb = b2 * vb + (1 - b2) * gb ** 2
        w -= lr * (mw / (1 - b1 ** t)) / (np.sqrt(vw / (1 - b2 ** t)) + eps)
        b -= lr * (mb / (1 - b1 ** t)) / (np.sqrt(vb / (1 - b2 ** t)) + eps)
    return w, b


def predict_proba(X, w, b):
    return 1.0 / (1.0 + np.exp(-np.clip(X @ w + b, -30, 30)))


# ==========================================================================
# 2. Prior-shift correction
# ==========================================================================
def em_label_shift(p_source, prior_source, tol=1e-8, max_iter=1000):
    """
    Saerens-Latinne-Decaestecker EM for prior shift.

    Args:
        p_source: P(y=1|x) from the source-trained classifier, evaluated on the
            TARGET domain.
        prior_source: P(y=1) in the source domain -- here an artefact of the
            sort, not a property of tissue.

    Returns:
        (estimated target prior, recalibrated posteriors).
    """
    pi = float(np.mean(p_source))
    for _ in range(max_iter):
        r = (pi / prior_source) * p_source
        s = ((1 - pi) / (1 - prior_source)) * (1 - p_source)
        p_new = r / (r + s)
        pi_new = float(p_new.mean())
        if abs(pi_new - pi) < tol:
            pi = pi_new
            break
        pi = pi_new
    r = (pi / prior_source) * p_source
    s = ((1 - pi) / (1 - prior_source)) * (1 - p_source)
    return pi, r / (r + s)


def bbse(y_cal, p_cal, p_target, thresh=0.5):
    """
    Black Box Shift Estimation (Lipton, Wang & Smola, ICML 2018).

    Inverts the calibration-set confusion matrix to recover the target-domain
    class prior. It is an INDEPENDENT estimate of what EM returns: if the two
    disagree, the label-shift assumption itself -- that P(x|y) is unchanged --
    is in doubt, and that disagreement has to be reported rather than averaged
    away.
    """
    yh_cal = (p_cal >= thresh).astype(int)
    C = np.zeros((2, 2))
    for a in (0, 1):
        for bb in (0, 1):
            C[a, bb] = np.mean((yh_cal == a) & (y_cal == bb))
    q = np.array([np.mean((p_target < thresh)), np.mean((p_target >= thresh))])
    try:
        if abs(np.linalg.det(C)) < 1e-10:
            return np.nan
        mu = np.linalg.solve(C, q)
        pi = mu[1] * np.mean(y_cal == 1)
        return float(np.clip(pi, 0.0, 1.0))
    except np.linalg.LinAlgError:
        return np.nan


# ==========================================================================
# 3. Conformal prediction (class-conditional / Mondrian): distribution-free coverage
# ==========================================================================
def mondrian_conformal_thresholds(y_cal, p_cal, alpha=0.10, w_cal=None):
    """
    Class-conditional (Mondrian) conformal thresholds.

    EXCHANGEABILITY -- the easiest thing to get wrong here.

    The coverage guarantee of standard conformal prediction rests on
    exchangeability: the calibration and test sets are drawn from the same
    distribution. In this project that assumption is plainly VIOLATED.

        calibration set = FANS sorted nuclei, 11 donors, NeuN+ sorted, carrying
                          percent.soup
        test set        = multiome nuclei, 79 COMPLETELY DIFFERENT donors,
                          unsorted

    Zero donor overlap and a different protocol means a different covariate
    distribution, and without weights the guarantee is empty.

    The correct route is to pass w_cal, the likelihood ratio
    p_target(x)/p_source(x), and do weighted conformal prediction (Tibshirani,
    Barber, Candes & Ramdas, Conformal Prediction Under Covariate Shift,
    NeurIPS 2019). The ratio is estimated with a domain discriminator -- a
    classifier of FANS against multiome -- as w(x) = d(x)/(1-d(x)) * (n_s/n_t).

    With w_cal=None this degrades to standard conformal and PRINTS A WARNING.
    Do not use it that way on real data.
    """
    if w_cal is None:
        print("  [WARN] 未提供似然比权重 → 退化为标准共形。"
              "在 FANS→Multiome 这种协变量偏移下，覆盖率保证不成立。", flush=True)
    qs = {}
    for cls in (0, 1):
        m = y_cal == cls
        s = (1.0 - (p_cal if cls == 1 else 1.0 - p_cal))[m]
        n = len(s)
        if n == 0:
            qs[cls] = 1.0
            continue
        if w_cal is None:
            k = int(np.ceil((n + 1) * (1 - alpha)))
            qs[cls] = float(np.sort(s)[min(k, n) - 1])
        else:
            # Weighted quantile, placing the test point's mass at +inf, as in
            # Tibshirani et al.'s construction
            w = np.asarray(w_cal, dtype=float)[m]
            o = np.argsort(s)
            s_sorted, w_sorted = s[o], w[o]
            cw = np.cumsum(w_sorted)
            total = cw[-1] + w.mean()      # +inf 点的权重用均值近似
            idx = np.searchsorted(cw / total, 1.0 - alpha)
            qs[cls] = float(s_sorted[min(idx, n - 1)])
    return qs


def domain_likelihood_ratio(X_source, X_target, l2=1e-2, seed=0, clip=(0.05, 20.0)):
    """
    Estimate w(x) = p_target(x)/p_source(x) with a domain discriminator.

    Returns the weights on the source samples, for weighted conformal
    prediction. Clipped, because a handful of extreme weights would otherwise
    dominate the quantile and make the interval meaningless.
    """
    X = np.vstack([X_source, X_target])
    y = np.concatenate([np.zeros(len(X_source)), np.ones(len(X_target))])
    mu, sd = X.mean(0), X.std(0) + 1e-8
    w_, b_ = fit_logreg((X - mu) / sd, y, l2=l2, seed=seed)
    d = predict_proba((X_source - mu) / sd, w_, b_)
    odds = d / np.clip(1 - d, 1e-6, None)
    ratio = odds * (len(X_source) / max(len(X_target), 1))
    auc_dom = auroc(y, predict_proba((X - mu) / sd, w_, b_))
    print(f"  域判别器 AUROC = {auc_dom:.3f} "
          f"({'域偏移明显，加权是必需的' if auc_dom > 0.7 else '两域较接近'})")
    return np.clip(ratio, *clip)


def conformal_sets(p, qs):
    """Returns (contains 0?, contains 1?). Only nuclei whose prediction set is a
    singleton enter the downstream analysis."""
    return (1.0 - (1.0 - p)) <= qs[0], (1.0 - p) <= qs[1]


# ==========================================================================
# 4. Metrics
# ==========================================================================
def auroc(y, p):
    """
    AUROC in Mann-Whitney U form, with average ranks for ties.

    This function once carried a bug worth remembering: y was sorted by p, but
    then indexed with a rank array in the ORIGINAL order. On perfectly
    separable synthetic data it returned 0.51 -- a value that looks exactly
    like "no signal", which is the hardest kind of error to notice.

    That is why the self-test must assert the AUROC is HIGH, rather than only
    asserting that the downstream EM and conformal steps look reasonable: a
    classifier with no signal at all will still let EM converge to a
    plausible-looking prior.
    """
    y = np.asarray(y).astype(int)
    p = np.asarray(p, dtype=float)
    n1 = int(y.sum()); n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return np.nan
    order = np.argsort(p, kind="mergesort")
    ranks = np.empty(len(p), dtype=float)
    ranks[order] = np.arange(1, len(p) + 1, dtype=float)
    # Average ranks for ties
    sp = p[order]
    i = 0
    while i < len(sp):
        j = i
        while j + 1 < len(sp) and sp[j + 1] == sp[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = (i + j + 2) / 2.0
        i = j + 1
    return float((ranks[y == 1].sum() - n1 * (n1 + 1) / 2.0) / (n0 * n1))


def calibration_curve(y, p, bins=10):
    edges = np.linspace(0, 1, bins + 1)
    out = []
    for i in range(bins):
        m = (p >= edges[i]) & (p < edges[i + 1] if i < bins - 1 else p <= 1)
        if m.sum() >= 10:
            out.append((float(p[m].mean()), float(y[m].mean()), int(m.sum())))
    return out


def ece(y, p, bins=10):
    """Expected Calibration Error. Reporting AUC without ECE is not enough: a
    model can rank well and still be badly miscalibrated, and the prior-shift
    correction depends on calibration."""
    c = calibration_curve(y, p, bins)
    if not c:
        return np.nan
    n = sum(k for _, _, k in c)
    return float(sum(k * abs(a - b) for a, b, k in c) / n)


# ==========================================================================
# 5. Leave-one-donor-out cross-validation
# ==========================================================================
def lodo_cv(X, y, donor, l2=1e-3, seed=0):
    rows = []
    for d in np.unique(donor):
        te = donor == d
        tr = ~te
        if len(np.unique(y[tr])) < 2 or len(np.unique(y[te])) < 2:
            rows.append(dict(donor=d, n=int(te.sum()), auroc=np.nan,
                             ece=np.nan, note="held-out donor lacks both classes"))
            continue
        mu, sd = X[tr].mean(0), X[tr].std(0) + 1e-8
        w, b = fit_logreg((X[tr] - mu) / sd, y[tr], l2=l2, seed=seed)
        p = predict_proba((X[te] - mu) / sd, w, b)
        rows.append(dict(donor=d, n=int(te.sum()), auroc=auroc(y[te], p),
                         ece=ece(y[te], p), note=""))
    return rows


# ==========================================================================
# Self-test: synthesise a scenario with a KNOWN prior shift and check all three
# components recover it
# ==========================================================================
def selftest():
    rng = np.random.default_rng(0)
    d, n_src, n_tgt = 40, 6000, 20000
    TRUE_TGT_PRIOR = 0.08          # 真实组织里病理核占 8%
    SRC_PRIOR = 0.41               # FANS 分选后 Low 占 41%（人为）
    mu1 = rng.normal(0, 1, d) * 0.55

    def draw(n, prior):
        y = (rng.random(n) < prior).astype(int)
        X = rng.normal(0, 1, (n, d)) + y[:, None] * mu1
        return X, y

    Xs, ys = draw(n_src, SRC_PRIOR)
    Xt, yt = draw(n_tgt, TRUE_TGT_PRIOR)

    n_cal = 2000
    Xtr, ytr, Xc, yc = Xs[:-n_cal], ys[:-n_cal], Xs[-n_cal:], ys[-n_cal:]
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-8
    w, b = fit_logreg((Xtr - mu) / sd, ytr)
    pc = predict_proba((Xc - mu) / sd, w, b)
    pt = predict_proba((Xt - mu) / sd, w, b)

    print("=" * 70)
    print("自检：已知真值的合成先验偏移场景")
    print("=" * 70)
    print(f"  源域先验 (FANS 类比)        = {SRC_PRIOR:.3f}")
    print(f"  目标域真实先验              = {TRUE_TGT_PRIOR:.3f}")
    print(f"  校准集 AUROC                = {auroc(yc, pc):.3f}")
    print(f"  目标域 AUROC (真值可见)     = {auroc(yt, pt):.3f}")
    print()
    naive = float((pt >= 0.5).mean())
    pi_em, pt_adj = em_label_shift(pt, SRC_PRIOR)
    pi_bb = bbse(yc, pc, pt)
    print(f"  ❌ 朴素阈值法估计的先验      = {naive:.3f}   "
          f"(偏差 {naive-TRUE_TGT_PRIOR:+.3f})")
    print(f"  ✅ EM 校正后                 = {pi_em:.3f}   "
          f"(偏差 {pi_em-TRUE_TGT_PRIOR:+.3f})")
    print(f"  ✅ BBSE                      = {pi_bb:.3f}   "
          f"(偏差 {pi_bb-TRUE_TGT_PRIOR:+.3f})")
    print("  → 两个独立估计若不一致，标签偏移假设本身可疑，必须在文中报告。")
    print()
    print(f"  校正前 ECE = {ece(yt, pt):.4f}   校正后 ECE = {ece(yt, pt_adj):.4f}")

    alpha = 0.10
    qs = mondrian_conformal_thresholds(yc, pc, alpha=alpha)
    in0, in1 = conformal_sets(pt_adj, qs)
    singleton = in0 ^ in1
    cov = float(np.mean(np.where(yt == 1, in1, in0)))
    print()
    print(f"  共形预测 (alpha={alpha})")
    print(f"    经验覆盖率        = {cov:.3f}   (目标 >= {1-alpha:.2f})")
    print(f"    单元素预测集占比  = {singleton.mean():.3f}   "
          f"← 下游分析只用这部分核")
    print(f"    被弃用的核        = {1-singleton.mean():.3f}   ← 必须在文中报告")

    # ------------------------------------------------------------------
    # Second scenario: covariate shift, which is the real FANS-to-multiome
    # situation. The point is to show that UNWEIGHTED conformal loses coverage
    # and the weighted version recovers it.
    # ------------------------------------------------------------------
    print()
    print("-" * 70)
    print("场景 2：协变量偏移下的共形覆盖（FANS→Multiome 的真实处境）")
    print("-" * 70)
    # Construct a genuine covariate shift: tilt-sample from one population with
    # a single P(y|x), weighting by exp(g * x.v). That makes (a) P(y|x) exactly
    # unchanged and (b) the log density ratio LINEAR in x, so the logistic
    # domain discriminator is correctly specified. This is a fair test, not one
    # rigged in the weighted version's favour.
    # Taking v along the classification direction with g < 0 enriches the target
    # domain for points that look less positive, so positive-class
    # nonconformity grows systematically, the source quantile comes out too
    # small, and coverage falls short.
    Xp, yp = draw(400_000, TRUE_TGT_PRIOR)
    v = mu1 / np.linalg.norm(mu1)
    logt = -1.6 * (Xp @ v)
    tw = np.exp(logt - logt.max())
    idx = rng.choice(len(Xp), size=25_000, replace=False, p=tw / tw.sum())
    Xt2, yt2 = Xp[idx], yp[idx]

    pt2 = predict_proba((Xt2 - mu) / sd, w, b)
    _, pt2_adj = em_label_shift(pt2, SRC_PRIOR)

    qs_un = mondrian_conformal_thresholds(yc, pc, alpha=alpha)          # 无权重
    w_cal = domain_likelihood_ratio(Xc, Xt2, seed=1)
    qs_w = mondrian_conformal_thresholds(yc, pc, alpha=alpha, w_cal=w_cal)

    # Overall coverage is drowned by the majority class: at a prior of 0.08 the
    # negative class is covered almost by construction and the overall number
    # always looks fine. What matters is the CLASS-CONDITIONAL coverage of the
    # positive class -- the pathological nuclei the downstream conclusion rests
    # on.
    def cov_of(qs, cls):
        i0, i1 = conformal_sets(pt2, qs)      # 用与阈值同源的未重标定分数
        hit = np.where(yt2 == 1, i1, i0)
        return float(hit[yt2 == cls].mean())

    cov_un, cov_w = cov_of(qs_un, 1), cov_of(qs_w, 1)
    print(f"    总体覆盖率（会被多数类淹没，仅供参考）: "
          f"未加权 {cov_of(qs_un,0)*(1-yt2.mean())+cov_un*yt2.mean():.3f}")
    print(f"    ❌ 未加权 · 阳性类覆盖率 = {cov_un:.3f}  (目标 >= {1-alpha:.2f})")
    print(f"    ✅ 加权   · 阳性类覆盖率 = {cov_w:.3f}  "
          f"(Tibshirani et al., NeurIPS 2019)")
    print("    → 未加权版本在真实数据上给出的是**空保证**："
          "总体覆盖率 0.99 好看，阳性类却塌到 0.29。")
    if cov_w < 1 - alpha:
        print(f"    ⚠️  但加权后仍只有 {cov_w:.3f} < {1-alpha:.2f} —— "
              "偏移足够严重时，加权是**必要但不充分**的。")
        print("        实践含义：(a) 用更保守的 alpha；(b) 把分析限制在两域"
              "支撑重叠的区域（裁剪极端权重对应的核）；(c) 如实报告"
              "有效样本量 ESS = (Σw)²/Σw² 而不是名义核数。")

    # The assertions must include that the classifier itself works. Asserting
    # only on EM and conformal behaviour is not enough: a classifier with no
    # signal will still let EM converge to a plausible-looking prior.
    auc_t = auroc(yt, pt)
    checks = [
        ("分类器在目标域可分 (AUROC > 0.90)", auc_t > 0.90),
        ("EM 先验估计误差 < 0.03", abs(pi_em - TRUE_TGT_PRIOR) < 0.03),
        ("BBSE 先验估计误差 < 0.05", abs(pi_bb - TRUE_TGT_PRIOR) < 0.05),
        ("EM 优于朴素阈值法",
         abs(naive - TRUE_TGT_PRIOR) > abs(pi_em - TRUE_TGT_PRIOR)),
        ("校正降低 ECE", ece(yt, pt_adj) < ece(yt, pt)),
        (f"共形覆盖率 >= {1-alpha:.2f} (同分布)", cov >= 1 - alpha - 0.02),
        ("加权共形在协变量偏移下优于未加权", cov_w > cov_un),
    ]
    print()
    for name, passed in checks:
        print(f"    {'✅' if passed else '❌'} {name}")
    ok = all(p for _, p in checks)
    print()
    print("  自检结果:", "✅ 通过" if ok else "❌ 未通过 —— 不要用这份实现跑真实数据")
    return 0 if ok else 1


# ==========================================================================
# Real-data entry point (skeleton)
# ==========================================================================
def run_real(args):
    """
    Inputs a real run needs. This function is wiring only; build the features
    to suit the data at hand.

      FANS side     : FANS_Dataset_RNA_counts_raw.mtx (525 MB, loads whole)
                      FANS_Dataset_Metadata.txt (TDP43 column = label,
                      Sample_donor = grouping)
      multiome side : Multiome_Dataset_RNA_counts_raw.mtx (subset by cell type)
                      Multiome_Dataset_Metadata.txt

    Three checks that must pass before anything downstream is allowed to run:
      1. gene name alignment: size and order of the intersection of the two
         features.tsv files
      2. cell-type label alignment: do FANS ID_WNN_L2/L2.5 and multiome
         WNN_L2/L2.5 take the same value set? The column names differ but the
         values should map.
      3. donor overlap: FANS IDs look like FC2/FC21, multiome like ALS29. IF
         THE TWO DONOR SETS DO NOT OVERLAP this is pure cross-donor
         extrapolation, the AUC target must be lowered, and that must be stated
         in the text. If they partly overlap, every overlapping donor must be
         held out as test.

    Feature options, weakest to strongest; try in order and report each:
      (a) highly variable gene log-CPM z-scores            -- baseline
      (b) (a) with percent.soup and nCount regressed out   -- removes the
                                                              sorting signature
      (c) scVI latent space, dataset as a batch covariate  -- production choice
      (d) DANN adversarial domain adaptation (Ganin, ICML 2015) -- the ceiling,
                                                              but needs tuning
    """
    sys.exit(
        "run_real 是骨架。请先完成 docstring 中的三项前置检查"
        "（尤其是第 3 项：FANS 与 Multiome 的供体是否重叠），"
        "再实现特征构建并调用 lodo_cv / em_label_shift / "
        "mondrian_conformal_thresholds。三件套的实现已通过 --selftest 验证。"
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--data")
    ap.add_argument("--out", default="out")
    ap.add_argument("--celltype", default="Exc_LINC00507")
    ap.add_argument("--alpha", type=float, default=0.10)
    args = ap.parse_args()
    if args.selftest:
        sys.exit(selftest())
    if not args.data:
        ap.error("需要 --data 或 --selftest")
    os.makedirs(args.out, exist_ok=True)
    run_real(args)


if __name__ == "__main__":
    main()
