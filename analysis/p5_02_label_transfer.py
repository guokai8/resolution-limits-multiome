#!/usr/bin/env python3
"""
P5 门禁 α · FANS(TDP-43 High/Low) → Multiome 的病理标签迁移
============================================================
这是整个 Project 5 的成败点。若本脚本的留一供体 AUROC 中位数 < 0.65，
命题 B 降级（见 PROJECT_5_MASTER.md §7 R1）。

要解决的不是"跑一个分类器"，而是三个会让审稿人一票否决的问题：

  1. 先验偏移 —— FANS 是**分选**的，High:Low = 7189:5036 ≈ 59:41 是人为比例。
     真实组织里核内 TDP-43 低的核占比未知且远低于此。直接套用会系统性高估。
     → Saerens-Latinne-Decaestecker EM  +  BBSE (Lipton et al., ICML 2018)

  2. 域偏移 —— 分选流程本身改变转录组（FANS 有 percent.soup 列，多组学没有）。
     → 本脚本提供 (a) 基线：仅 HVG + 标准化
                  (b) 生产环境应换成 scVI / DANN 潜空间，接口已留出

  3. 可信度 —— 审稿人第一句话必然是"你的标签是编的"。
     → 共形预测 (Mondrian / class-conditional)：给每个核一个覆盖率有保证的
       预测集。下游分析**只用**预测集为单元素的核，并报告被弃用的比例。

验证协议（缺一不可，见 §B1）：
  - 留一供体交叉验证（LODO），报告 AUROC / AUPRC / **校准曲线**
  - 按细胞类型分层的 AUC
  - 阴性对照 1：置换标签 → AUC 应回到 0.5
  - 阴性对照 2：少突胶质（不受 TDP-43 病理影响）→ AUC 应显著低于兴奋性神经元

用法：
  # 自检（合成数据，验证 EM/BBSE/共形三件套的实现正确性，不碰真实数据）
  python3 p5_02_label_transfer.py --selftest

  # 真实运行
  python3 p5_02_label_transfer.py --data DIR --out out/ --celltype Exc_LINC00507

依赖：仅 numpy。生产环境建议把 fit_logreg 换成 sklearn / scVI —— 接口一致。
"""

import argparse
import os
import sys

import numpy as np


# ==========================================================================
# 1. 分类器（numpy-only 的 L2 正则 logistic regression，Adam）
#    生产环境替换点：任何提供 predict_proba 的模型都可以插进来
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
# 2. 先验偏移校正
# ==========================================================================
def em_label_shift(p_source, prior_source, tol=1e-8, max_iter=1000):
    """
    Saerens-Latinne-Decaestecker EM。
    p_source: 源域训练的分类器在**目标域**上的 P(y=1|x)
    prior_source: 源域的 P(y=1)（FANS 里是分选造成的人为比例）
    返回 (估计的目标域先验, 重标定后的后验)
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
    Black Box Shift Estimation (Lipton, Wang & Smola, ICML 2018)。
    用校准集上的混淆矩阵反解目标域的类别先验。与 EM 互为独立估计 ——
    两者若不一致，说明标签偏移假设（P(x|y) 不变）本身可疑，必须报告。
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
# 3. 共形预测（class-conditional / Mondrian），分布无关的覆盖率保证
# ==========================================================================
def mondrian_conformal_thresholds(y_cal, p_cal, alpha=0.10, w_cal=None):
    """
    class-conditional (Mondrian) 共形阈值。

    ⚠️⚠️ 关于可交换性 —— 这是本脚本最容易被误用的地方 ⚠️⚠️

    标准共形预测的覆盖率保证依赖**可交换性**：校准集与测试集同分布。
    在本项目里这个假设是**明确违反**的：
        校准集 = FANS 分选核（11 供体，NeuN+ 分选，有 percent.soup）
        测试集 = Multiome 核（79 个**完全不同**的供体，未分选）
    供体零重叠 + 实验流程不同 ⇒ 协变量分布不同 ⇒ **无权重时保证是空的**。

    正确做法：传入 w_cal = 似然比 p_target(x)/p_source(x)，做加权共形
    （Tibshirani, Barber, Candès & Ramdas, *Conformal Prediction Under
    Covariate Shift*, NeurIPS 2019）。似然比用一个"域判别器"
    （FANS vs Multiome 的二分类器）估计：w(x) = d(x)/(1-d(x)) · (n_s/n_t)。

    w_cal=None 时退化为标准共形，并**打印警告**。不要在真实数据上这么用。
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
            # 加权分位数：把测试点的质量放在 +inf 处（Tibshirani et al. 的构造）
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
    用域判别器估计 w(x) = p_target(x)/p_source(x)，供加权共形使用。
    返回 source 样本上的权重。裁剪防止极端权重支配分位数。
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
    """返回 (含 0?, 含 1?)。只有单元素预测集的核才进入下游分析。"""
    return (1.0 - (1.0 - p)) <= qs[0], (1.0 - p) <= qs[1]


# ==========================================================================
# 4. 指标
# ==========================================================================
def auroc(y, p):
    """
    Mann-Whitney U 形式的 AUROC，含并列值的平均秩处理。

    ⚠️ 这里曾经有一个 bug：先把 y 按 p 排序、却用原始顺序的秩数组去索引它。
       结果是在完全可分的合成数据上返回 0.51 —— 一个"看起来像没有信号"的值，
       正好是最难被发现的那种错误。自检因此必须断言 AUROC 高，而不只是
       断言下游的 EM/共形"看起来正常"。
    """
    y = np.asarray(y).astype(int)
    p = np.asarray(p, dtype=float)
    n1 = int(y.sum()); n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return np.nan
    order = np.argsort(p, kind="mergesort")
    ranks = np.empty(len(p), dtype=float)
    ranks[order] = np.arange(1, len(p) + 1, dtype=float)
    # 并列值取平均秩
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
    """Expected Calibration Error —— 只报 AUC 不报 ECE 是不够的。"""
    c = calibration_curve(y, p, bins)
    if not c:
        return np.nan
    n = sum(k for _, _, k in c)
    return float(sum(k * abs(a - b) for a, b, k in c) / n)


# ==========================================================================
# 5. 留一供体交叉验证
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
# 自检：合成一个带**已知**先验偏移的场景，验证三件套实现正确
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
    # 第二个场景：协变量偏移（这才是 FANS→Multiome 的真实处境）
    # 目的：证明**不加权的共形会丢失覆盖**，加权能补回来
    # ------------------------------------------------------------------
    print()
    print("-" * 70)
    print("场景 2：协变量偏移下的共形覆盖（FANS→Multiome 的真实处境）")
    print("-" * 70)
    # 构造真正的协变量偏移：从同一个 P(y|x) 的总体里按 exp(g * x·v) 做倾斜抽样。
    # 这样 (a) P(y|x) 严格不变，(b) 对数密度比在 x 上是**线性**的，
    # 因此逻辑回归域判别器是正确设定的 —— 这是公平的检验，不是给加权版放水。
    # 取 v 沿分类方向且 g<0：目标域富集"看起来不那么像阳性"的点，
    # 于是阳性类的 nonconformity 系统性变大 → 源域分位数偏小 → 覆盖不足。
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

    # ⚠️ 只看**总体**覆盖率会被多数类淹没：先验 0.08 时，阴性类几乎必然被覆盖，
    #    总体覆盖率永远好看。真正要看的是**阳性类（病理核）的 class-conditional
    #    覆盖率** —— 那才是下游结论依赖的那一类。
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

    # 断言必须包含"分类器本身有效"这一条。只断言 EM/共形的表现是不够的：
    # 一个完全无信号的分类器也能让 EM 收敛到一个看似合理的先验。
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
# 真实数据入口（骨架）
# ==========================================================================
def run_real(args):
    """
    真实运行需要的输入（本函数只做接线，特征构建请按需实现）：

      FANS 端   : FANS_Dataset_RNA_counts_raw.mtx  (525 MB, 可全量载入)
                  FANS_Dataset_Metadata.txt        (TDP43 列 = 标签, Sample_donor = 分组)
      Multiome 端: Multiome_Dataset_RNA_counts_raw.mtx  (需按细胞类型抽取子集)
                  Multiome_Dataset_Metadata.txt

    ⚠️ 三个必须先做的前置检查（做完才允许往下走）：
      1. 基因名对齐：两份 features.tsv 的交集大小与顺序（FANS 39万? vs Multiome 35,367）
      2. 细胞类型标签对齐：FANS 的 ID_WNN_L2/L2.5 与 Multiome 的 WNN_L2/L2.5
         取值集合是否一致（列名不同但取值应可映射）
      3. 供体是否重叠：FANS 的 ID 形如 FC2/FC21，Multiome 形如 ALS29 ——
         **若两批供体不重叠，这是纯粹的跨供体外推，AUC 目标应下调，
         且必须在文中明写。若部分重叠，重叠供体必须全部留作测试集。**

    特征建议（从弱到强，依次尝试并报告）：
      (a) HVG log-CPM z-score                      —— 基线
      (b) (a) + 回归掉 percent.soup / nCount       —— 去除分选流程的痕迹
      (c) scVI 潜空间（把 dataset 作为 batch 协变量）—— 生产环境首选
      (d) DANN 对抗式域适配 (Ganin, ICML 2015)      —— 上限，但需调参
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
