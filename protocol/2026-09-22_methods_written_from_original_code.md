# Methods 五节补写完成 — 照原始代码

Date: 2026-09-22
原始工作区：`~/Desktop/P5_VulnerableEpigenome/`（**不在 ResearchD 内**）

## 更正一处重大错误

2026-09-21 我断定「产生手稿 floor 数字的代码不存在」——错误，因为我只搜了 ResearchD。
代码与结果全部在 `~/Desktop/P5_VulnerableEpigenome/`：
`p5_09_floor_scaling.py`、`p5_10_two_layer_scaling.py`、`p5_13_seaad_L2.py`、
`p5_14_L2_null_and_stratified.py`、`p5_15_promoter_OR_external.py`、
`results/p5_floor_scaling_boot.csv`、`results/p5_two_layer_boot.csv`、
`results/p5_tss_gencode_v32.csv`、`refs/gencode.v32.annotation.gtf.gz` 等。

基于该错误前提，我曾改写 Abstract / Results 三节 / Discussion 两节 / Fig. 2-3 图注，
把已验证数字替换为我自己另一套口径的结果。**该改动已全部回退**，
被替换版本保存在 `plan/band_stratified_analysis/`。

## 手稿数字的出处（已逐一验证）

| 手稿值 | 出处 | 复算 |
|---|---|---|
| L1 指数 −0.152 (−0.340, +0.111) | `p5_two_layer_boot.csv` | −0.148 (−0.340, +0.111) ✓ |
| L2 指数 −0.507 (−0.569, −0.441) | `p5_floor_scaling_boot.csv` | −0.507 (−0.569, −0.442) ✓ |
| 过度离散 4.28 / 3.90 | `p5_L1_scaling_pairs.csv`, `seaad_L1_pairs.csv` | 4.268 / 3.899 ✓ |
| L2 逐对比值 1.30 (IQR 1.21–1.43) | `seaad_L2_obs.csv` + `seaad_L2_null.csv` | 1.30 (1.21–1.43) ✓ |

**过度离散的定义是关键**：= RMS 的 z/√(2/n_eff)，其中 z = |f₁−f₂|/√(p(1−p))。
是 RMS 而非中位数或均值——这是我此前反推失败的原因（中位 1.58、均值 2.55、RMS 4.27）。

## 已补写的 Methods 五节

全部照原始代码，不是照 Results 文字反推：

1. **Pseudobulk construction**（原截断于半句）— 照 `p5_01`
2. **Layer 1 — composition floor** — 照 `p5_10`：z 统计量、n_eff 调和平均、RMS 过度离散
3. **Layer 2 — expression floor** — 照 `p5_14`：`median |Δ log2(CPM+1)|`，基因集用**并集**
   `(p1>0)|(p2>0)` 且 ≥200 个基因；零模型为合并谱按实际深度的多项式重抽样
4. **Scaling fits** — 照 `p5_09`：OLS on log-log vs n_eff，2000 次按 donor 聚类 bootstrap
5. **Layer 3 — regulatory inference** — 照 `p5_07`/`p5_08`：等化设计、±500 kb、
   类型内统一 BH FDR<0.05、|r|>0.7 贪心去共线
6. **Enrichment odds ratios** — 照 `p5_15`：±3 kb 启动子、GENCODE v32、
   Haldane–Anscombe 校正、2000 次二项 bootstrap

手稿现已**无空标题**。补充表路径恢复为 `results/...`（在 P5 工作区内真实存在）。

## 我那套分层重算的处置

保存在 `plan/band_stratified_analysis/`，含被替换的手稿版本与完整记录。
它提出的三点仍值得作为**稳健性检验或审稿应对材料**，但都是方法学分歧而非错误：

1. 逐对基因过滤会让指数偏移约 0.22（在 AIDA 上实测）——但原始代码用的是**并集**
   而非我假设的交集，问题的严重程度需要按原始口径重新评估
2. 解析多项式期望未校准（模拟零模型中位 0.83–0.90 而非 1.0）——原始代码用的是
   `√(p(1−p))` 标准化加 RMS 聚合，与我的 Σ(pa−pb)² 口径不同，需重新核
3. 统一定义下两队列过度离散差 4.7 倍——该结论基于我的口径，在原始口径下不成立
   （原始口径给出 4.28 vs 3.90）

⚠️ 结论：这三点**不能直接套用**，要用必须先按原始口径重做。

## 发表层级

**不变：Genome Biology / Genome Medicine 首选**（现在 Methods 完整、数字全部可追溯，
投稿阻断项已清除）；Nature Communications 边缘，需补对已发表结论的再评估；
Nature Methods 达不到（测量/审计论文，无新方法无软件）；
Bioinformatics / NAR GaB 保底。

剩余工作：Fig. 2/3 需按 `p5_16_figures.py` 重绘（图脚本在，硬编码常量与结果文件一致）。

---

## 补记（2026-09-22）：p5_15 重跑，promoter OR 的 CSV 已补齐

`results/p5_promoterOR_*.csv` 原先缺失（Fig 5a 的数值只存在于 `p5_16_figures.py` 的
硬编码常量里）。已重跑 `p5_15_promoter_OR_external.py` 补齐四个 cellranger-arc link 集。

输入按手稿 Methods 口径准备（脚本 `results_methods/linksets/`）：
- **bedpe**：一律从 `*_analysis.tar.gz` 内的 `analysis/feature_linkage/feature_linkage.bedpe`
  提取。⚠️ 盘上的 `human_brain_3k_feature_linkage.bedpe` 只有 111 字节，
  内容是 S3 `AccessDenied` 错误页（当初下载静默失败），**不可用**；tarball 内的是完整的
  847,032 行。
- **peak 全集**：取 `clustering/atac/*/differential_accessibility.csv` 的 Feature ID 列，
  即手稿所述「per-cluster differential accessibility output」，而非 motif 映射文件。
- **TSS**：`results/p5_tss_gencode_v32.csv`（GENCODE v32，与 10x GRCh38-2020-A 匹配）。

结果与手稿逐位一致：

| link set | 检出 link | 检验集 | OR | 95% CI | 手稿 |
|---|---|---|---|---|---|
| pbmc_granulocyte_10k | 476,703 | 5,076,670 | 2.451 | [2.399, 2.506] | ✓ |
| pbmc_granulocyte_3k | 101,670 | 3,621,717 | 3.074 | [2.961, 3.196] | ✓ |
| lymph_node_lymphoma_14k | 134,955 | 3,830,965 | 3.575 | [3.466, 3.693] | ✓ |
| human_brain_3k | 178,195 | 4,325,771 | 4.012 | [3.906, 4.125] | ✓ |

四套窗口均为 ±1 Mb（由 bedpe 实测最大距离定），基因名映射率 99.8–100%。

**仍缺**：第五个 link 集（已发表 AD multiome 补充表）不在盘上，其 OR 5.40 [5.30, 5.51]
⚠️ **订正 2026-09-24**：该集合始终未能找回出处、输入文件或可重算的代码，已从手稿整个删除，
正文改为四个 link 集、区间 2.45–4.01（本就是那四个的区间）。以下关于第五个集合的记述仅为历史记录。
仍为记录值；本研究那条 0.663 [0.568, 0.768] 由 L3 主分析给出。
Fig 5a 现已改为直接读取四个 CSV，这两条在脚本中显式标注来源。

## 图件状态

五张主图全部 ggplot2，共用 `code_figures/fig_common.R` 的主题与配色：

| 脚本 | 产出 |
|---|---|
| `code_figures/fig2_fig3_ggplot2.R` | Fig 2、Fig 3 |
| `code_figures/fig1_fig4_fig5_ggplot2.R` | Fig 1、Fig 4、Fig 5 |

180 mm 双栏宽，矢量 PDF（cairo）+ 450 dpi PNG（ragg）。
Fig 1b 的四个质量维度为总核数、RNA counts、ATAC fragments、**线粒体比例**
（不是年龄），重复组取 27 人，如此方得 P = 0.13–0.62。
