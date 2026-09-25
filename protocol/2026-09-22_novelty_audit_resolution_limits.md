# Novelty audit: the Resolution limits manuscript

Date: 2026-09-22. Subject: `MANUSCRIPT.md` —— "Resolution limits of single-nucleus multiome analysis".
两个定向检索：① 启动子邻近富集是否已被用作 peak–gene link 的质控诊断；
② 技术下限 / 噪声地板在单核数据中是否已有已发表的测量。

⚠️ 此前四份 novelty 核查（2026-08-10、09-20、09-21×2）全部针对**另一个项目**
（donors-vs-cells、信度框架、临床基线 gap），从未核查过本手稿的主张。本文档是第一次。

---

## 总判定

**四条主张里，两条 novelty 站得住，一条需要重新定位，一条是资源性的。**
没有发现任何一条被完全 scoop。但发现一篇**必须引用的最接近先行工作**，
而手稿目前**完全没有参考文献列表**。

---

## 检索 ①：启动子邻近富集作为 link 集诊断

### 结论：未发现直接先行工作；但"启动子邻近是混杂"这一认识已存在，方向相反

**已有的做法是把启动子邻近当作要控制的混杂，不是当作诊断信号。**

- *Multiome Peak-Gene Reranking Benchmark* (Zenodo 22032460) 明确处理
  "promoter-proximity confounding"，理由是 "Raw top-N metrics over the unrestricted
  500 kb universe reward promoter collapse: the distance-only control wins at top-50
  with a median distance of 3.5 bp"。其对策是**按距离分箱**以隔离该效应。
  → 同一个观察（启动子邻近对最容易检出），**相反的用途**：它要消除这个效应，
    手稿要把它当作信号强度的读数。

- peak–gene link 方法的验证**普遍依赖外部 ground truth**：
  SCARlink (Nat Genet 2024) 用 promoter capture Hi-C，并报告 fine-mapped eQTL
  富集 11–15×、GWAS 富集 5–12×；SCEG-HiC (NAR 2026) 用 bulk Hi-C 先验；
  Enhlink (Genome Biol 2024) 亦然。
  → 手稿诊断的卖点正是**不需要 ground truth、不需要功效计算、不需要原始数据**。
    这一点在检索到的方法里没有对应物。

- scATAC 质控的公认指标是 TSS enrichment、FRiP、unique fragments、线粒体比例、
  片段长度分布（见 PEAKQC 2025 的综述性陈述）。**TSS enrichment 是针对 read 分布的，
  不是针对 link 集的**——两者同名但不同对象，写作时须明确区分，否则审稿人会误读成
  "你们重新发明了 TSS enrichment"。

### 判定：**novelty 成立，但必须精确定位**

可以主张的：把启动子邻近富集 OR 作为**link 集层面**、**无需 ground truth** 的
可报告诊断，并给出跨公开 link 集的参考区间。⚠️ 订正 2026-09-24：第五个 link 集因无出处已删除，参考区间为四个集合的 2.45–4.01。

⚠️ 不可主张的：不能说"首次注意到启动子邻近与 link 检出有关"——这个观察已存在。
也不能与 scATAC 的 TSS enrichment 混为一谈。

---

## 检索 ②：技术下限的已发表测量

### 最接近先行工作（**必须引用**）

**Zhang et al., "Precision and Accuracy in Quantitative Measurement of Gene Expression
from Single-cell/nucleus RNA Sequencing Data", *Genomics, Proteomics & Bioinformatics*
2025; PMID 40857558; bioRxiv 2024.04.12.589216.**

23 个数据集、3,682,576 个细胞、339 个样本。与本手稿 Layer 2 撞得最近。

已核实其方法（读了已发表版 PMC12603356 的正文）：

| 维度 | Zhang et al. 2025 | 本手稿 |
|---|---|---|
| 重复如何构造 | **随机切分**同一个体同一细胞类型的细胞为三组 | 同一 donor 跨两块芯片的**独立解离 + 上机 + 建库** |
| 报告的量 | 基因水平 CV（脑 0.68 ± 0.24） | 下限 = 基因间 median \|Δlog₂CPM\| |
| 对细胞数的关系 | **仅报告相关系数**（−0.66、−0.78），定性描述 plateau | 拟合幂律指数 −0.507 / −0.505，含 bootstrap CI |
| 与采样零模型比较 | **无**（全文未涉及 multinomial / Poisson 重采样） | 匹配的多项式零模型 −0.482，超额 1.30× |
| 组成层 | 未涉及 | 有 |
| 调控层 | 未涉及 | 有 |

**结构性差异（可写进 cover letter 与 Discussion）**：随机切分同一管细胞
**在原理上只能看到采样噪声**——它无法观测解离与建库之间的差异，因为那一层根本没有被复制。
而本手稿测到的 1.30× 超额，恰恰就是那一层。两者不是"精度不同"，是**测的不是同一个量**。

⚠️ 这个论证站得住，但写作时必须同时承认：该文确实建立了"精度随细胞数变化"这一事实，
本手稿的增量是**量化其函数形式、给出采样零模型基线、并把可迁移的部分与不可迁移的部分分开**。

### 其他相邻工作（已在前次核查中确认，未被 scoop）

- **scPower**：前瞻性、参数化功效框架，手稿已在 Introduction 点名并正确区分。
- **scPOST (2021)**：differential abundance 的功效模拟，建模 "cluster frequency
  variation across samples"，但不从技术重复测量，也不做标度。
- **sccomp (PNAS 2023) / voomCLR (2024) / propeller**：组成差异检验与偏倚校正，
  处理的是 sample-to-sample 变异的建模，不是技术下限的测量。
- **Squair 2021 (Nat Commun) / Zimmerman 2021**：pseudoreplication，同一层级结构的
  不同用途（用混合模型处理 donor 内依赖）。
- **BEARscc**：用 spike-in 构造**模拟**技术重复评估聚类稳健性——概念相邻，
  对象是聚类不是表达/组成下限。

### 判定：**Layer 2 novelty 从"测量精度"下调为"给出采样基线并区分可迁移性"**

---

## 逐条 novelty 强度（核查后）

| 主张 | 核查前估计 | 核查后 | 理由 |
|---|---|---|---|
| ① 调控层失败 + 启动子富集诊断 | 最强 | **最强（维持）** | 无 ground-truth 的 link 集诊断未见先行；阳性对照 2.45–4.01 完整 |
| ② 组成层：量可复现、标度不可复现 | 中等偏强 | **中等偏强（维持）** | 未见有人用真实技术重复测组成下限并做标度 |
| ③ 表达层标度 | 最弱 | **最弱（进一步下调）** | Zhang 2025 已建立精度–细胞数关系；增量限于零模型基线与函数形式 |
| ④ 重复结构稀缺性调查 | 弱（资源性） | **弱（维持）** | 负面调查，支撑"为何无人做过"，本身不是发现 |

**故事的重心应放在 ① 和 ②**，③ 作为必要的定量铺垫而非卖点。
手稿现有写法（作者自己声明指数不是经验发现）与此判定一致，不需要改结论，只需补引用。

---

## 必须做的三件事

1. **补参考文献列表。** 手稿目前**一条正式引用都没有**，正文仅以 "scPower" 等名字提及。
   这是投稿的硬阻断项，与 novelty 无关。
2. **引用并正面讨论 Zhang et al. 2025**，用上表的结构性差异说明增量。
   不引用会被审稿人直接指出，且显得像是不知道。
3. **精确定位启动子富集诊断**：与 scATAC 的 TSS enrichment 明确区分；
   承认启动子邻近作为混杂已被 benchmark 认识，本文的贡献是把它反向用作信号读数
   并给出参考区间。

## 不可写的措辞

按 [[novelty-check-before-claiming]]：本次检索覆盖的是英文文献的公开检索层面，
不等于穷尽。**不要写 "to our knowledge, no one has..."**。
可写 "we are not aware of a reported diagnostic that..."，或直接陈述做了什么、
与已有工作差在哪，让读者自己判断。

---

## 对发表层级的影响

核查结果**不改变**此前判断：

- **Genome Biology / Genome Medicine —— 首选。** ① 与 ② 均未被 scoop，
  "三层不等价"的框架完整，补上引用后可投。
- **Nature Communications —— 边缘偏可试。** ① 的 novelty 经核查后仍然成立，
  这是此前判断中最大的不确定性，现已解除。若再补一次对已发表组成差异结论的
  再评估（用本文的下限重算它们的效应量），有实质机会。
- **Nature Methods —— 达不到。** 与 novelty 无关，是体裁问题：测量/审计论文，
  无新方法、无软件工具。
- **Bioinformatics / NAR GaB —— 保底。**

---

## 补记 2026-09-24：归属已写入正文

审查意见 B5（"pooled vs donor-level 那节把已有工作当成自己的发现"）已处理。
`### What the floor excludes` 一节开头新增一段，明确声明这个诊断不是本文的观察，
并逐一点名：

| 文献 | 贡献 | 正文编号 |
|---|---|---|
| Cao et al. *BMC Bioinformatics* 2019 (scDC) | 自助法估计 subject 层不确定性 | [12] |
| Büttner et al. *Nat Commun* 2021 (scCODA) | Dirichlet-multinomial，构造上容纳过散度 | [13] |
| Phipson et al. *Bioinformatics* 2022 (propeller) | 最接近本文表述：比例相对二项过散 | [14] |
| Mangiola et al. *PNAS* 2023 (sccomp) | 直接建模均值–变异关系 | [15] |
| Squair et al. *Nat Commun* 2021 | 同一错误在表达层的规模量化（约九成研究） | [16] |

五条引文均于 2026-09-24 对照出版记录核实（卷、期、页、DOI）。

正文现在把本文的贡献限定为：这些方法确立了组间变异必须建模，但**都没有测量技术分量
有多大**，因为都没有同一供体的重复文库。本文加的是那个测量，以及由它支持的两个问题
——有多少对比落在下限以下（与用哪个检验无关），以及有多少落在**多项界**以下（任何
建模选择都救不回来）。

⚠️ 正文全文无 "first time / to our knowledge / novel" 一类未加限定的新颖性措辞，
已检查确认。
