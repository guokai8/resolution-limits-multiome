# C1 前置文献核查：聚合式 linker 与核数依赖

核查日期 2026-10-03。在声明 C1 的任何新颖性之前完成。
下列四条 DOI 已于同日用 Crossref 解析核对，标题、刊名、卷、年份、作者数均一致：
- Leblanc & Lettre, *Sci Rep* **13** (2023), 2 位作者 — doi:10.1038/s41598-023-31040-w
- Poirion et al., *Genome Biol* **25** (2024), 10 位作者 — doi:10.1186/s13059-024-03374-9
- Zhang et al., *PLoS Comput Biol* **21** (2025), 2 位作者 — doi:10.1371/journal.pcbi.1013697
- Liu et al., *Nat Commun* **16** (2025), 2 位作者 — doi:10.1038/s41467-025-63626-5

结论先说：**C1 的核心机制已有文献，不能当新发现报。** 可以报的是在固定操作点上
的受控对照，以及链接集合层面的位置信息诊断量。

## 已发表的、与 C1 重叠的结论

### 1. 聚合会损害 peak–gene linkage 的准确性 —— 已发表

Enhlink（Poirion et al., *Genome Biol* **25**, 2024；doi:10.1186/s13059-024-03374-9）直接写道，改动
ArchR 的嵌入步骤后发现 **"the first embedding step of ArchR was actually
detrimental to its accuracy"**，并且减少近邻数会提高 ArchR 的准确性。他们还报告
基于模拟的准确性"mostly dependent on average promoter accessibility across cells
and number of cells"，f1 > 0.8 需要 **> 5000 细胞**。

→ 「聚合有害」与「准确性依赖细胞数」两点都不是新的。
→ 其 >5000 细胞阈值与本文的 crossing point（原队列 ≈562、外部队列 ≈2,731 核）
  方向一致，是**外部支持**，应当引用。

### 2. trans 染色体零假设 —— 是 Signac 自己的惯例，不是本文发明

Signac 的 `LinkPeaks` 构造零分布的方式就是取**排除被检基因所在染色体**的、
GC 含量与覆盖度匹配的 peak（trans-links），再把 Pearson r 标准化成 z 分数。

→ C1 用的 trans 配对是同一逻辑，只是作用在**整个链接集合**而不是单个链接上。
  必须写成"沿用 Signac 已发表的 trans 背景逻辑"，不得写成新零假设。
  好处是：这是本领域自己接受的零假设，对照因此更有说服力。

### 3. 细胞数 / 丰度影响 peak–gene 统计 —— 已发表，且是最近的先例

Leblanc, F. J. A. & Lettre, G. "Major cell-types in multiomic single-nucleus
datasets impact statistical modeling of links between regulatory sequences and
target genes." *Scientific Reports* **13**, 3924 (2023).
doi:10.1038/s41598-023-31040-w

其结论：最丰富的细胞类型反而**损失**检出效力，因为 trans 零分布在主导细胞类型
里变成双峰，把 z 分数压低。把 mononuclear phagocytes 从 3,788 降采样到 500 细胞
反而显著提高了 z 分数（CD14 单核细胞 t 检验 P = 1.7 × 10⁻⁸¹）。全数据集
15,113 条链接，只用 CD14 细胞则 2,499 条。并建议至少要有 15 个细胞在两个模态
上都非零。

→ 这是本文 Layer 3 最近的先例，**目前正文未引用，必须补上**。
→ 与本文的关系：方向相反而不冲突。他们说的是 Signac z 分数法在高丰度细胞类型
  里的**校准**问题；本文说的是在一般队列能给出的核数下，启动子富集相对**基因组
  背景**不存在。两者都指向"核数与丰度决定 Layer 3 能不能读"，应当在 Discussion
  里并列。
→ 同时它也是对本文 κ 随丰度上升（r = 0.874）这一观察的独立旁证：丰度在两个层
  面都不是中性的。

### 4. 元细胞构造本身会抬高假阳性相关 —— 已发表

"Library size-stabilized metacells construction enhances co-expression network
analysis in single-cell data"（*PLOS Computational Biology*；
doi:10.1371/journal.pcbi.1013697）报告未受控的 library size 方差会
"inflates false-positive correlations and distorts co-expression networks"。
mcRigor（*Nat Commun* 2025；doi:10.1038/s41467-025-63626-5）给元细胞划分的
可靠性提供统计检验。ArchR 社区自己也在讨论用置换导出经验 FDR 取代参数 FDR
（GreenleafLab/ArchR discussion #1149）。

→ 「聚合体之间的参数 FDR 被高估」在工具社区内部已是共识，不能当新发现。

## C1 还能主张什么

1. **固定操作点上的受控对照。** 同一批核、同一降采样深度、同一 ±500 kb 窗口、
   同一 3 kb 启动子定义、同一 BH FDR、被检验的 gene×peak 集合逐一锁定相同——
   两臂之间**唯一**的差别是聚合。已有文献的对比都跨越了多个口径差异。
2. **链接集合层面的位置信息诊断量。** positional OR = OR(观测) / OR(trans)。
   已有工作用 trans 背景校准单个链接；把它用成"整个检出集合里还剩多少位置信息"
   的读数，核查范围内未见。
3. **定量结论的形式。** 聚合臂给出看起来显著的启动子富集（promoter OR 1.10–1.24，
   区间不含 1），而 trans 零假设**逐位复现**同一个数（positional OR 0.99–1.14）；
   单核臂的位置信息随核数单调上升。这个"看似阳性实为零"的具体形态未见报道。

## 对正文的影响

- 必须新增引用：Leblanc & Lettre 2023；Enhlink 2024。
- C1 在正文里的措辞必须是"确认并定位"而不是"首次发现"。
- Enhlink 的 >5,000 细胞阈值应写进 Discussion，作为 crossing point 的外部旁证。

---

## C1 跑完之后的结论（2026-10-03 晚）

210 个配置，0 失败。两个细胞类型 × 三个 linker 臂 × 六档核数 × 两个零假设，
外加 ±250 kb 窗口与聚合体大小 k 的敏感性。表 20。

**可以主张的**

1. 核数依赖在三个臂里都成立，两个细胞类型都成立，±250 kb 与 ±500 kb 都成立。
   这是审稿人最关心的那个问题的答案：**限制在数据里，不在"未聚合的相关"这个
   算法选择里**。这一条此前没人做过——Enhlink 报的是跨方法的准确性差异，不是
   同一批核上固定操作点的核数响应。

2. 过渡点随 linker 移动，就像它随队列移动一样。单核相关在 400–900 之间穿过 1；
   ArchR 在 150 就已在 1 以上。所以"值不迁移"这个模式又多了一个维度。

3. 在等化设计点（150 核）上，**没有哪个臂给出携带超过弱位置信息的链接集合**：
   四个 臂×细胞类型 组合的位置 OR 是 2.200 (1.20–4.04)、1.396 (0.74–2.65)、
   1.033 (0.90–1.19)、1.263 (1.08–1.47)，两个不可与 1 区分，另两个很小，臂之间
   没有一致的高低。到 900 核四个都明确带位置信息（3.058 / 4.311 / 4.928 / 4.139）。

4. 只用 promoter OR 会看错。ArchR 在 150 核给 1.305，看起来"救回了富集"，但
   trans 零假设逐位复现它（位置 OR 1.033）。这说明本文提出的诊断量必须配零假设
   才可读——这一点加强了诊断量那一节，而不是削弱它。

5. 聚合臂若不用 ArchR 自己的取链接阈值，而用本文别处的 FDR 0.05，会返回
   25,814–230,188 条链接、位置 OR 0.982–1.329，即近乎与位置无关的巨大集合。
   **对照必须用工具自己的默认阈值**，否则就是在打稻草人。这是我第一版犯的错，
   已改。

**不能主张的**

- 不能说"聚合有害"是本文发现——Enhlink 已发表（见上）。本文的贡献是在固定
  操作点上定位它，并给出链接集合层面的位置信息读数。
- 不能说 trans 零假设是本文发明——它是 Signac 的背景逻辑（见上）。
- 没测 Signac 本体，也没测 SCENT / scMultiMap 这类直接建模计数的方法。
  正文限制一节已明说。
- k 是自由参数，影响很大（n=2400 时 k=10/25/50/100 给 2/67/541/2,611 条链接，
  promoter OR 19.343/16.898/5.968/3.273）。固定在 ArchR 默认的 25，没有调参，
  正文按敏感性报告而不是按推荐报告。
