#!/usr/bin/env Rscript
# Fig 1、Fig 4、Fig 5 重绘（ggplot2），风格与 Fig 2/3 统一
#
# Fig 1  资源与设计
#   a  供体 × 芯片排布，跨芯片重复的 27 人排在上方
#   b  重复 vs 非重复供体的质量维度（各维度除以中位数后同轴）
#   c  公开脑单核资源的系统检索：四套皆不合格
#   d  三层推断的递增需求
#
# Fig 4  L3 下限
#   a  等化设计（25 个亚型各 150 核）
#   b  每基因独立 link 数：主变量在所有亚型恒为 0
#   c  功效曲线：33,227 条 link 足以检出任何合理幅度的富集
#   d  实测 OR 显著 < 1（富集的反面）
#
# Fig 5  诊断的参考区间
#   a  六个 link 集的 OR 森林图
#   b  窗口敏感性
#   c  SEA-AD 28 个 multiome 文库：linkage 未饱和
#
# 用法： Rscript fig1_fig4_fig5_ggplot2.R [results_dir] [out_dir] [meta_dir]

args <- commandArgs(trailingOnly = TRUE)
RES  <- path.expand(if (length(args) >= 1) args[1] else "~/Desktop/P5_VulnerableEpigenome/results")
OUT  <- path.expand(if (length(args) >= 2) args[2] else "~/Desktop/ResearchD/figures")
META <- path.expand(if (length(args) >= 3) args[3] else
  "~/Desktop/ResearchD/Data/Other_Datasets/Multiome_Dataset/files")

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                  value = TRUE)[1])), "fig_common.R"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) read_csv(file.path(RES, f), show_col_types = FALSE)

## ======================================================================
## Fig 1
## ======================================================================
meta <- read.delim(file.path(META, "Multiome_Dataset_Metadata.txt"),
                   colClasses = "character", quote = "", comment.char = "")
meta$Chip <- sub("Well[0-9]+$", "", meta$Well)

cc <- meta %>% count(ID, Chip, name = "n")
nchip <- cc %>% count(ID, name = "n_chip")
ord <- nchip %>%
  left_join(cc %>% group_by(ID) %>% summarise(tot = sum(n), .groups = "drop"), by = "ID") %>%
  arrange(desc(n_chip), desc(tot))
n_rep <- sum(ord$n_chip >= 2)

grid <- expand.grid(ID = ord$ID, Chip = sort(unique(meta$Chip)),
                    stringsAsFactors = FALSE) %>%
  left_join(cc, by = c("ID", "Chip")) %>%
  mutate(n = tidyr::replace_na(n, 0),
         ID = factor(ID, levels = rev(ord$ID)),
         Chip = factor(Chip, levels = paste0("Chip", sort(as.integer(
           sub("Chip", "", unique(meta$Chip)))))))

p1a <- ggplot(grid, aes(Chip, ID, fill = log10(n + 1))) +
  geom_tile() +
  geom_hline(yintercept = nrow(ord) - n_rep + 0.5, colour = C_SEA, linewidth = .6) +
  annotate("text", x = 6.5, y = nrow(ord) - n_rep + 1.6, size = 2.3, colour = C_SEA,
           lineheight = .95, vjust = 0,
           label = paste0("Donor-matched process replicates:\n", n_rep,
                          " donors on two chips (26 with \u226550 nuclei on both)")) +
  scale_fill_gradient(low = "#F2F7FC", high = C_ULM) +
  scale_x_discrete(labels = function(x) sub("Chip", "", x)) +
  labs(x = "Chip", y = "Donor (ordered)") +
  theme_pub() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 6.5))

mito <- meta %>%
  mutate(pmito = suppressWarnings(as.numeric(percent_mito))) %>%
  group_by(ID) %>% summarise(pmito = median(pmito, na.rm = TRUE), .groups = "drop")
don <- rd("p5_fig1_donor_table.csv") %>%
  left_join(nchip, by = "ID") %>% left_join(mito, by = "ID") %>%
  mutate(rep = n_chip >= 2)
# 手稿的四个质量维度：总核数、RNA counts、ATAC fragments、线粒体比例
## 按 GB 要求精简到三个与概念信息直接相关的维度，线粒体比例移出主图
dims <- c(n_nuclei = "Nuclei", median_nCount_RNA = "RNA\ncomplexity",
          median_nCount_ATAC = "ATAC\ncomplexity")
qual <- lapply(names(dims), function(v) {
  x <- suppressWarnings(as.numeric(don[[v]]))
  data.frame(dim = dims[[v]], rep = don$rep, val = x / median(x, na.rm = TRUE))
}) %>% bind_rows() %>% filter(is.finite(val), !is.na(rep)) %>%
  mutate(dim = factor(dim, levels = unname(dims)))

pv <- sapply(names(dims), function(v) {
  x <- suppressWarnings(as.numeric(don[[v]]))
  suppressWarnings(wilcox.test(x[don$rep], x[!don$rep])$p.value)
})

set.seed(0)
p1b <- ggplot(qual, aes(dim, val, colour = rep)) +
  geom_point(position = position_jitterdodge(jitter.width = .22, dodge.width = .55),
             size = .75, alpha = .75, stroke = 0) +
  scale_colour_manual(values = c(`TRUE` = C_ULM, `FALSE` = C_GREY)) +
  annotate("text", x = .55, y = Inf, hjust = 0, vjust = 1.25, size = 2.3,
           colour = "grey25", lineheight = 1,
           label = sprintf("blue = replicated\ngrey = not\nP = %.2f-%.2f",
                           min(pv), max(pv))) +
  # 顶部留白，否则 Nuclei 列最高的点会压在图例文字上
  scale_y_continuous(expand = expansion(mult = c(.05, .30))) +
  labs(x = NULL, y = "Value / median") +
  theme_pub()

## ---- 1c  三个统计单元的示意图（按 GB 要求重做）------------------------------
## 要传达的是：同一批数据支持三种断言，各自的统计单元不同，因此分辨率不同。
box <- data.frame(
  x    = c(1, 2.5, 4),
  lab  = c("Composition", "Expression", "Regulation"),
  unit = c("donor\nproportion", "donor\npseudobulk", "nucleus,\nfeature\ncovariance"),
  col  = c(L_COMP, L_EXPR, L_REG))

p1d <- ggplot(box) +
  ## 顶部来源
  annotate("label", x = 2.5, y = 4.42, label = "Single-nucleus multiome data",
           size = 2.3, family = BASE_FAMILY, fill = "grey95",
           label.size = 0, label.padding = unit(2.6, "pt"), colour = "grey15") +
  ## 分叉
  annotate("segment", x = 2.5, xend = 2.5, y = 4.18, yend = 3.88,
           colour = "grey55", linewidth = .4) +
  annotate("segment", x = 1, xend = 4, y = 3.88, yend = 3.88,
           colour = "grey55", linewidth = .4) +
  annotate("segment", x = box$x, xend = box$x, y = 3.88, yend = 3.52,
           colour = "grey55", linewidth = .4,
           arrow = arrow(length = unit(3, "pt"), type = "closed")) +
  ## 三个断言
  geom_label(aes(x, 3.26, label = lab, fill = I(col)), colour = "white",
             size = 2.15, fontface = "bold", family = BASE_FAMILY,
             label.size = 0, label.padding = unit(2.2, "pt")) +
  ## 到统计单元
  annotate("segment", x = box$x, xend = box$x, y = 3.0, yend = 2.62,
           colour = "grey55", linewidth = .4,
           arrow = arrow(length = unit(3, "pt"), type = "closed")) +
  geom_text(aes(x, 2.22, label = unit, colour = I(col)),
            size = 2.05, family = BASE_FAMILY, lineheight = .95) +
  ## 结论
  annotate("segment", x = .5, xend = 4.5, y = 1.62, yend = 1.62,
           colour = "grey75", linewidth = .35) +
  annotate("text", x = 2.5, y = 1.3, size = 2.35, colour = "grey15",
           family = BASE_FAMILY, lineheight = 1.05,
           label = "Different statistical units\n\u2192 different resolution limits") +
  scale_x_continuous(limits = c(0.3, 4.7)) +
  scale_y_continuous(limits = c(0.95, 4.68)) +
  theme_blank()


fig1 <- p1a + p1b + p1d +
  plot_layout(widths = c(1.1, 0.9, 1.15)) +
  plot_annotation(tag_levels = "a")
save_fig(fig1, OUT, "Fig1_resource", height = 68)

## ======================================================================
## Fig 4
## ======================================================================
cells <- rd("A2_cells_seed0.csv")
fe    <- rd("p5_claimA_features.csv")

nsub <- cells %>% count(celltype, name = "n") %>% arrange(desc(n)) %>%
  mutate(i = row_number())

p4a <- ggplot(nsub, aes(i, n)) +
  geom_col(fill = C_ULM, width = .85) +
  geom_hline(yintercept = 150, linetype = "22", linewidth = .4, colour = "grey20") +
  annotate("text", x = nrow(nsub) / 2, y = 157, size = 2.05, colour = "grey25",
           family = BASE_FAMILY, label = "150 nuclei per subtype") +
  scale_y_continuous(limits = c(0, 228), expand = expansion(mult = c(0, .02))) +
  # QC 文字一律放到柱体之外，避免白字压在柱间空隙上不可读
  ann(x = .5, y = 226, hjust = 0, size = 2.0, colour = "grey30",
      lab = paste("Test: equalised nucleus design",
                  "(reference link sets at full depth)",
                  "depth ratio 1.0000, no nucleus below target", sep = "\n")) +
  labs(y = "Nuclei", x = paste0(nrow(nsub), " neuronal subtypes")) +
  theme_pub() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# b: 每个 peak 在该亚型 150 个核中被检出的数量（中位，逐亚型）
pk <- rd("panel_b_nuclei_per_peak_by_subtype.csv") %>%
  arrange(median_nuclei_per_peak) %>% mutate(i = row_number())
med_pk <- median(pk$median_nuclei_per_peak)

p4b <- ggplot(pk, aes(i, median_nuclei_per_peak)) +
  geom_linerange(aes(ymin = q25, ymax = q75), colour = C_LITE, linewidth = 1.1) +
  geom_point(colour = C_ULM, size = 1.2) +
  geom_hline(yintercept = med_pk, colour = C_SEA, linewidth = .55) +
  annotate("text", x = 1, y = 7.2, hjust = 0, size = 2.1, colour = C_SEA,
           lineheight = 1,
           label = sprintf("median %.0f of 150 nuclei\nin every subtype", med_pk)) +
  annotate("text", x = nrow(pk), y = 1.15, hjust = 1, size = 2.0, colour = "grey35",
           label = sprintf("%.0f%% of peaks seen in one nucleus",
                           100 * median(pk$frac_peaks_in_1_nucleus))) +
  scale_y_continuous(limits = c(0.8, 8.2), breaks = c(1, 2, 4, 6, 8)) +
  ann(x = 1, y = 8.1, lab = "Sparse detection", size = 2.5, colour = "grey15") +
  labs(x = "Neuronal subtypes", y = "Nuclei per detected peak") +
  theme_pub() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# c: 实测 OR —— 显著耗竭（正文引 Fig. 4c）
orv <- data.frame(or = 0.663, lo = 0.568, hi = 0.768)
p4c <- ggplot(orv) +
  geom_errorbar(aes(1, ymin = lo, ymax = hi), width = .06, linewidth = .5,
                colour = C_SEA) +
  geom_point(aes(1, or), colour = C_SEA, size = 2.6) +
  geom_hline(yintercept = 1, linetype = "22", linewidth = .45, colour = "grey20") +
  annotate("text", x = 1.22, y = 0.663, hjust = .5, size = 2.3, colour = C_SEA,
           family = BASE_FAMILY, lineheight = 1,
           label = "0.663\n[0.568, 0.768]") +
  annotate("text", x = .72, y = 0.04, hjust = 0, size = 2.0, colour = "grey30",
           family = BASE_FAMILY, lineheight = 1,
           label = "reduced recovery,\nnot absence of regulation") +
  coord_flip(xlim = c(.6, 1.45), ylim = c(0, 1.3)) +
  labs(x = NULL, y = "Promoter-enrichment odds ratio") +
  theme_pub() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())


## 检出 z 分析按 GB 要求移到补充图（SuppFig6），主图保留三个面板
fig4 <- p4a + p4b + p4c +
  plot_layout(widths = c(1, 1.2, .95)) + plot_annotation(tag_levels = "a")
save_fig(fig4, OUT, "Fig4_regulatory_floor", height = 68)

## ======================================================================
## Fig 5
## ======================================================================
# 四个 cellranger-arc link 集直接读 p5_15 的产出；
# ⚠️ 另两行无对应 CSV：已发表 AD multiome 的补充表不在盘上，
#    本研究那一条由 L3 主分析给出，故按记录的数值补入。
# ⚠️ 必须用与主分析同窗口（±500 kb）的那一版：data/derived_results/window_500000/。
#    旧的 ±1 Mb 版本（且窗口不自洽）已移入 superseded_inconsistent_window/，与 0.663 不可比。
W500 <- file.path(dirname(RES), basename(RES), "window_500000")
if (!dir.exists(W500)) W500 <- file.path("data/derived_results", "window_500000")
ext <- list.files(W500, "^p5_promoterOR_.*\\.csv$", full.names = TRUE) %>%
  lapply(read_csv, show_col_types = FALSE) %>% bind_rows() %>%
  transmute(lab = sub("_sorted", "", label), or = OR, lo, hi)
stopifnot(nrow(ext) == 4)
# 四个 cellranger-arc 参照集直接读 p5_15 的产出；本研究那一条由 L3 主分析给出。
or_sets <- bind_rows(
  ext %>% arrange(or),
  data.frame(lab = "This study (150 nuclei)", or = 0.663, lo = 0.568, hi = 0.768)) %>%
  mutate(lab = factor(lab, levels = rev(lab)), dep = or < 1)

p5a <- ggplot(or_sets, aes(or, lab, colour = dep)) +
  annotate("rect", xmin = 2.26, xmax = 3.56, ymin = -Inf, ymax = Inf,
           fill = C_ULM, alpha = .07) +   # 四个参照 link 集的实际区间
  # 本研究的工作点单独加底色条，按 GB 要求在面板内直接突出
  # 注意：x 为 log10 标度，-Inf/Inf 会变成 NaN 而整块被丢弃，必须给有限边界
  annotate("rect", xmin = .451, xmax = 8.95, ymin = .42, ymax = 1.5,
           fill = C_SEA, alpha = .10) +
  geom_vline(xintercept = 1, linetype = "22", linewidth = .45, colour = "grey20") +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = .9) +
  geom_point(size = 1.5) +
  geom_text(aes(x = hi * 1.13, label = sprintf("%.2f", or)), hjust = 0, size = 2.2) +
  ann(x = .465, y = 1.42, lab = "operating point", size = 2.1, colour = C_SEA,
      hjust = 0, vjust = .5) +
  ann(x = 2.84, y = 5.42, lab = "reference range", size = 2.1, colour = C_ULM,
      hjust = .5, vjust = .5) +
  scale_colour_manual(values = c(`TRUE` = C_SEA, `FALSE` = C_ULM)) +
  scale_x_log10(limits = c(.45, 9), breaks = c(.5, 1, 2, 5),
                labels = c("0.5", "1", "2", "5")) +
  scale_y_discrete(expand = expansion(add = c(.6, .9))) +
  labs(y = NULL,
       x = "Promoter-enrichment odds ratio (log scale)") +
  theme_pub() +
  theme(axis.text.y = element_text(size = 6.5))

# 同一个 link 集 (human_brain_3k) 在三个窗口下的 OR，由 p5_15 --window 重算，
# 背景集与检出集同步收紧（见 code/analysis/p5_19_recompute_refs_at_window.sh）。
win <- data.frame(w = factor(c("\u00b1500 kb", "\u00b11 Mb", "\u00b11.7 Mb"),
                             levels = c("\u00b1500 kb", "\u00b11 Mb", "\u00b11.7 Mb")),
                  or = c(3.558, 4.054, 6.593))
p5b <- ggplot(win, aes(w, or, fill = w)) +
  geom_col(width = .55) +
  geom_text(aes(label = sprintf("%.2f", or)), vjust = -0.5, size = 2.6,
            colour = "grey15") +
  scale_fill_manual(values = c(C_ULM, C_MID, C_LITE)) +
  scale_y_continuous(limits = c(0, 7.6), expand = expansion(mult = c(0, .02))) +
  scale_x_discrete(labels = c("\u00b1500 kb" = "\u00b1500\nkb", "\u00b11 Mb" = "\u00b11\nMb",
                              "\u00b11.7 Mb" = "\u00b11.7\nMb")) +
  labs(x = NULL, y = "Odds ratio, same link set") +
  theme_pub() +
  theme(axis.text.x = element_text(size = 6.8))

lib <- rd("seaad_L3_libraries.csv") %>%
  filter(!is.na(n_cells), !is.na(link), !is.na(atac_frag)) %>%
  mutate(x = n_cells * atac_frag)
b_sat <- 0.421
a_sat <- exp(median(log(lib$link) - b_sat * log(lib$x)))
sat <- data.frame(x = range(lib$x)) %>% mutate(y = a_sat * x^b_sat)

p5c <- ggplot(lib, aes(x, link)) +
  geom_point(colour = C_SEA, size = 1.1, alpha = .85, stroke = 0) +
  geom_line(data = sat, aes(x, y), colour = "grey15", linewidth = .7) +
  annotate("text", x = min(lib$x) * 1.15, y = max(lib$link) * .95, hjust = 0,
           size = 2.3, colour = "grey20",
           label = sprintf("b = +%.2f  (no plateau)", b_sat)) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "bl", size = .25, short = unit(1, "pt"),
                      mid = unit(1.6, "pt"), long = unit(2.4, "pt")) +
  labs(
       x = "Nuclei \u00d7 fragments per nucleus", y = "Feature linkages detected") +
  theme_pub()

fig5 <- p5a + p5b + p5c +
  plot_layout(widths = c(1.25, .7, 1)) + plot_annotation(tag_levels = "a")
save_fig(fig5, OUT, "Fig5_diagnostic", height = 62)

cat(sprintf("Fig1: %d donors, %d chips, %d replicated; P = %.2f-%.2f\n",
            nrow(ord), n_distinct(meta$Chip), n_rep, min(pv), max(pv)))
cat(sprintf("Fig4: %d subtypes x %d nuclei; redundancy median %.1f, means %.3f-%.3f\n",
            nrow(nsub), unique(nsub$n)[1], median(fe$redundancy),
            min(fe$redundancy_mean), max(fe$redundancy_mean)))
cat(sprintf("Fig5: %d link sets, %d SEA-AD libraries\n", nrow(or_sets), nrow(lib)))
cat("-> ", OUT, "\n")
