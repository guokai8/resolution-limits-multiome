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
           label = paste0(n_rep, " donors on two chips\n(26 with >=50 nuclei on both)")) +
  scale_fill_gradient(low = "#F2F7FC", high = C_ULM) +
  scale_x_discrete(labels = function(x) sub("Chip", "", x)) +
  labs(title = "Donor x chip layout", x = "Chip", y = "Donor (ordered)") +
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
dims <- c(n_nuclei = "Nuclei", median_nCount_RNA = "RNA",
          median_nCount_ATAC = "ATAC", pmito = "Mito %")
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
  labs(title = "Quality is unbiased", x = NULL, y = "Value / median") +
  theme_pub()

lay <- data.frame(
  lab  = c("L1  Composition", "L2  Expression", "L3  Regulation"),
  need = c("cell assignment", "aggregated profile", "per-element, per-cell"),
  fill = c(C_LITE, C_MID, C_ULM),
  txt  = c(C_ULM, "white", "white"),
  ymin = c(2.1, 1.1, 0.1)) %>% mutate(ymax = ymin + .8)

p1d <- ggplot(lay) +
  geom_rect(aes(xmin = 0, xmax = 1.12, ymin = ymin, ymax = ymax, fill = I(fill))) +
  geom_text(aes(.06, ymin + .52, label = lab, colour = I(txt)),
            hjust = 0, size = 2.15, fontface = "bold") +
  geom_text(aes(.06, ymin + .22, label = need, colour = I(txt)),
            hjust = 0, size = 1.9) +
  annotate("segment", x = 1.26, xend = 1.26, y = 2.9, yend = 0.15,
           arrow = arrow(length = unit(4, "pt"), ends = "last", type = "closed"),
           colour = "grey45", linewidth = .45) +
  annotate("text", x = 1.40, y = 1.5, label = "increasing\ndemand", angle = 90,
           size = 2.1, colour = "grey35", lineheight = 1) +
  scale_x_continuous(limits = c(0, 1.52)) +
  scale_y_continuous(limits = c(0, 3.2)) +
  labs(title = "Three layers") +
  theme_blank()

fig1 <- p1a + p1b + p1d +
  plot_layout(widths = c(1.15, 1, 0.95)) +
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
  annotate("text", x = nrow(nsub) / 2, y = 162, size = 2.2, colour = "grey20",
           label = "150 nuclei per subtype") +
  annotate("text", x = .5, y = 118, hjust = 0, size = 2.0, colour = "white",
           lineheight = 1, label = "depth ratio 1.0000\nno nucleus below target") +
  scale_y_continuous(limits = c(0, 190), expand = expansion(mult = c(0, .02))) +
  labs(title = "Equalised design", y = "Nuclei",
       x = paste0(nrow(nsub), " neuronal subtypes")) +
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
  labs(title = "Peaks are seen in ~2 nuclei",
       x = "Neuronal subtypes", y = "Nuclei per detected peak") +
  theme_pub() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# c: 实测 OR —— 显著耗竭（正文引 Fig. 4c）
orv <- data.frame(or = 0.663, lo = 0.568, hi = 0.768)
p4c <- ggplot(orv) +
  geom_col(aes(1, or), fill = C_SEA, width = .3) +
  geom_errorbar(aes(1, ymin = lo, ymax = hi), width = .09, linewidth = .5,
                colour = "grey20") +
  geom_hline(yintercept = 1, linetype = "22", linewidth = .45, colour = "grey20") +
  annotate("text", x = 1.26, y = 0.663, hjust = .5, size = 2.3, colour = C_SEA,
           label = "0.663 [0.568, 0.768]") +
  annotate("text", x = .72, y = 0.04, hjust = 0, size = 2.0, colour = "grey30",
           lineheight = 1, label = "significant depletion,\nnot absence of signal") +
  coord_flip(xlim = c(.6, 1.45), ylim = c(0, 1.3)) +
  labs(title = "Observed: depletion", x = NULL, y = "Odds ratio") +
  theme_pub() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# d: 功效曲线（正文引 Fig. 4d）
p0     <- 0.0095
n_test <- 3.03e6
n_obs  <- 33227
pw <- lapply(c(3.0, 2.0, 1.5), function(orr) {
  nl <- 10^seq(2, 5, length.out = 160)
  p1 <- orr * p0 / (1 - p0 + orr * p0)
  se <- sqrt(p1 * (1 - p1) / nl + p0 * (1 - p0) / n_test)
  data.frame(nl = nl, z = (p1 - p0) / se, or = factor(orr, levels = c(3, 2, 1.5)))
}) %>% bind_rows()
ztar <- c(`3` = 9.0, `2` = 6.2, `1.5` = 3.6)
tip <- pw %>% group_by(or) %>%
  slice_min(abs(z - ztar[as.character(or)]), n = 1) %>% ungroup()

p4d <- ggplot(pw, aes(nl, z, colour = or)) +
  geom_line(linewidth = .8) +
  geom_hline(yintercept = 1.96, linetype = "22", linewidth = .4, colour = "grey20") +
  geom_vline(xintercept = n_obs, colour = C_SEA, linewidth = .8) +
  geom_text(data = tip, aes(nl, z, label = paste("OR", or), colour = or),
            hjust = 1.15, vjust = -0.35, size = 2.0, show.legend = FALSE) +
  annotate("text", x = 120, y = 2.4, hjust = 0, size = 2.0, colour = "grey30",
           label = "P = 0.05") +
  annotate("text", x = n_obs * 1.25, y = 1.0, hjust = 0, size = 2.0, colour = C_SEA,
           lineheight = 1, label = "33,227\nlinks") +
  scale_colour_manual(values = c(`3` = C_LITE, `2` = C_MID, `1.5` = C_ULM)) +
  scale_x_log10(breaks = 10^(2:5), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(expand = expansion(mult = c(0, .02))) +
  coord_cartesian(ylim = c(0, 11)) +
  labs(title = "Power was ample", x = "Detected links", y = "Detection z-score") +
  theme_pub()

fig4 <- p4a + p4b + p4c + p4d +
  plot_layout(widths = c(.95, 1.15, .95, 1.1)) + plot_annotation(tag_levels = "a")
save_fig(fig4, OUT, "Fig4_regulatory_floor", height = 68)

## ======================================================================
## Fig 5
## ======================================================================
# 四个 cellranger-arc link 集直接读 p5_15 的产出；
# ⚠️ 另两行无对应 CSV：已发表 AD multiome 的补充表不在盘上，
#    本研究那一条由 L3 主分析给出，故按记录的数值补入。
ext <- grep("_w[0-9]+\\.csv$", list.files(RES, "^p5_promoterOR_.*\\.csv$", full.names = TRUE),
              value = TRUE, invert = TRUE) %>%
  lapply(read_csv, show_col_types = FALSE) %>% bind_rows() %>%
  transmute(lab = sub("_sorted", "", label), or = OR, lo, hi)
# 四个 cellranger-arc 参照集直接读 p5_15 的产出；本研究那一条由 L3 主分析给出。
or_sets <- bind_rows(
  ext %>% arrange(or),
  data.frame(lab = "This study (150 nuclei)", or = 0.663, lo = 0.568, hi = 0.768)) %>%
  mutate(lab = factor(lab, levels = rev(lab)), dep = or < 1)

p5a <- ggplot(or_sets, aes(or, lab, colour = dep)) +
  annotate("rect", xmin = 2.45, xmax = 4.01, ymin = -Inf, ymax = Inf,
           fill = C_ULM, alpha = .07) +   # 四个参照 link 集的实际区间
  geom_vline(xintercept = 1, linetype = "22", linewidth = .45, colour = "grey20") +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = .9) +
  geom_point(size = 1.5) +
  geom_text(aes(x = hi * 1.09, label = sprintf("%.2f", or)), hjust = 0, size = 2.2) +
  scale_colour_manual(values = c(`TRUE` = C_SEA, `FALSE` = C_ULM)) +
  scale_x_log10(limits = c(.45, 9), breaks = c(.5, 1, 2, 5),
                labels = c("0.5", "1", "2", "5")) +
  labs(title = "Reference range across link sets", y = NULL,
       x = "Promoter-enrichment odds ratio (log scale)") +
  theme_pub() +
  theme(axis.text.y = element_text(size = 6.5))

# 同一个 link 集 (human_brain_3k) 在两个窗口下的 OR，来自 p5_15 --window
win <- data.frame(w = factor(c("+/-1 Mb", "+/-1.7 Mb"),
                             levels = c("+/-1 Mb", "+/-1.7 Mb")),
                  or = c(4.012, 6.591))
p5b <- ggplot(win, aes(w, or, fill = w)) +
  geom_col(width = .55) +
  geom_text(aes(label = sprintf("%.2f", or)), vjust = -0.5, size = 2.6,
            colour = "grey15") +
  scale_fill_manual(values = c(C_ULM, C_LITE)) +
  scale_y_continuous(limits = c(0, 7.6), expand = expansion(mult = c(0, .02))) +
  labs(title = "Window matters", x = NULL, y = "Odds ratio, same link set") +
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
  labs(title = "No saturation",
       x = "Nuclei x fragments per nucleus", y = "Feature linkages detected") +
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
