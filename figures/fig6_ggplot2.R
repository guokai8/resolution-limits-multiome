#!/usr/bin/env Rscript
# Fig 6 · 框架图：分辨率由断言的统计单元决定。
# 用法: Rscript fig6_ggplot2.R <derived_results_dir> <out_dir>
#
# 按 Genome Biology 的要求重做：这张图是概念核心，不再是三组互不相干的曲线。
#   a  生物学问题 → 统计单元 → 分辨率估计量 → 实验需求
#   b  三类断言各自的统计单元与分辨率量
#   c  三者在同一条「每细胞类型核数」轴上的紧凑对照（不画全部曲线）
args <- commandArgs(trailingOnly = TRUE)
OUT  <- ifelse(length(args) >= 2, args[2], "figures")
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])),
                 "fig_common.R"))

## ---- a  流程 ---------------------------------------------------------------
step <- data.frame(
  y   = c(4, 3, 2, 1),
  lab = c("Biological question", "Statistical unit",
          "Resolution estimator", "Experimental requirement"),
  eg  = c("does this cell type change?", "donor proportion",
          "overdispersion (\u03ba)", "nuclei per cell type"))

pa <- ggplot(step) +
  geom_label(aes(1, y, label = lab), fill = "grey95", colour = "grey10",
             size = 2.4, fontface = "bold", family = BASE_FAMILY,
             label.size = 0, label.padding = unit(3, "pt")) +
  geom_text(aes(1, y - 0.33, label = eg), colour = "grey40",
            size = 2.0, family = BASE_FAMILY, fontface = "italic") +
  annotate("segment", x = 1, xend = 1, y = step$y[-4] - 0.46, yend = step$y[-1] + 0.2,
           colour = "grey55", linewidth = .45,
           arrow = arrow(length = unit(3.2, "pt"), type = "closed")) +
  scale_x_continuous(limits = c(0.3, 1.7)) +
  scale_y_continuous(limits = c(0.55, 4.45)) +
  labs(tag = "a") +
  theme_blank()

## ---- b  三类断言的对照表 ----------------------------------------------------
tab <- data.frame(
  y    = c(3, 2, 1),
  claim = c("Composition", "Expression", "Regulation"),
  unit  = c("donor", "donor pseudobulk", "nucleus + method"),
  meas  = c("\u03ba", "sampling floor", "promoter enrichment"),
  col   = c(L_COMP, L_EXPR, L_REG))

pb <- ggplot(tab) +
  annotate("segment", x = 0.04, xend = 4.85, y = 3.52, yend = 3.52,
           colour = "grey35", linewidth = .4) +
  annotate("text", x = c(0.04, 1.58, 3.22), y = 3.74, hjust = 0, size = 1.95,
           colour = "grey25", fontface = "bold", family = BASE_FAMILY,
           label = c("Biological claim", "Statistical unit", "Resolution measure")) +
  geom_text(aes(0.04, y, label = claim, colour = I(col)), hjust = 0,
            size = 2.05, fontface = "bold", family = BASE_FAMILY) +
  geom_text(aes(1.58, y, label = unit), hjust = 0, colour = "grey25",
            size = 2.05, family = BASE_FAMILY) +
  geom_text(aes(3.22, y, label = meas), hjust = 0, colour = "grey25",
            size = 2.05, family = BASE_FAMILY) +
  annotate("segment", x = 0.04, xend = 4.85, y = 0.6, yend = 0.6,
           colour = "grey75", linewidth = .3) +
  scale_x_continuous(limits = c(0, 4.90)) +
  scale_y_continuous(limits = c(0.5, 3.95)) +
  labs(tag = "b") +
  theme_blank()

## ---- c  同一轴上的紧凑对照 --------------------------------------------------
## 只标各层的操作区间，不画任何曲线。数值全部取自正文。
band <- data.frame(
  layer = factor(c("Expression", "Composition", "Regulation"),
                 levels = c("Regulation", "Composition", "Expression")),
  lo    = c(17,     44435,  400),
  hi    = c(91,    499575, 5000),
  mid   = c(39,    125949,  900),
  note  = c("measured support, median effective n = 39",
            "1 percentage point at measured \u03ba",
            "promoter-enrichment transition"),
  col   = c(L_EXPR, L_COMP, L_REG))

pc <- ggplot(band, aes(y = layer, colour = I(col))) +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 2.6, alpha = .32) +
  geom_point(aes(x = mid), size = 2.2) +
  geom_text(aes(x = sqrt(lo * hi), label = note, hjust = c(0, 1, 0.5)),
            vjust = -1.6, size = 1.95, family = BASE_FAMILY, colour = "grey30") +
  scale_x_log10(limits = c(10, 1.1e6),
                breaks = c(10, 100, 1000, 10000, 1e5, 1e6),
                labels = c("10", "100", "1,000", "10,000", expression(10^5), expression(10^6))) +
  annotation_logticks(sides = "b", size = .25,
                      short = unit(1, "pt"), mid = unit(1.6, "pt"), long = unit(2.4, "pt")) +
  labs(tag = "c", x = "Nuclei per cell type", y = NULL) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 7.2, colour = "grey15"))

fig <- (pa | pb | pc) + plot_layout(widths = c(0.72, 1.18, 1.18))
save_fig(fig, OUT, "Fig6_resolution_limits", width = 180, height = 62)
cat("Fig6: framework (flow / claim table / common axis)\n")
