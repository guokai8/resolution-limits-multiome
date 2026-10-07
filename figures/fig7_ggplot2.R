#!/usr/bin/env Rscript
# Fig 7 · 核数依赖是否随 linker 改变。
# 用法: Rscript fig7_ggplot2.R <supplementary_tables_dir> <out_dir>
#
# a 启动子富集 OR 对核数，三个 linker 臂 × 两个细胞类型
# b 位置 OR（相对 trans 配对零假设），标出 >=200 链接的可读下界
# c 聚合体大小 k 的敏感性：同一批数据、同一核数，链接数与 OR 随 k 变化
# d trans 零假设的 9 次独立重复与留一染色体 jackknife
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(readr)
                                library(patchwork); library(scales); library(ragg)})
args <- commandArgs(trailingOnly = TRUE)
ST  <- ifelse(length(args) >= 1, args[1], "supplementary_tables")
OUT <- ifelse(length(args) >= 2, args[2], "figures")
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])),
                 "fig_common.R"))

t20 <- read_csv(file.path(ST, "Supplementary_Table_20_linker_comparison.csv"),
                show_col_types = FALSE)
t21 <- read_csv(file.path(ST, "Supplementary_Table_21_trans_null_validation.csv"),
                show_col_types = FALSE)

## 三个臂的命名与配色：正文口径的单核相关、ArchR 自己的阈值、以及
## 只用 FDR 的松阈值聚合（后者用来说明阈值本身的作用）
arm_of <- function(linker, archr) {
  ifelse(linker == "single_cell", "Single-nucleus correlation",
         ifelse(archr, "Aggregation + ArchR-like retrieval", "Aggregation + FDR only"))
}
ARMC <- c("Single-nucleus correlation"              = C_ULM,
          "Aggregation + ArchR-like retrieval"      = C_SEA,
          "Aggregation + FDR only"                  = C_GREY)
CT <- c("Exc_LINC00507_FREM3" = 16, "Oligodendrocytes" = 17)

main <- t20 |>
  filter(window == 500000, k %in% c(0, 25)) |>
  mutate(arm = arm_of(linker, archr_defaults),
         arm = factor(arm, levels = names(ARMC)),
         celltype = factor(celltype, levels = names(CT)))

## ---- a · 启动子富集 OR 对核数 ----------------------------------------------
pa <- ggplot(main, aes(n, promoter_OR, colour = arm, shape = celltype)) +
  geom_hline(yintercept = 1, linetype = "22", linewidth = .35, colour = "grey45") +
  geom_vline(xintercept = 150, linewidth = .3, colour = "grey55") +
  geom_line(aes(group = interaction(arm, celltype)), linewidth = .45, alpha = .9) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = ARMC, name = NULL,
                      limits = names(ARMC), drop = FALSE) +
  scale_shape_manual(values = CT, name = NULL,
                     limits = names(CT), drop = FALSE) +
  scale_x_log10(breaks = c(150, 400, 900, 2400, 7500), labels = label_comma(accuracy = 1)) +
  scale_y_log10(breaks = c(.5, 1, 2, 5, 10, 30)) +
  ann(x = 165, y = 34, lab = "Same nuclei,\ndifferent procedures", size = 2.4,
      colour = "grey15") +
  labs(x = "Nuclei per cell type", y = "Promoter-enrichment odds ratio") +
  guides(colour = guide_legend(nrow = 1, order = 1),
         shape  = guide_legend(nrow = 1, order = 2)) +
  theme_pub()

## ---- b · 位置 OR（相对 trans 零假设）---------------------------------------
pb <- main |>
  filter(links_obs >= 200) |>
  ggplot(aes(n, positional_OR, colour = arm, shape = celltype)) +
  geom_hline(yintercept = 1, linetype = "22", linewidth = .35, colour = "grey45") +
  geom_vline(xintercept = 150, linewidth = .3, colour = "grey55") +
  geom_linerange(aes(ymin = positional_lo, ymax = positional_hi), linewidth = .4) +
  geom_line(aes(group = interaction(arm, celltype)), linewidth = .4, alpha = .8) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = ARMC, name = NULL,
                      limits = names(ARMC), drop = FALSE) +
  scale_shape_manual(values = CT, name = NULL,
                     limits = names(CT), drop = FALSE) +
  scale_x_log10(breaks = c(150, 400, 900, 2400, 7500), labels = label_comma(accuracy = 1)) +
  scale_y_log10(breaks = c(.5, 1, 2, 5, 10, 20)) +
  ann(x = 160, y = 17, size = 2.25, colour = "grey30",
      lab = "Aggregation reduces sparsity but\ndoes not guarantee positional structure") +
  labs(x = "Nuclei per cell type", y = "Positional odds ratio") +
  guides(colour = "none", shape = "none") +
  theme_pub()

## ---- c · 聚合体大小 k 的敏感性 ---------------------------------------------
ks <- t20 |>
  filter(window == 500000, archr_defaults, n == 2400,
         celltype == "Exc_LINC00507_FREM3") |>
  select(k, links_obs, promoter_OR) |> arrange(k)
pc <- ggplot(ks, aes(links_obs, promoter_OR)) +
  geom_line(colour = C_SEA, linewidth = .45) +
  geom_point(size = 1.6, colour = C_SEA) +
  geom_text(aes(label = paste0("k=", k)), hjust = -0.18, vjust = -0.5,
            size = 2.2, colour = "grey30") +
  scale_x_log10(breaks = c(2, 10, 70, 500, 3000), labels = label_comma(accuracy = 1)) +
  scale_y_log10(breaks = c(3, 6, 10, 20)) +
  expand_limits(x = 6000, y = 26) +
  labs(x = "Links recovered", y = "Promoter-enrichment odds ratio") +
  theme_pub()

## ---- d · trans 零假设的稳定性 ----------------------------------------------
## 表 21 的 arm 列用中文记了 ArchR 臂，这里映射成图上用的英文标签
rep9 <- t21 |> filter(!is.na(positional_OR_median)) |>
  mutate(armlab = ifelse(grepl("ArchR", arm, fixed = TRUE),
                         "Aggregation + ArchR-like retrieval", "Single-nucleus correlation"),
         lab = paste0(armlab, ", n=", comma(n))) |>
  arrange(n, armlab) |>
  mutate(lab = factor(lab, levels = rev(unique(lab))))
pd <- ggplot(rep9, aes(positional_OR_median, lab)) +
  geom_vline(xintercept = 1, linetype = "22", linewidth = .35, colour = "grey45") +
  geom_linerange(aes(xmin = positional_OR_min, xmax = positional_OR_max),
                 linewidth = .45, colour = C_MID) +
  geom_point(size = 1.5, colour = C_ULM) +
  scale_x_log10(breaks = c(.5, 1, 2, 5, 20)) +
  ann(x = 1.15, y = 6.6, lab = "Procedure changes\npositional structure",
      size = 2.3, colour = "grey25") +
  labs(x = "Positional odds ratio", y = NULL) +
  theme_pub() + theme(axis.text.y = element_text(size = 6))

fig <- (pa | pb) / (pc | pd) +
  plot_layout(heights = c(1, .9), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.margin = margin(2, 0, 0, 0),
        legend.text = element_text(size = 6.2),
        legend.key.width = unit(10, "pt"),
        legend.key.height = unit(8, "pt"))
save_fig(fig, OUT, "Fig7_linker_comparison", width = 180, height = 160)
cat(sprintf("Fig7: %d arm-by-celltype series, k sweep %d points, %d validation rows\n",
            nrow(distinct(main, arm, celltype)), nrow(ks), nrow(rep9)))
