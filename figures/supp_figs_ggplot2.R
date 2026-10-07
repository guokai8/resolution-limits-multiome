# Supplementary Figures 1-5
# 用法: Rscript supp_figs_ggplot2.R <derived_results_dir> <out_dir>
args <- commandArgs(trailingOnly = TRUE)
DR  <- if (length(args) >= 1) args[1] else "data/derived_results"
OUT <- if (length(args) >= 2) args[2] else "figures"
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
      value = TRUE)[1])), "fig_common.R"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## 基础 pdf() 设备不认 Arial，矢量输出必须走 cairo_pdf（与主图 save_fig 一致）
save2 <- function(pl, name, w, h) {
  agg_png(file.path(OUT, paste0(name, ".png")), width = w, height = h,
          units = "mm", res = 450); print(pl); invisible(dev.off())
  cairo_pdf(file.path(OUT, paste0(name, ".pdf")), width = w/25.4, height = h/25.4)
  print(pl); invisible(dev.off())
}

## ---- SF1 · 核数梯队按细胞类型拆开 -------------------------------------------
pri <- read_csv(file.path(DR, "primary_L3_nucleus_ladder.csv"), show_col_types = FALSE) |>
  mh_repeats(c("celltype", "n")) |> mutate(cohort = "Primary (motor cortex)")
ext <- read_csv(file.path(DR, "nabec_L3_nucleus_ladder.csv"), show_col_types = FALSE) |>
  mh_repeats(c("celltype", "n")) |> mutate(cohort = "External (prefrontal cortex)")
lad <- bind_rows(pri, ext) |>
  mutate(cohort = factor(cohort, levels = c("Primary (motor cortex)",
                                            "External (prefrontal cortex)")))
sf1 <- ggplot(lad, aes(n, OR, colour = celltype)) +
  geom_hline(yintercept = 1, linetype = "22", linewidth = .35, colour = "grey45") +
  geom_line(linewidth = .45) +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = .35) +
  geom_point(size = 1.3) +
  facet_wrap(~cohort) +
  scale_colour_manual(values = c(C_ULM, C_SEA, C_MID, C_NULL)) +
  scale_x_log10(breaks = c(150, 400, 900, 2400, 7500), labels = label_comma(accuracy = 1)) +
  scale_y_log10(breaks = c(.25, .5, 1, 2, 4, 8)) +
  labs(x = "Nuclei per cell type", y = "Promoter-enrichment odds ratio") +
  theme_pub() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 6.5),
        strip.background = element_blank(), strip.text = element_text(size = 8.5))
save2(sf1, "SuppFig1_ladder_by_celltype", 180, 78)

## ---- SF2 · ATAC 深度序列 ----------------------------------------------------
## 直接取 p5_17 的 B8 沉淀值：那里的种子已按重复计算抽样合并，与正文一致。
## 外部复现表把种子编码在文件名里而没有 seed 列，就地重算会把种子当独立层。
dd <- read_csv(file.path(DR, "p5_17_B8_depth_series.csv"), show_col_types = FALSE) |>
  transmute(fold, OR = OR_seeds_as_repeats,
            lo = lo_seeds_as_repeats, hi = hi_seeds_as_repeats)
sf2 <- ggplot(dd, aes(fold, OR)) +
  geom_hline(yintercept = 1, linetype = "22", linewidth = .35, colour = "grey45") +
  geom_line(colour = C_ULM, linewidth = .45) +
  geom_linerange(aes(ymin = lo, ymax = hi), colour = C_ULM, linewidth = .4) +
  geom_point(size = 1.6, colour = C_ULM) +
  scale_x_log10(breaks = c(1, 2, 4, 8), labels = c("1x", "2x", "4x", "8x")) +
  ## 区间下界到 0.265、上界到 1.405，范围必须容纳按重复计算抽样合并后的宽度
  scale_y_log10(limits = c(.24, 1.6), breaks = c(.25, .5, 1, 1.5)) +
  ann(x = 1.05, y = 1.5, size = 2.4, colour = "grey15", hjust = 0,
      lab = "ATAC depth alone does not\nrecover regulatory structure") +
  labs(x = "ATAC depth, relative to the equalised design\n(5,564 fragments per nucleus)",
       y = "Promoter-enrichment odds ratio") +
  theme_pub()
save2(sf2, "SuppFig2_depth_series", 110, 72)

## ---- SF3 · 口径敏感性 -------------------------------------------------------
# ⚠️ 旧版把 `convention` 的两个水平当作口径对比，但它们的 depth_atac 相差 5–22 倍
#    （unmatched 27,973–120,798 对 p5_08 matched 5,564），那是**深度**对比，不是口径对比。
#    真正的口径对比在同一行内、同一深度上：OR（把不可检验的 peak 计入背景与 BH 分母）
#    对 OR_varpeaks（排除它们）。只取深度固定的那一组。
cv <- read_csv(file.path(DR, "nabec_L3_convention_sensitivity.csv"), show_col_types = FALSE) |>
  filter(convention == "p5_08 matched") |>
  transmute(celltype,
            incl_OR = OR,  incl_lo = lo,  incl_hi = hi,
            excl_OR = OR_varpeaks, excl_lo = lo_varpeaks, excl_hi = hi_varpeaks) |>
  pivot_longer(-celltype,
               names_to = c("conv", ".value"),
               names_pattern = "(incl|excl)_(.*)") |>
  mutate(convention = factor(conv, levels = c("excl", "incl"),
                             labels = c("background and FDR denominator\nexclude non-callable peaks",
                                        "both include them (as in this study)")))
sf3 <- ggplot(cv, aes(reorder(celltype, OR), OR, colour = convention)) +
  geom_hline(yintercept = 1, linetype = "22", linewidth = .35, colour = "grey45") +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = .45,
                 position = position_dodge(width = .55)) +
  geom_point(size = 1.8, position = position_dodge(width = .55)) +
  coord_flip() +
  scale_colour_manual(values = c(C_SEA, C_ULM)) +
  scale_y_log10(breaks = c(.5, 1, 2)) +
  labs(x = NULL, y = "Promoter-enrichment odds ratio\n(same nuclei, same depth)") +
  theme_pub() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 6.2), legend.key.height = unit(13, "pt"))
save2(sf3, "SuppFig3_convention_sensitivity", 120, 82)

## ---- SF4 · 外部复现 11 个分层森林图 -----------------------------------------
## 逐分层的 OR 来自沉淀表自身，不跨种子合并，所以这里直接读取。
ext_rep <- read_csv(file.path(DR, "nabec_L3_external_replication.csv"),
                    show_col_types = FALSE)
fo <- ext_rep |> filter(n == 150, !grepl("da[0-9]+", file)) |>
  mutate(gran = ifelse(grepl("fine", file), "fine subtype", "cell class"),
         lab  = paste0(celltype, ifelse(gran == "fine subtype", " (cluster)", "")))
sf4 <- ggplot(fo, aes(reorder(lab, OR), OR, colour = gran)) +
  geom_hline(yintercept = 1, linetype = "22", linewidth = .35, colour = "grey45") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.568, ymax = 0.768,
           fill = C_LITE, alpha = .45) +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = .45) +
  geom_point(size = 1.7) + coord_flip() +
  scale_colour_manual(values = c("cell class" = C_NULL, "fine subtype" = C_ULM)) +
  scale_y_log10(breaks = c(.3, .5, 1, 2)) +
  scale_x_discrete(expand = expansion(add = c(.6, 1.75))) +
  # 阴影带就是主分析 0.663 [0.568, 0.768] 的区间，必须在图内直接说明
  ann(x = nrow(fo) + .85, y = .663, lab = "primary cohort,\n150 nuclei (95% CI)",
      size = 2.0, colour = "grey30", hjust = .5, vjust = .5) +
  labs(x = NULL, y = "Promoter-enrichment odds ratio at 150 nuclei") +
  theme_pub() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 6.5))
save2(sf4, "SuppFig4_external_replication_forest", 110, 90)

## ---- SF5 · 190 个组成对比：合并检验 vs 供体层 vs 下限 ------------------------
fc <- read_csv(file.path(DR, "floor_reassessment_contrasts.csv"), show_col_types = FALSE)
fc2 <- fc |>
  filter(!is.na(delta_over_floor_multinomial), delta_over_floor_multinomial > 0) |>
  mutate(ratio = delta_over_floor_multinomial,
         cls = case_when(p_naive_pooled < .05 & p_value < .05 ~ "significant by both",
                         p_naive_pooled < .05 ~ "cell-pooled only",
                         TRUE ~ "neither"))
sf5 <- ggplot(fc2, aes(ratio, fill = cls)) +
  geom_vline(xintercept = 1, linetype = "22", linewidth = .35, colour = "grey30") +
  geom_histogram(bins = 42, colour = "white", linewidth = .15) +
  scale_fill_manual(values = c("cell-pooled only" = C_SEA, "significant by both" = C_ULM,
                               "neither" = C_GREY)) +
  scale_x_log10(breaks = c(.2, 1, 5, 20, 100)) +
  scale_y_continuous(expand = expansion(mult = c(0, .14))) +
  labs(x = "Donor-level difference / multinomial expectation",
       y = "Contrasts", fill = NULL) +
  theme_pub() +
  theme(legend.position = c(.02, .98), legend.justification = c(0, 1),
        legend.text = element_text(size = 6.5))
save2(sf5, "SuppFig5_floor_reassessment", 110, 74)

cat("Supplementary Figures 1-6 written\n")

## ---- SF6 · 检出 z 分析（按 GB 要求从主图 Fig 4 移来）-------------------------
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

sf6 <- ggplot(pw, aes(nl, z, colour = or)) +
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
  labs(x = "Detected links", y = "Detection z-score") +
  theme_pub()

save2(sf6, "SuppFig6_detection_power", 95, 72)
