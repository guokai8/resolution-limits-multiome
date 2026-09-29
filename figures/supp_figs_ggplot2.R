# Supplementary Figures 1-5
# 用法: Rscript supp_figs_ggplot2.R <derived_results_dir> <out_dir>
args <- commandArgs(trailingOnly = TRUE)
DR  <- if (length(args) >= 1) args[1] else "data/derived_results"
OUT <- if (length(args) >= 2) args[2] else "figures"
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
      value = TRUE)[1])), "fig_common.R"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

mh_pool <- function(df, by) {
  df |>
    mutate(a = prox_links, b = n_links - prox_links,
           cc = prox_tested_all, d = n_tested_all - prox_tested_all,
           tot = a + b + cc + d, R = a * d / tot, S = b * cc / tot,
           P = (a + d) / tot, Q = (b + cc) / tot) |>
    group_by(across(all_of(by))) |>
    summarise(Rs = sum(R), Ss = sum(S), PR = sum(P * R), PS = sum(P * S),
              QR = sum(Q * R), QS = sum(Q * S), .groups = "drop") |>
    mutate(OR = Rs / Ss,
           se = sqrt(PR / (2 * Rs^2) + (PS + QR) / (2 * Rs * Ss) + QS / (2 * Ss^2)),
           lo = OR * exp(-1.96 * se), hi = OR * exp(1.96 * se))
}
save2 <- function(pl, name, w, h) {
  agg_png(file.path(OUT, paste0(name, ".png")), width = w, height = h,
          units = "mm", res = 450); print(pl); invisible(dev.off())
  pdf(file.path(OUT, paste0(name, ".pdf")), width = w/25.4, height = h/25.4)
  print(pl); invisible(dev.off())
}

## ---- SF1 · 核数梯队按细胞类型拆开 -------------------------------------------
pri <- read_csv(file.path(DR, "primary_L3_nucleus_ladder.csv"), show_col_types = FALSE) |>
  mh_pool(c("celltype", "n")) |> mutate(cohort = "Primary (motor cortex)")
ext <- read_csv(file.path(DR, "nabec_L3_nucleus_ladder.csv"), show_col_types = FALSE) |>
  mh_pool(c("celltype", "n")) |> mutate(cohort = "External (prefrontal cortex)")
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
  labs(x = "nuclei per cell type", y = "enrichment odds ratio") +
  theme_pub() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 6.5),
        strip.background = element_blank(), strip.text = element_text(size = 8.5))
save2(sf1, "SuppFig1_ladder_by_celltype", 180, 78)

## ---- SF2 · ATAC 深度序列 ----------------------------------------------------
rep <- read_csv(file.path(DR, "nabec_L3_external_replication.csv"), show_col_types = FALSE)
dep <- rep |> filter(grepl("da[0-9]+", file), n == 150) |>
  mutate(depth = as.numeric(sub(".*da([0-9]+).*", "\\1", file)))
base <- rep |> filter(!grepl("da[0-9]+|fine", file), n == 150,
                      celltype %in% c("Oligo", "ExN")) |> mutate(depth = 5564)
dd <- bind_rows(dep, base) |> mh_pool("depth") |> mutate(fold = depth / 5564)
sf2 <- ggplot(dd, aes(fold, OR)) +
  geom_hline(yintercept = 1, linetype = "22", linewidth = .35, colour = "grey45") +
  geom_line(colour = C_ULM, linewidth = .45) +
  geom_linerange(aes(ymin = lo, ymax = hi), colour = C_ULM, linewidth = .4) +
  geom_point(size = 1.6, colour = C_ULM) +
  scale_x_log10(breaks = c(1, 2, 4, 8), labels = c("1x", "2x", "4x", "8x")) +
  scale_y_log10(limits = c(.3, 2), breaks = c(.4, .6, 1, 1.5)) +
  labs(x = "ATAC depth, relative to the equalised design (5,564 fragments)",
       y = "enrichment odds ratio") +
  theme_pub()
save2(sf2, "SuppFig2_depth_series", 90, 70)

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
  labs(x = NULL, y = "enrichment odds ratio, same nuclei, same depth (5,564 ATAC fragments)") +
  theme_pub() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 6.2), legend.key.height = unit(13, "pt"))
save2(sf3, "SuppFig3_convention_sensitivity", 120, 82)

## ---- SF4 · 外部复现 11 个分层森林图 -----------------------------------------
fo <- rep |> filter(n == 150, !grepl("da[0-9]+", file)) |>
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
  labs(x = NULL, y = "enrichment odds ratio at 150 nuclei",
       subtitle = "shading: primary cohort, 0.663 [0.568, 0.768]") +
  theme_pub() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 6.5),
        plot.subtitle = element_text(size = 6.5, colour = "grey35", margin = margin(b = 4)))
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
  labs(x = "donor-level difference / multinomial floor",
       y = "contrasts", fill = NULL) +
  theme_pub() +
  theme(legend.position = c(.98, .98), legend.justification = c(1, 1),
        legend.text = element_text(size = 6.5))
save2(sf5, "SuppFig5_floor_reassessment", 110, 74)

cat("Supplementary Figures 1-5 written\n")
