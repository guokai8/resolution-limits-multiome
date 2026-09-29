# Fig 6 · 三层推断的分辨率极限，放在同一根「每细胞类型核数」轴上
# 用法: Rscript fig6_ggplot2.R <derived_results_dir> <out_dir>
args <- commandArgs(trailingOnly = TRUE)
DR  <- if (length(args) >= 1) args[1] else "data/derived_results"
OUT <- if (length(args) >= 2) args[2] else "figures"
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
      value = TRUE)[1])), "fig_common.R"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- a. 调控层：两个队列的核数梯队 ------------------------------------------
mh_pool <- function(df, by) {
  df |>
    mutate(a = prox_links, b = n_links - prox_links,
           cc = prox_tested_all, d = n_tested_all - prox_tested_all,
           tot = a + b + cc + d,
           R = a * d / tot, S = b * cc / tot,
           P = (a + d) / tot, Q = (b + cc) / tot) |>
    group_by(across(all_of(by))) |>
    summarise(Rs = sum(R), Ss = sum(S), PR = sum(P * R), PS = sum(P * S),
              QR = sum(Q * R), QS = sum(Q * S), .groups = "drop") |>
    mutate(OR = Rs / Ss,
           se = sqrt(PR / (2 * Rs^2) + (PS + QR) / (2 * Rs * Ss) + QS / (2 * Ss^2)),
           lo = OR * exp(-1.96 * se), hi = OR * exp(1.96 * se))
}

pri <- read_csv(file.path(DR, "primary_L3_nucleus_ladder.csv"), show_col_types = FALSE) |>
  mh_pool("n") |> mutate(cohort = "Primary (motor cortex)")
ext <- read_csv(file.path(DR, "nabec_L3_nucleus_ladder.csv"), show_col_types = FALSE) |>
  mh_pool("n") |> mutate(cohort = "External (prefrontal cortex)")
lad <- bind_rows(pri, ext) |>
  mutate(cohort = factor(cohort, levels = c("Primary (motor cortex)",
                                            "External (prefrontal cortex)")))
COHC <- c("Primary (motor cortex)" = C_ULM, "External (prefrontal cortex)" = C_SEA)

pa <- ggplot(lad, aes(n, OR, colour = cohort)) +
  annotate("rect", xmin = 130, xmax = 9000, ymin = 2.26, ymax = 3.56,
           fill = C_LITE, alpha = .3) +
  annotate("text", x = 158, y = 3.15, label = "full-depth link sets",
           hjust = 0, size = 2.3, colour = C_ULM) +
  geom_hline(yintercept = 1, linetype = "22", linewidth = .35, colour = "grey45") +
  geom_vline(xintercept = 150, linewidth = .3, colour = "grey55") +
  geom_line(linewidth = .5) +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = .4) +
  geom_point(size = 1.4) +
  annotate("text", x = 150, y = 11.5, label = "equalised\ndesign", hjust = -0.1,
           vjust = 1, size = 2.3, colour = "grey35", lineheight = .95) +
  scale_colour_manual(values = COHC) +
  scale_x_log10(limits = c(130, 9000), breaks = c(150, 400, 900, 2400, 7500),
                labels = label_comma(accuracy = 1)) +
  scale_y_log10(breaks = c(.5, 1, 2, 4, 8), labels = c("0.5", "1", "2", "4", "8")) +
  labs(tag = "a", title = "Regulatory: promoter enrichment",
       x = "nuclei per cell type", y = "enrichment odds ratio") +
  theme_pub() +
  theme(legend.position = c(.97, .03), legend.justification = c(1, 0),
        legend.title = element_blank(), legend.text = element_text(size = 6.2),
        legend.key.width = unit(11, "pt"))

## ---- b. 表达层：下限对核数 --------------------------------------------------
# 50 核配对准则同样管表达层：ALS26 第二块芯片只有 37 个核，不合格
ex <- read_csv(file.path(DR, "p5_floor_scaling_pairs.csv"), show_col_types = FALSE) |>
  filter(donor != "ALS26") |>
  filter(median_abs_log2FC > 0, n_eff > 0)
cfe  <- coef(lm(log(median_abs_log2FC) ~ log(n_eff), data = ex))
linee <- tibble(n = exp(seq(log(8), log(2000), length.out = 100))) |>
  mutate(f = exp(cfe[1] + cfe[2] * log(n)))

pb <- ggplot(ex, aes(n_eff, median_abs_log2FC)) +
  geom_point(size = .5, alpha = .28, colour = C_GREY) +
  geom_line(data = linee, aes(n, f), colour = C_ULM, linewidth = .5, inherit.aes = FALSE) +
  scale_x_log10(breaks = c(10, 50, 200, 1000), labels = label_comma(accuracy = 1)) +
  scale_y_log10(breaks = c(.1, .25, .5, 1, 2)) +
  labs(tag = "b", title = "Expression: technical-noise threshold",
       x = "effective nuclei per cell type",
       y = "threshold, absolute log2 fold change") +
  annotate("text", x = 11, y = .13, hjust = 0, size = 2.5, colour = C_ULM,
           label = sprintf("floor = %.2f n^%.3f", exp(cfe[1]), cfe[2])) +
  theme_pub()

## ---- c. 组成层：最小可检测差异 ----------------------------------------------
kap <- read_csv(file.path(DR, "p5_17_B7_cohort_overdispersion.csv"),
                show_col_types = FALSE)
k_sea <- round(kap$kappa[kap$cohort == "Seattle atlas"], 2)
k_mc  <- round(kap$kappa[kap$cohort == "Motor cortex"], 2)
gridc <- expand.grid(N = 10^seq(log10(2e3), log10(2e6), length.out = 80),
                     kappa = c(1, k_sea, k_mc)) |>
  mutate(F = 100 * 1.96 * kappa * sqrt(.1 * .9) * sqrt(2 / N),
         lab = factor(kappa, levels = c(1, k_sea, k_mc),
                      labels = c("kappa = 1 (multinomial)",
                                 sprintf("kappa = %.2f (Seattle)", k_sea),
                                 sprintf("kappa = %.2f (motor cortex)", k_mc))))
pc <- ggplot(gridc, aes(N, F, colour = lab, linetype = lab)) +
  geom_line(linewidth = .5) +
  geom_hline(yintercept = 1, linetype = "22", linewidth = .35, colour = "grey45") +
  scale_colour_manual(values = c(C_NULL, C_SEA, C_ULM)) +
  scale_linetype_manual(values = c("42", "solid", "solid")) +
  scale_x_log10(breaks = c(1e4, 1e5, 1e6), labels = c("10k", "100k", "1M")) +
  scale_y_log10(breaks = c(.1, .3, 1, 3, 10)) +
  labs(tag = "c", title = "Composition: technical-noise threshold",
       x = "nuclei per group",
       y = "threshold, percentage points") +
  theme_pub() +
  theme(legend.position = c(.03, .06), legend.justification = c(0, 0),
        legend.title = element_blank(), legend.text = element_text(size = 6.4),
        legend.key.width = unit(13, "pt"))

fig <- (pa | pb | pc) + plot_layout(widths = c(1.12, 1, 1.02))
agg_png(file.path(OUT, "Fig6_resolution_limits.png"), width = 180, height = 64,
        units = "mm", res = 450); print(fig); invisible(dev.off())
pdf(file.path(OUT, "Fig6_resolution_limits.pdf"), width = 180/25.4, height = 64/25.4)
print(fig); invisible(dev.off())

cross <- function(d, t) {
  d <- d[order(d$n), ]
  out <- NA_real_
  for (i in seq_len(nrow(d) - 1)) {          # 取最后一次上穿（曲线可非单调）
    if (d$OR[i] < t && d$OR[i + 1] >= t) {
      b <- log(d$OR[i + 1] / d$OR[i]) / log(d$n[i + 1] / d$n[i])
      out <- exp(log(d$n[i]) + log(t / d$OR[i]) / b)
    }
  }
  out
}
for (co in levels(lad$cohort)) {
  d <- lad[lad$cohort == co, ]
  cat(sprintf("%-30s OR=1 at n=%.0f;  OR=2.26 at n=%.0f\n", co, cross(d, 1), cross(d, 2.26)))
}
cat(sprintf("expression exponent %.3f\n", cfe[2]))
