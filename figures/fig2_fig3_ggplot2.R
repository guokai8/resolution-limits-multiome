#!/usr/bin/env Rscript
# Fig 2 与 Fig 3 重绘（ggplot2）
#
# Fig 2 全文核心视觉论证：同一版式，两个相反的结果
#   a  L2 两队列拟合线【平行】—— 斜率相同、截距不同（不是重合）
#   b  L1 两队列拟合线【分开】—— 斜率不同
#   c  过度离散幅度两队列一致（复现的那一半）
#   两个面板都叠一条 n^(-1/2) 参考斜率，使斜率对比不受版式与坐标跨度影响
#
# Fig 3
#   a  实测 vs 匹配的多项式零模型
#   b  逐对 实测/零模型 比值分布（中位 1.30x）
#   c  类型内指数，按 n_eff 跨度标注可determined 程度
#
# 数据来自 ~/Desktop/P5_VulnerableEpigenome/results/（原始分析产物）
# 用法： Rscript fig2_fig3_ggplot2.R [results_dir] [out_dir]

args    <- commandArgs(trailingOnly = TRUE)
RES     <- path.expand(if (length(args) >= 1) args[1] else "~/Desktop/P5_VulnerableEpigenome/results")
OUT     <- path.expand(if (length(args) >= 2) args[2] else "~/Desktop/ResearchD/figures")

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                  value = TRUE)[1])), "fig_common.R"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- 读数据 -----------------------------------------------------------
rd <- function(f) read_csv(file.path(RES, f), show_col_types = FALSE)

als_L2 <- rd("p5_floor_scaling_pairs.csv") %>%
  transmute(x = n_eff, y = median_abs_log2FC, cohort = "Motor cortex") %>%
  filter(x > 0, y > 0)
sea_L2 <- rd("seaad_L2_obs.csv") %>%
  transmute(x = n_eff, y = floor, cohort = "Seattle atlas") %>% filter(x > 0, y > 0)

als_L1 <- rd("p5_L1_scaling_pairs.csv") %>%
  transmute(x = n_eff, y = z, cohort = "Motor cortex") %>% filter(x > 0, y > 0)
sea_L1 <- rd("seaad_L1_pairs.csv") %>%
  transmute(x = n_eff, y = z, cohort = "Seattle atlas") %>% filter(x > 0, y > 0)

obs <- rd("seaad_L2_obs.csv")  %>% transmute(donor, ct, x = n_eff, y = floor)
nul <- rd("seaad_L2_null.csv") %>% transmute(donor, ct, x = n_eff, y = floor)
byct <- rd("seaad_L2_by_celltype.csv")

## ======================================================================
## Fig 2
## ======================================================================
L2 <- bind_rows(als_L2, sea_L2)
L1 <- bind_rows(als_L1, sea_L1)
xr <- range(c(L2$x, L1$x))                          # a/b 共用 x 轴范围

f2a_fit <- bind_rows(
  fit_line(als_L2, range(als_L2$x)) %>% mutate(cohort = "Motor cortex"),
  fit_line(sea_L2, range(sea_L2$x)) %>% mutate(cohort = "Seattle atlas"))
f2b_fit <- bind_rows(
  fit_line(als_L1, range(als_L1$x)) %>% mutate(cohort = "Motor cortex"),
  fit_line(sea_L1, range(sea_L1$x)) %>% mutate(cohort = "Seattle atlas"))

leg <- function(d, labs, xpos, ypos) {              # 面板内文字图例
  data.frame(cohort = names(labs), lab = unname(labs),
             x = xpos, y = ypos * (0.55^(seq_along(labs) - 1)))
}

f2a_leg <- leg(NULL, c("Motor cortex" = "Motor cortex  b = -0.507 [-0.569, -0.441]",
                       "Seattle atlas" = "Seattle atlas  b = -0.505 [-0.524, -0.483]"),
               xr[1] * 1.5, 0.075)
f2b_leg <- leg(NULL, c("Motor cortex" = "Motor cortex  b = -0.152 [-0.340, +0.111]",
                       "Seattle atlas" = "Seattle atlas  b = -0.584 [-0.838, -0.369]"),
               xr[1] * 1.5, 2.2e-4)

p2a <- ggplot(L2, aes(x, y, colour = cohort)) +
  geom_point(alpha = .13, size = .5, stroke = 0, show.legend = FALSE) +
  geom_line(data = pow_guide(25, 2.6), aes(x, y), inherit.aes = FALSE,
            colour = "grey55", linetype = "22", linewidth = .4) +
  annotate("text", x = 25 * 3.4, y = 2.6 * (3.4)^-0.5 * 1.5,
           label = "n^{-1/2}", parse = TRUE, size = 2.4, colour = "grey45") +
  geom_line(data = f2a_fit, aes(x, y, colour = cohort), linewidth = .85) +
  geom_text(data = f2a_leg, aes(x, y, label = lab, colour = cohort),
            hjust = 0, size = 2.35, show.legend = FALSE) +
  annotate("text", x = xr[2] * .92, y = 4.2, hjust = 1, vjust = 1, size = 2.5,
           colour = "grey25", lineheight = .95,
           label = "parallel:\nsame slope,\ndifferent intercept") +
  scale_colour_manual(values = COH) +
  scale_x_log10(limits = xr, breaks = 10^(0:4),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = 10^(-6:1),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "bl", size = .25,
                      short = unit(1, "pt"), mid = unit(1.6, "pt"), long = unit(2.4, "pt")) +
  labs(title = "Expression layer", x = expression("Effective number of nuclei ("*italic(n)[eff]*")"),
       y = "Expression floor") +
  theme_pub()

p2b <- ggplot(L1, aes(x, y, colour = cohort)) +
  geom_point(alpha = .13, size = .5, stroke = 0, show.legend = FALSE) +
  geom_line(data = pow_guide(60, 3e-3, span = 1.5), aes(x, y), inherit.aes = FALSE,
            colour = "grey55", linetype = "22", linewidth = .4) +
  annotate("text", x = 60 * 3.2, y = 3e-3 * (3.2)^-0.5 * 2.0,
           label = "n^{-1/2}", parse = TRUE, size = 2.4, colour = "grey45") +
  geom_line(data = f2b_fit, aes(x, y, colour = cohort), linewidth = .85) +
  geom_text(data = f2b_leg, aes(x, y, label = lab, colour = cohort),
            hjust = 0, size = 2.35, show.legend = FALSE) +
  annotate("text", x = xr[2] * .92, y = 2.2, hjust = 1, vjust = 1, size = 2.5,
           colour = "grey25", lineheight = .95,
           label = "divergent:\nslopes differ") +
  scale_colour_manual(values = COH) +
  scale_x_log10(limits = xr, breaks = 10^(0:4),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = 10^(-6:1),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "bl", size = .25,
                      short = unit(1, "pt"), mid = unit(1.6, "pt"), long = unit(2.4, "pt")) +
  labs(title = "Composition layer",
       x = expression("Effective number of nuclei ("*italic(n)[eff]*")"),
       y = "Composition discrepancy") +
  theme_pub()

# 四个队列，按重复所跨的步骤分组
lv <- c("Motor cortex", "Seattle atlas", "PsychAD MSSM", "PsychAD RADC")
od <- data.frame(cohort = factor(lv, levels = lv),
                 od = c(4.28, 3.90, 1.357, 1.220),
                 grp = c("incl. tissue + dissociation", "incl. separate sample",
                         "loading / library only", "loading / library only"))
p2c <- ggplot(od, aes(cohort, od, fill = grp)) +
  geom_col(width = .62) +
  geom_hline(yintercept = 1, linetype = "22", linewidth = .4, colour = "grey20") +
  geom_text(aes(label = sprintf("%.2f", od)), vjust = -0.6, size = 2.5,
            colour = "grey15") +
  annotate("segment", x = 2.62, xend = 4.38, y = 1.78, yend = 1.78,
           colour = "grey45", linewidth = .35) +
  annotate("text", x = 3.5, y = 1.83, vjust = 0, size = 2.0, colour = "grey35",
           label = "loading / library only") +
  scale_fill_manual(values = c(`incl. tissue + dissociation` = C_ULM,
                               `incl. separate sample` = C_MID,
                               `loading / library only` = C_LITE)) +
  scale_x_discrete(labels = c("Motor\ncortex", "Seattle\natlas", "PsychAD\nMSSM", "PsychAD\nRADC")) +
  scale_y_continuous(limits = c(0, 4.85), expand = expansion(mult = c(0, .02))) +
  labs(title = "Overdispersion", x = NULL,
       y = "Composition overdispersion") +
  theme_pub() +
  theme(axis.text.x = element_text(size = 6.4))

fig2 <- p2a + p2b + p2c + plot_layout(widths = c(1, 1, .95)) +
  plot_annotation(tag_levels = "a")
save_fig(fig2, OUT, "Fig2_transferability", height = 68)

## ======================================================================
## Fig 3
## ======================================================================
f3dat <- bind_rows(obs %>% mutate(set = "Observed"),
                   nul %>% mutate(set = "Sampling null")) %>% filter(x > 0, y > 0)
f3fit <- bind_rows(
  fit_line(obs %>% filter(x > 0, y > 0), range(obs$x)) %>% mutate(set = "Observed"),
  fit_line(nul %>% filter(x > 0, y > 0), range(nul$x)) %>% mutate(set = "Sampling null"))
SET <- c("Observed" = C_SEA, "Sampling null" = C_NULL)

f3a_leg <- data.frame(
  set = c("Observed", "Sampling null"),
  lab = c("Observed  b = -0.505", "Sampling null  b = -0.482"),
  x = min(f3dat$x) * 1.3, y = c(0.075, 0.056))

p3a <- ggplot(f3dat, aes(x, y, colour = set)) +
  geom_point(alpha = .11, size = .5, stroke = 0) +
  geom_line(data = f3fit, aes(x, y, colour = set), linewidth = .85) +
  geom_text(data = f3a_leg, aes(x, y, label = lab, colour = set),
            hjust = 0, size = 2.35) +
  scale_colour_manual(values = SET) +
  scale_x_log10(breaks = 10^(0:4),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = 10^(-6:1),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "bl", size = .25,
                      short = unit(1, "pt"), mid = unit(1.6, "pt"), long = unit(2.4, "pt")) +
  labs(title = "Observed vs sampling null",
       x = expression("Effective number of nuclei ("*italic(n)[eff]*")"),
       y = "Expression floor") +
  theme_pub()

ratio <- inner_join(obs, nul, by = c("donor", "ct"), suffix = c("_o", "_m")) %>%
  mutate(r = y_o / y_m) %>% filter(is.finite(r), r > 0)
med <- median(ratio$r)

p3b <- ggplot(ratio, aes(r)) +
  geom_histogram(bins = 46, fill = C_SEA, colour = NA, alpha = .85) +
  geom_vline(xintercept = 1, linetype = "22", linewidth = .4, colour = "grey20") +
  geom_vline(xintercept = med, colour = C_ULM, linewidth = .55) +
  annotate("text", x = med * 1.06, y = Inf, vjust = 1.45, hjust = 0, size = 2.5,
           colour = C_ULM, lineheight = .95,
           label = sprintf("median\n%.2f-fold", med)) +
  scale_y_continuous(expand = expansion(mult = c(0, .06))) +
  labs(title = "Excess over pure sampling",
       x = "Observed / null floor, per pair", y = "Pairs") +
  theme_pub()

byct <- byct %>%
  mutate(ct = factor(ct, levels = ct[order(b)]),
         wide = span >= quantile(span, .80))
med_b <- median(byct$b)

p3c <- ggplot(byct, aes(b, ct)) +
  geom_vline(xintercept = -0.5, linetype = "22", linewidth = .4, colour = "grey20") +
  geom_vline(xintercept = med_b, colour = C_ULM, linewidth = .5, alpha = .65) +
  geom_point(aes(size = wide, colour = wide)) +
  scale_size_manual(values = c(`FALSE` = 1.25, `TRUE` = 2.5)) +
  scale_colour_manual(values = c(`FALSE` = C_SEA, `TRUE` = C_ULM)) +
  annotate("text", x = -0.30, y = 2.2, hjust = 1, vjust = 0, size = 2.3,
           colour = "grey25", lineheight = 1,
           label = sprintf("median %.3f\nlarge = widest range in\neffective nucleus number", med_b)) +
  labs(title = "Not driven by abundance",
       x = "Within-cell-type exponent", y = NULL) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 6.2))

fig3 <- p3a + p3b + p3c + plot_layout(widths = c(1, 1, 1.05)) +
  plot_annotation(tag_levels = "a")
save_fig(fig3, OUT, "Fig3_null_and_stratified", height = 72)

cat(sprintf("Fig2: L2 n=%d, L1 n=%d\n", nrow(L2), nrow(L1)))
cat(sprintf("Fig3: pairs n=%d, ratio median %.3f (IQR %.2f-%.2f), cell types %d, median b %.3f\n",
            nrow(f3dat), med, quantile(ratio$r, .25), quantile(ratio$r, .75),
            nrow(byct), med_b))
cat("-> ", OUT, "\n")
