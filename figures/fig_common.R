# 七张主图与五张补充图共用的配色与主题。
# 按 Genome Biology 的图件要求：图内不放句子式标题，只留面板字母、轴标签与
# 必要的就地注释；概念信息交给图注。字体统一 Arial/Helvetica。

## Rscript 在 C locale 下会把非 Latin-1 字符（如 \u03ba）静默转写成 ".."，
## PNG 和 PDF 都会出错且不报任何警告。在这里锁定 UTF-8，脚本不依赖调用方环境。
if (!isTRUE(l10n_info()$`UTF-8`)) {
  for (lc in c("en_US.UTF-8", "C.UTF-8", "UTF-8")) {
    if (suppressWarnings(Sys.setlocale("LC_CTYPE", lc)) != "") break
  }
}
stopifnot(isTRUE(l10n_info()$`UTF-8`))

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr); library(dplyr); library(readr)
  library(patchwork); library(scales); library(ragg)
})

## 字体：Arial 优先，缺失时退回 Helvetica，再退回系统无衬线
BASE_FAMILY <- local({
  fams <- tryCatch(systemfonts::system_fonts()$family, error = function(e) character())
  if ("Arial" %in% fams) "Arial" else if ("Helvetica" %in% fams) "Helvetica" else "sans"
})

## ---- 配色 -------------------------------------------------------------------
## 两条规则并存：
##  (1) 同一队列在所有图里同色（队列身份优先）；
##  (2) 在没有队列区分的概念性面板里，三层各用一个色族。
C_ULM  <- "#1F4E79"   # 运动皮层队列（原队列）
C_SEA  <- "#C0504D"   # Seattle atlas
C_PSY  <- "#2E8B74"   # PsychAD split-aliquot
C_EXT  <- "#D98A2B"   # NABEC/HBCC 外部队列
C_NULL <- "#7F7F7F"   # 抽样零模型 / trans 零假设
C_MID  <- "#5B9BD5"
C_LITE <- "#9DC3E6"
C_GREY <- "#BBBBBB"   # 非重复 / 次要
COH    <- c("Motor cortex" = C_ULM, "Seattle atlas" = C_SEA)

## 三个推断层的色族（概念图与框架图用）
L_COMP <- "#2E8B74"   # 组成层
L_EXPR <- "#1F4E79"   # 表达层
L_REG  <- "#C0504D"   # 调控层
LAYERC <- c("Composition" = L_COMP, "Expression" = L_EXPR, "Regulation" = L_REG)

# theme 的 base_family 只管主题元素；geom_text / annotate("text") 仍走设备默认字体，
# 会在 PDF 里混进 Bitstream Vera。把 geom 默认字体也设成 Arial，全图只嵌一种字体。
update_geom_defaults("text",  list(family = BASE_FAMILY))
update_geom_defaults("label", list(family = BASE_FAMILY))

theme_pub <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = BASE_FAMILY) +
    theme(
      axis.line           = element_line(linewidth = .35, colour = "grey20"),
      axis.ticks          = element_line(linewidth = .35, colour = "grey20"),
      axis.text           = element_text(colour = "grey20"),
      axis.title          = element_text(colour = "grey10"),
      ## 图内不放标题：GB 要求概念信息由图注承担
      plot.title          = element_blank(),
      plot.subtitle       = element_blank(),
      plot.tag            = element_text(face = "bold", size = base_size + 2,
                                         family = BASE_FAMILY),
      plot.tag.position   = c(0, 1),
      axis.title.y        = element_text(margin = margin(r = 2)),
      legend.position     = "none",
      legend.background   = element_blank(),
      legend.key.height   = unit(9, "pt"),
      plot.margin         = margin(9, 7, 4, 3)
    )
}

## 面板内的就地注释（替代被删掉的标题）
ann <- function(x, y, lab, size = 2.5, colour = "grey25", hjust = 0, vjust = 1, ...) {
  annotate("text", x = x, y = y, label = lab, size = size, colour = colour,
           hjust = hjust, vjust = vjust, family = BASE_FAMILY, lineheight = .95, ...)
}

# 无坐标轴的画布（用于示意图 / 文字表）
theme_blank <- function(base_size = 9) {
  theme_pub(base_size) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_blank())
}

## 幂律参考线：过 (x0,y0)、斜率 b
pow_guide <- function(x0, y0, b = -0.5, span = 1.6) {
  xs <- 10^(log10(x0) + c(-span, span) / 2)
  data.frame(x = xs, y = y0 * (xs / x0)^b)
}

## log-log OLS 拟合线
fit_line <- function(d, xr) {
  m  <- lm(log(y) ~ log(x), data = d)
  xs <- 10^seq(log10(xr[1]), log10(xr[2]), length.out = 120)
  data.frame(x = xs, y = exp(coef(m)[1]) * xs^coef(m)[2])
}

## 统一的保存：矢量 PDF + 450 dpi PNG
save_fig <- function(p, out, name, width = 180, height = 68) {
  ggsave(file.path(out, paste0(name, ".pdf")), p,
         width = width, height = height, units = "mm", device = cairo_pdf)
  ggsave(file.path(out, paste0(name, ".png")), p,
         width = width, height = height, units = "mm", dpi = 450,
         device = ragg::agg_png)
  invisible(NULL)
}

## ---- 种子按重复计算抽样合并（与 Methods 及 p5_17 的 B3 一致）----------------
## 三个降采样种子是同一批核的重抽样，不是独立 strata。正确做法是：先在每个
## 种子内部按细胞类型做 Mantel-Haenszel 合并，再跨种子合成，方差为
## 「种子内抽样方差的均值 + 种子间下采样方差」。把种子当独立层会把自由度
## 虚增三倍，正是本文在组成层检验里批评的那种 pseudoreplication。
mh_repeats <- function(df, by) {
  per <- df |>
    mutate(a = prox_links, b = n_links - prox_links,
           cc = prox_tested_all, d = n_tested_all - prox_tested_all,
           tot = a + b + cc + d,
           R = a * d / tot, S = b * cc / tot,
           P = (a + d) / tot, Q = (b + cc) / tot) |>
    group_by(across(all_of(c(by, "seed")))) |>
    summarise(Rs = sum(R), Ss = sum(S), PR = sum(P * R), PS = sum(P * S),
              QR = sum(Q * R), QS = sum(Q * S), .groups = "drop") |>
    mutate(logOR = log(Rs / Ss),
           var_w = PR / (2 * Rs^2) + (PS + QR) / (2 * Rs * Ss) + QS / (2 * Ss^2))
  per |>
    group_by(across(all_of(by))) |>
    summarise(mu = mean(logOR),
              within = mean(var_w),
              between = if (dplyr::n() > 1) stats::var(logOR) else 0,
              n_seeds = dplyr::n(), .groups = "drop") |>
    mutate(se = sqrt(within + between),
           OR = exp(mu),
           lo = exp(mu - 1.96 * se), hi = exp(mu + 1.96 * se))
}
