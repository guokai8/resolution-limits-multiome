# 五张主图共用的配色与主题（Fig 1-5 风格一致）
# 被 fig2_fig3_ggplot2.R 与 fig1_fig4_fig5_ggplot2.R source

suppressPackageStartupMessages({
  library(ggplot2)
library(tidyr); library(dplyr); library(readr)
  library(patchwork); library(scales); library(ragg)
})

## 全文一致的配色
C_ULM  <- "#1F4E79"   # 运动皮层队列 / 强调
C_SEA  <- "#C0504D"   # Seattle atlas / 反例
C_NULL <- "#7F7F7F"   # 抽样零模型
C_MID  <- "#5B9BD5"   # 中间层
C_LITE <- "#9DC3E6"   # 浅层
C_GREY <- "#BBBBBB"   # 非重复 / 次要
COH    <- c("Motor cortex" = C_ULM, "Seattle atlas" = C_SEA)

theme_pub <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line           = element_line(linewidth = .35, colour = "grey20"),
      axis.ticks          = element_line(linewidth = .35, colour = "grey20"),
      axis.text           = element_text(colour = "grey20"),
      axis.title          = element_text(colour = "grey10"),
      plot.title          = element_text(face = "plain", size = base_size, hjust = 0,
                                         margin = margin(b = 4, l = 13)),
      plot.title.position = "plot",
      plot.tag            = element_text(face = "bold", size = base_size + 2),
      plot.tag.position   = c(0, 1),
      axis.title.y        = element_text(margin = margin(r = 2)),
      legend.position     = "none",
      legend.background   = element_blank(),
      legend.key.height   = unit(9, "pt"),
      plot.margin         = margin(6, 7, 4, 3)
    )
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
