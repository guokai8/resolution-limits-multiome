# fig_common.R -- shared palette, theme and helpers for every figure script.
#
# Sourced by fig1_fig4_fig5, fig2_fig3, fig6, fig7 and supp_figs. It defines
# four things and nothing else:
#
#   1. the locale and font setup that the output devices need,
#   2. the colour constants, so a cohort keeps one colour across every panel,
#   3. theme_pub() / theme_blank() / ann() / save_fig(), the house style,
#   4. two statistics helpers, fit_line() and mh_repeats().
#
# House style follows the journal's figure requirements: no panel carries a
# sentence-style title. A panel shows only its letter, its axis labels and any
# annotation that has to sit next to the data. Everything conceptual belongs in
# the figure legend.

## ---- Locale ----------------------------------------------------------------
## Rscript inherits the caller's locale, and under the C locale R transliterates
## any character outside Latin-1 to ".." on its way to the graphics device. A
## kappa becomes two periods in both the PNG and the PDF, with no warning at any
## point. Pin LC_CTYPE to UTF-8 here so the scripts do not depend on how they
## were invoked, and assert it so a failure is loud rather than silent.
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

## ---- Font ------------------------------------------------------------------
## Arial first, Helvetica if Arial is missing, the system sans as a last resort.
## Either of the first two satisfies the journal.
BASE_FAMILY <- local({
  fams <- tryCatch(systemfonts::system_fonts()$family, error = function(e) character())
  if ("Arial" %in% fams) "Arial" else if ("Helvetica" %in% fams) "Helvetica" else "sans"
})

## ---- Colour ----------------------------------------------------------------
## Two rules run side by side.
##   (1) A cohort keeps its colour wherever it appears. Cohort identity wins
##       over everything else, so one mapping carries across the whole set.
##   (2) A conceptual panel that draws no cohort distinction (Figure 1c,
##       Figure 6) gives each of the three inference layers its own colour.
C_ULM  <- "#1F4E79"   # motor cortex, the primary cohort
C_SEA  <- "#C0504D"   # Seattle atlas
C_PSY  <- "#2E8B74"   # PsychAD split-aliquot replicates
C_EXT  <- "#D98A2B"   # NABEC/HBCC, the external cohort
C_NULL <- "#7F7F7F"   # sampling null and trans-pairing null
C_MID  <- "#5B9BD5"   # mid and light blue, for ordinal shading inside one panel
C_LITE <- "#9DC3E6"
C_GREY <- "#BBBBBB"   # non-replicated donors and other secondary series
COH    <- c("Motor cortex" = C_ULM, "Seattle atlas" = C_SEA)

## Rule (2): one colour per inference layer.
L_COMP <- "#2E8B74"   # composition
L_EXPR <- "#1F4E79"   # expression
L_REG  <- "#C0504D"   # regulation
LAYERC <- c("Composition" = L_COMP, "Expression" = L_EXPR, "Regulation" = L_REG)

## A theme's base_family only reaches theme elements -- axis text, titles, tags.
## Text drawn by a geom, including everything from annotate("text", ...), keeps
## the device default instead, which leaves a second font embedded in the PDF.
## Setting the geom defaults too means one font family per figure.
update_geom_defaults("text",  list(family = BASE_FAMILY))
update_geom_defaults("label", list(family = BASE_FAMILY))

## ---- Theme -----------------------------------------------------------------
## Titles and subtitles are blanked rather than merely left unset, so a stray
## labs(title = ...) in a figure script cannot put a title back into a panel.
theme_pub <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = BASE_FAMILY) +
    theme(
      axis.line           = element_line(linewidth = .35, colour = "grey20"),
      axis.ticks          = element_line(linewidth = .35, colour = "grey20"),
      axis.text           = element_text(colour = "grey20"),
      axis.title          = element_text(colour = "grey10"),
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

## An in-panel annotation. These carry what a title would otherwise have
## carried, so they are used freely. Defaults place the text top-left.
ann <- function(x, y, lab, size = 2.5, colour = "grey25", hjust = 0, vjust = 1, ...) {
  annotate("text", x = x, y = y, label = lab, size = size, colour = colour,
           hjust = hjust, vjust = vjust, family = BASE_FAMILY, lineheight = .95, ...)
}

## A bare canvas with no axes, for the schematic panels (Figure 1c, Figure 6a
## and 6b) that position every element themselves.
theme_blank <- function(base_size = 9) {
  theme_pub(base_size) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_blank())
}

## ---- Small helpers ---------------------------------------------------------
## A power-law reference line through (x0, y0) with slope b, drawn over `span`
## decades centred on x0. Used for the slope -1/2 sampling guides.
pow_guide <- function(x0, y0, b = -0.5, span = 1.6) {
  xs <- 10^(log10(x0) + c(-span, span) / 2)
  data.frame(x = xs, y = y0 * (xs / x0)^b)
}

## Ordinary least squares on log(y) ~ log(x), returned as 120 points ready to
## draw. The exponent it recovers is the b quoted in the text.
fit_line <- function(d, xr) {
  m  <- lm(log(y) ~ log(x), data = d)
  xs <- 10^seq(log10(xr[1]), log10(xr[2]), length.out = 120)
  data.frame(x = xs, y = exp(coef(m)[1]) * xs^coef(m)[2])
}

## Save once as vector PDF and once as 450 dpi PNG, 180 mm wide by default,
## which is the journal's two-column width. cairo_pdf rather than the base pdf()
## device: the base device cannot embed Arial and fails with "invalid font type".
save_fig <- function(p, out, name, width = 180, height = 68) {
  ggsave(file.path(out, paste0(name, ".pdf")), p,
         width = width, height = height, units = "mm", device = cairo_pdf)
  ggsave(file.path(out, paste0(name, ".png")), p,
         width = width, height = height, units = "mm", dpi = 450,
         device = ragg::agg_png)
  invisible(NULL)
}

## ---- Pooling the nucleus ladder --------------------------------------------
## Mantel-Haenszel pooling with the three downsampling seeds treated as repeated
## measures, matching the Methods and the B3 table from p5_17.
##
## Why the order matters: the seeds are Monte Carlo re-draws of the same nuclei,
## not independent strata. Pooling them as strata triples the apparent degrees
## of freedom -- the pseudoreplication this paper criticises in composition
## testing -- and narrows the intervals about twofold.
##
## The correct order is: pool cell types by Mantel-Haenszel *within* a seed,
## then combine across seeds with
##     variance = mean(within-seed sampling variance) + between-seed variance.
##
## Expected columns: prox_links / n_links, the promoter-proximal and total
## detected links; prox_tested_all / n_tested_all, the same counts over every
## tested pair, which is the background; plus `seed` and whatever `by` names.
mh_repeats <- function(df, by) {
  ## Within each seed and stratum: the Mantel-Haenszel R and S terms, plus the
  ## P and Q terms the Robins-Breslow-Greenland variance needs.
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

  ## Then across seeds. With one seed there is no between-seed term, so it
  ## contributes zero and the interval rests on the within-seed variance alone.
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
