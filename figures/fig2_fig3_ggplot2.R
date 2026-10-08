#!/usr/bin/env Rscript
# Figures 2 and 3.
#
# Figure 2 is the paper's central visual argument: the same layout applied to
# two layers, giving two opposite answers.
#   a  Expression floors, two cohorts. The fitted lines are PARALLEL -- same
#      exponent, different intercept. Parallel, not coincident: the shape
#      transfers between datasets, the magnitude does not.
#   b  Composition discrepancies, same two cohorts. The fitted lines DIVERGE --
#      the exponents themselves differ.
#   c  Overdispersion in four replicate sets, ordered by what each replicate
#      pair spans.
# Panels a and b share one x range and both carry a slope -1/2 guide, so the
# comparison of slopes cannot be an artefact of panel width or axis span.
#
# Figure 3 asks whether the expression floor is just multinomial sampling.
#   a  Observed floors against a matched multinomial null.
#   b  Per-pair ratio of the two (median 1.30-fold).
#   c  Within-cell-type exponents, sized by how well each one is determined.
#
# Usage: Rscript fig2_fig3_ggplot2.R [results_dir] [out_dir] [derived_results_dir]
#
# results_dir holds the primary analysis outputs and is not redistributed here;
# derived_results_dir holds the deposited tables. Panel 2c reads only the
# latter, so it regenerates from this repository alone.

args    <- commandArgs(trailingOnly = TRUE)
RES     <- path.expand(if (length(args) >= 1) args[1] else "~/Desktop/P5_VulnerableEpigenome/results")
OUT     <- path.expand(if (length(args) >= 2) args[2] else "~/Desktop/ResearchD/figures")
DR_REPO <- if (length(args) >= 3) args[3] else "data/derived_results"

## Find fig_common.R next to this script, whatever directory Rscript was run in.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                  value = TRUE)[1])), "fig_common.R"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- Load ------------------------------------------------------------------
rd <- function(f) read_csv(file.path(RES, f), show_col_types = FALSE)

## L2 is the expression layer, L1 the composition layer.
##
## ALS26 is dropped for the same reason it is dropped everywhere else: its
## second chip yielded 37 nuclei, below the 50-nucleus floor the pairing rule
## requires. The rule is applied to the expression layer too, not only to
## composition, so the two layers rest on the same set of pairs.
als_L2 <- rd("p5_floor_scaling_pairs.csv") %>% filter(donor != "ALS26") %>%
  transmute(x = n_eff, y = median_abs_log2FC, cohort = "Motor cortex") %>%
  filter(x > 0, y > 0)
sea_L2 <- rd("seaad_L2_obs.csv") %>%
  transmute(x = n_eff, y = floor, cohort = "Seattle atlas") %>% filter(x > 0, y > 0)

als_L1 <- rd("p5_L1_scaling_pairs.csv") %>%
  transmute(x = n_eff, y = z, cohort = "Motor cortex") %>% filter(x > 0, y > 0)
sea_L1 <- rd("seaad_L1_pairs.csv") %>%
  transmute(x = n_eff, y = z, cohort = "Seattle atlas") %>% filter(x > 0, y > 0)

## Figure 3 needs the observed floors, the matched null, and the per-cell-type
## exponent fits.
obs <- rd("seaad_L2_obs.csv")  %>% transmute(donor, ct, x = n_eff, y = floor)
nul <- rd("seaad_L2_null.csv") %>% transmute(donor, ct, x = n_eff, y = floor)
byct <- rd("seaad_L2_by_celltype.csv")

## ======================================================================
## Figure 2
## ======================================================================
L2 <- bind_rows(als_L2, sea_L2)
L1 <- bind_rows(als_L1, sea_L1)
xr <- range(c(L2$x, L1$x))     # one x range for a and b, so slopes stay comparable

## Each cohort is fitted over its own x support, not over the shared range, so
## no line is extrapolated past the nuclei that cohort actually has.
f2a_fit <- bind_rows(
  fit_line(als_L2, range(als_L2$x)) %>% mutate(cohort = "Motor cortex"),
  fit_line(sea_L2, range(sea_L2$x)) %>% mutate(cohort = "Seattle atlas"))
f2b_fit <- bind_rows(
  fit_line(als_L1, range(als_L1$x)) %>% mutate(cohort = "Motor cortex"),
  fit_line(sea_L1, range(sea_L1$x)) %>% mutate(cohort = "Seattle atlas"))

## In-panel text legend: one coloured line per cohort, stacked downwards on a
## log axis by repeatedly multiplying by 0.55.
leg <- function(d, labs, xpos, ypos) {
  data.frame(cohort = names(labs), lab = unname(labs),
             x = xpos, y = ypos * (0.55^(seq_along(labs) - 1)))
}

## The exponents and intervals quoted here come from the bootstrap in p5_09 and
## p5_10 and are the values reported in the text.
f2a_leg <- leg(NULL, c("Motor cortex" = "Motor cortex  b = -0.497 [-0.567, -0.432]",
                       "Seattle atlas" = "Seattle atlas  b = -0.505 [-0.524, -0.483]"),
               xr[1] * 1.5, 0.075)
f2b_leg <- leg(NULL, c("Motor cortex" = "Motor cortex  b = -0.152 [-0.340, +0.111]",
                       "Seattle atlas" = "Seattle atlas  b = -0.584 [-0.838, -0.369]"),
               xr[1] * 1.5, 2.2e-4)

## Panel a: expression floors. Points are single donor-by-cell-type
## observations and are drawn faintly; the fitted lines carry the message.
p2a <- ggplot(L2, aes(x, y, colour = cohort)) +
  geom_point(alpha = .13, size = .5, stroke = 0, show.legend = FALSE) +
  geom_line(data = pow_guide(25, 2.6), aes(x, y), inherit.aes = FALSE,
            colour = "grey55", linetype = "22", linewidth = .4) +
  ## The guide's label sits below and right of the dashed line, clear of the
  ## two-line annotation in the upper right.
  annotate("text", x = 25 * 7, y = 2.6 * (7)^-0.5 * 0.80,
           label = "slope \u22121/2", size = 2.4, colour = "grey45",
           family = BASE_FAMILY) +
  geom_line(data = f2a_fit, aes(x, y, colour = cohort), linewidth = .85) +
  geom_text(data = f2a_leg, aes(x, y, label = lab, colour = cohort),
            hjust = 0, size = 2.35, show.legend = FALSE) +
  annotate("text", x = xr[2] * .92, y = 4.2, hjust = 1, vjust = 1, size = 2.5,
           colour = "grey25", lineheight = .95,
           label = "Same scaling exponent,\ndataset-specific magnitude") +
  scale_colour_manual(values = COH) +
  scale_x_log10(limits = xr, breaks = 10^(0:4),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = 10^(-6:1),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "bl", size = .25,
                      short = unit(1, "pt"), mid = unit(1.6, "pt"), long = unit(2.4, "pt")) +
  labs(x = expression("Effective number of nuclei ("*italic(n)[eff]*")"),
       y = "Expression floor") +
  theme_pub()

## Panel b: composition discrepancies, built exactly like panel a so that the
## difference between the panels is the data and nothing else.
p2b <- ggplot(L1, aes(x, y, colour = cohort)) +
  geom_point(alpha = .13, size = .5, stroke = 0, show.legend = FALSE) +
  geom_line(data = pow_guide(60, 3e-3, span = 1.5), aes(x, y), inherit.aes = FALSE,
            colour = "grey55", linetype = "22", linewidth = .4) +
  annotate("text", x = 60 * 3.2, y = 3e-3 * (3.2)^-0.5 * 2.0,
           label = "slope \u22121/2", size = 2.4, colour = "grey45",
           family = BASE_FAMILY) +
  geom_line(data = f2b_fit, aes(x, y, colour = cohort), linewidth = .85) +
  geom_text(data = f2b_leg, aes(x, y, label = lab, colour = cohort),
            hjust = 0, size = 2.35, show.legend = FALSE) +
  annotate("text", x = xr[2] * .92, y = 2.2, hjust = 1, vjust = 1, size = 2.5,
           colour = "grey25", lineheight = .95,
           label = "Divergent scaling") +
  scale_colour_manual(values = COH) +
  scale_x_log10(limits = xr, breaks = 10^(0:4),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = 10^(-6:1),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "bl", size = .25,
                      short = unit(1, "pt"), mid = unit(1.6, "pt"), long = unit(2.4, "pt")) +
  labs(x = expression("Effective number of nuclei ("*italic(n)[eff]*")"),
       y = "Composition discrepancy") +
  theme_pub()

## Panel c: four replicate sets, ordered by what the replicate pair spans --
## the more of the protocol a pair repeats, the more overdispersion it sees.
## kappa is recomputed by p5_17 from the deposited per-pair tables rather than
## written in as a literal, so the panel cannot drift from the tables.
lv <- c("Motor cortex", "Seattle atlas", "PsychAD MSSM", "PsychAD RADC")
od <- read_csv(file.path(DR_REPO, "p5_17_B7_cohort_overdispersion.csv"),
               show_col_types = FALSE) |>
  transmute(cohort = factor(cohort, levels = lv), od = kappa,
            grp = c("incl. tissue + dissociation", "incl. separate sample",
                    "loading / library only", "loading / library only")[
                      match(cohort, lv)])
## Bars are coloured by cohort, matching panels a and b. What each replicate
## spans is labelled directly instead -- inside the two tall bars, and with a
## bracket over the two short ones.
p2c <- ggplot(od, aes(cohort, od, fill = cohort)) +
  geom_col(width = .62) +
  geom_hline(yintercept = 1, linetype = "22", linewidth = .4, colour = "grey20") +
  geom_text(aes(label = sprintf("%.2f", od)), vjust = -0.6, size = 2.5,
            colour = "grey15") +
  annotate("text", x = 1, y = .25, angle = 90, hjust = 0, vjust = .5,
           size = 2.0, colour = "white", label = "incl. tissue + dissociation") +
  annotate("text", x = 2, y = .25, angle = 90, hjust = 0, vjust = .5,
           size = 2.0, colour = "white", label = "incl. separate sample") +
  annotate("segment", x = 2.62, xend = 4.38, y = 1.78, yend = 1.78,
           colour = "grey45", linewidth = .35) +
  annotate("text", x = 3.5, y = 1.83, vjust = 0, size = 2.0, colour = "grey35",
           label = "loading / library only") +
  scale_fill_manual(values = c(`Motor cortex` = C_ULM, `Seattle atlas` = C_SEA,
                               `PsychAD MSSM` = C_PSY, `PsychAD RADC` = C_EXT)) +
  scale_x_discrete(labels = c("Motor\ncortex", "Seattle\natlas", "PsychAD\nMSSM", "PsychAD\nRADC")) +
  scale_y_continuous(limits = c(0, 4.85), expand = expansion(mult = c(0, .02))) +
  labs(x = NULL,
       y = "Composition overdispersion") +
  theme_pub() +
  theme(axis.text.x = element_text(size = 6.4))

fig2 <- p2a + p2b + p2c + plot_layout(widths = c(1, 1, .95)) +
  plot_annotation(tag_levels = "a")
save_fig(fig2, OUT, "Fig2_transferability", height = 68)

## ======================================================================
## Figure 3
## ======================================================================
## The null resamples both members of each pair multinomially from their pooled
## profile at the observed depths, so it has the pair's own depth and gene set
## and differs from the observation only in having no process variation.
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

## Panel a: the two exponents nearly coincide, which is the point -- the shape
## is sampling. The offset between the lines is what panel b measures.
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
  ann(x = max(f3dat$x) * .85, y = 2.8,
      lab = "Observed \u2248 sampling expectation", hjust = 1, size = 2.5) +
  labs(x = expression("Effective number of nuclei ("*italic(n)[eff]*")"),
       y = "Expression floor") +
  theme_pub()

## Panel b: the same pairs, observed over null. The median is the calibration
## constant -- sampling gets the shape right and the magnitude wrong by 1.30x.
ratio <- inner_join(obs, nul, by = c("donor", "ct"), suffix = c("_o", "_m")) %>%
  mutate(r = y_o / y_m) %>% filter(is.finite(r), r > 0)
med <- median(ratio$r)

p3b <- ggplot(ratio, aes(r)) +
  geom_histogram(bins = 46, fill = C_SEA, colour = NA, alpha = .85) +
  geom_vline(xintercept = 1, linetype = "22", linewidth = .4, colour = "grey20") +
  geom_vline(xintercept = med, colour = C_ULM, linewidth = .55) +
  ## Put the median label inside the distribution, at about 60% of peak height,
  ## rather than on the top or bottom edge where it reads as a caption.
  annotate("label", x = med * 1.15, y = max(hist(ratio$r, breaks = 46, plot = FALSE)$counts) * 0.62,
           size = 2.5, hjust = 0,
           colour = C_ULM, fill = "white", label.size = 0, label.padding = unit(1.6, "pt"),
           family = BASE_FAMILY,
           label = sprintf("median %.2f-fold", med)) +
  scale_y_continuous(expand = expansion(mult = c(0, .06))) +
  labs(x = "Observed / null floor, per pair", y = "Pairs") +
  theme_pub()

## Panel c: one exponent per cell type, ordered by value. `span` is the range of
## effective nucleus number backing that fit, so the largest points are the
## best-determined exponents, not the most extreme ones.
byct <- byct %>%
  mutate(ct = factor(ct, levels = ct[order(b)]),
         wide = span >= quantile(span, .80))
med_b <- median(byct$b)

## Label only the three steepest and three flattest. Axis tick labels crowd
## adjacent rows together at this panel height, so the six names are drawn next
## to their own points; the complete list is in the supplementary tables.
lab6 <- byct %>%
  filter(ct %in% levels(ct)[c(1:3, (nlevels(ct) - 2):nlevels(ct))])

p3c <- ggplot(byct, aes(b, ct)) +
  geom_vline(xintercept = -0.5, linetype = "22", linewidth = .4, colour = "grey20") +
  geom_vline(xintercept = med_b, colour = C_ULM, linewidth = .5, alpha = .65) +
  geom_point(aes(size = wide, colour = wide)) +
  scale_size_manual(values = c(`FALSE` = 1.25, `TRUE` = 2.5)) +
  scale_colour_manual(values = c(`FALSE` = C_SEA, `TRUE` = C_ULM)) +
  geom_text(data = lab6, aes(b, ct, label = ct), inherit.aes = FALSE,
            hjust = 1, nudge_x = -0.013, size = 2.0, colour = "grey25",
            family = BASE_FAMILY) +
  annotate("text", x = -0.30, y = 2.2, hjust = 1, vjust = 0, size = 2.3,
           colour = "grey25", lineheight = 1,
           label = sprintf("median %.3f\nlarge = widest range in\neffective nucleus number", med_b)) +
  ## Extra room on the left so the direct labels are not clipped.
  scale_x_continuous(expand = expansion(mult = c(.26, .05))) +
  labs(x = "Within-cell-type exponent", y = NULL) +
  theme_pub() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

fig3 <- p3a + p3b + p3c + plot_layout(widths = c(1, 1, 1.05)) +
  plot_annotation(tag_levels = "a")
save_fig(fig3, OUT, "Fig3_null_and_stratified", height = 72)

## Echo the counts behind each panel, so a run that silently lost rows is
## visible in the console rather than only in the figure.
cat(sprintf("Fig2: L2 n=%d, L1 n=%d\n", nrow(L2), nrow(L1)))
cat(sprintf("Fig3: pairs n=%d, ratio median %.3f (IQR %.2f-%.2f), cell types %d, median b %.3f\n",
            nrow(f3dat), med, quantile(ratio$r, .25), quantile(ratio$r, .75),
            nrow(byct), med_b))
cat("-> ", OUT, "\n")
