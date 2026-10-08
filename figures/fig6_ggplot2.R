#!/usr/bin/env Rscript
# Figure 6 -- the framework figure: resolution follows from the statistical unit
# the claim is defined on.
#
# Usage: Rscript fig6_ggplot2.R <derived_results_dir> <out_dir>
#
# This is the conceptual centre of the paper rather than a results figure, so it
# carries no curves at all. Three panels:
#   a  the chain, biological question -> statistical unit -> resolution
#      estimator -> experimental requirement,
#   b  that chain filled in for each of the three classes of claim,
#   c  the three resulting requirements placed on one shared axis of nuclei per
#      cell type, which is the quantity an experimenter actually plans.
#
# Panels a and b are drawn by hand on a blank canvas; every number in panel c is
# quoted from the Results rather than recomputed here. The first argument is
# accepted for symmetry with the other figure scripts but is unused.
args <- commandArgs(trailingOnly = TRUE)
OUT  <- ifelse(length(args) >= 2, args[2], "figures")
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])),
                 "fig_common.R"))

## ---- a  The chain ----------------------------------------------------------
## Four boxes stacked top to bottom, each with an italic example underneath it.
## The example follows one claim -- a change in cell-type proportion -- the whole
## way down, so the reader sees a single path rather than four abstractions.
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
  ## Arrows between consecutive boxes: from just under box i to just above
  ## box i+1, hence y[-4] and y[-1].
  annotate("segment", x = 1, xend = 1, y = step$y[-4] - 0.46, yend = step$y[-1] + 0.2,
           colour = "grey55", linewidth = .45,
           arrow = arrow(length = unit(3.2, "pt"), type = "closed")) +
  scale_x_continuous(limits = c(0.3, 1.7)) +
  scale_y_continuous(limits = c(0.55, 4.45)) +
  labs(tag = "a") +
  theme_blank()

## ---- b  The same chain for all three claims --------------------------------
## A three-column table drawn with text geoms. Column positions are tuned so
## the widest entry in each column clears the next column at print size.
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
  ## The claim column is coloured by layer, matching Figure 1c.
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

## ---- c  The three requirements on one axis ---------------------------------
## Band = the range given in the Results, point = the central estimate. These
## are ranges, not confidence intervals, and the legend says so.
##
##   Expression  17-91     the measured support of the floor (median effective
##                         n = 39), so larger n is extrapolation, not evidence
##   Composition 44,435-   nuclei per group to resolve one percentage point at
##               499,575   the measured overdispersion, median to most
##                         overdispersed cell type; cost scales as kappa squared
##   Regulation  400-5,000 where promoter enrichment crosses unity, primary
##                         cohort to external cohort
##
## Three decades separate them. That gap is the paper's practical point, and it
## only shows up once all three sit on one axis.
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
  ## Labels sit at the geometric centre of each band, which is the midpoint on
  ## a log axis. hjust differs per row to keep them inside the panel.
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
