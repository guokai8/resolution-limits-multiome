#!/usr/bin/env Rscript
# Figure 7 -- does the nucleus dependence survive a change of linking procedure?
#
# Usage: Rscript fig7_ggplot2.R <supplementary_tables_dir> <out_dir>
#
# The reviewer's objection this figure answers: single-nucleus correlation is a
# weak linker, so the nucleus dependence reported for Layer 3 might be an
# artefact of that choice rather than a property of the data. The test holds the
# nuclei, the depth and the tested gene-peak set fixed and varies only the
# procedure.
#
#   a  Promoter enrichment against nucleus number, three procedures x two cell
#      types. The dependence is present in every arm.
#   b  The same runs scored against the trans-pairing null, which asks whether a
#      link set carries positional information at all. Only series with at least
#      200 links are shown; below that the null is uninformative.
#   c  Aggregate size k at fixed data and 2,400 nuclei -- one analyst-chosen
#      hyperparameter moves the link count tenfold.
#   d  Stability of the trans-pairing null over nine independent reassignments
#      and a leave-one-chromosome-out jackknife.
#
# Everything is read from the deposited supplementary tables, so this figure
# regenerates from the repository alone.
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

## The three arms. "ArchR-like retrieval" means ArchR's own published defaults
## (corCutOff 0.45, FDRCutOff 1e-4, varCutOff 0.25), not merely aggregation; the
## FDR-only arm keeps the aggregation but drops those thresholds, which
## separates the effect of aggregating from the effect of the thresholds.
arm_of <- function(linker, archr) {
  ifelse(linker == "single_cell", "Single-nucleus correlation",
         ifelse(archr, "Aggregation + ArchR-like retrieval", "Aggregation + FDR only"))
}
ARMC <- c("Single-nucleus correlation"              = C_ULM,
          "Aggregation + ArchR-like retrieval"      = C_SEA,
          "Aggregation + FDR only"                  = C_GREY)
CT <- c("Exc_LINC00507_FREM3" = 16, "Oligodendrocytes" = 17)

## k = 0 is the unaggregated arm; k = 25 is the aggregate size used for the main
## comparison. The k sweep in panel c uses the other values.
main <- t20 |>
  filter(window == 500000, k %in% c(0, 25)) |>
  mutate(arm = arm_of(linker, archr_defaults),
         arm = factor(arm, levels = names(ARMC)),
         celltype = factor(celltype, levels = names(CT)))

## ---- a  Promoter enrichment against nucleus number -------------------------
## The vertical rule at 150 marks the equalised design of Figure 4, so the
## reader can see where that result sits on this ladder.
pa <- ggplot(main, aes(n, promoter_OR, colour = arm, shape = celltype)) +
  geom_hline(yintercept = 1, linetype = "22", linewidth = .35, colour = "grey45") +
  geom_vline(xintercept = 150, linewidth = .3, colour = "grey55") +
  geom_line(aes(group = interaction(arm, celltype)), linewidth = .45, alpha = .9) +
  geom_point(size = 1.5) +
  ## limits + drop = FALSE so the legend lists all three arms and both cell
  ## types even where a panel happens not to draw one of them.
  scale_colour_manual(values = ARMC, name = NULL,
                      limits = names(ARMC), drop = FALSE) +
  scale_shape_manual(values = CT, name = NULL,
                     limits = names(CT), drop = FALSE) +
  scale_x_log10(breaks = c(150, 400, 900, 2400, 7500), labels = label_comma(accuracy = 1)) +
  scale_y_log10(breaks = c(.5, 1, 2, 5, 10, 30)) +
  ann(x = 165, y = 34, lab = "Same nuclei,\ndifferent procedures", size = 2.4,
      colour = "grey15") +
  labs(x = "Nuclei per cell type", y = "Promoter-enrichment odds ratio") +
  ## Only this panel emits a legend. patchwork's guides = "collect" cannot merge
  ## guides whose key glyphs differ, so letting every panel emit one produces
  ## duplicates at the bottom of the figure.
  guides(colour = guide_legend(nrow = 1, order = 1),
         shape  = guide_legend(nrow = 1, order = 2)) +
  theme_pub()

## ---- b  Positional odds ratio against the trans-pairing null ---------------
## Each gene is tested against the window of a gene on a different chromosome,
## which is Signac's own published background logic. The ratio of observed to
## trans enrichment asks whether a detected set carries positional information
## beyond what the windows alone would give.
##
## The 200-link floor is not cosmetic: the null's own validation (panel d and
## Supplementary Note 5) shows it becomes unstable below roughly that count, so
## series under it are dropped rather than plotted and caveated.
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

## ---- c  Sensitivity to aggregate size --------------------------------------
## Same data, same 2,400 nuclei, same thresholds; only k changes. Plotting the
## odds ratio against the link count rather than against k shows the trade
## directly: more links, weaker enrichment.
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
  expand_limits(x = 6000, y = 26) +   # room for the k labels at both ends
  labs(x = "Links recovered", y = "Promoter-enrichment odds ratio") +
  theme_pub()

## ---- d  Is the trans-pairing null itself stable? ---------------------------
## Point = median over reassignments, bar = full range, so the bar is a spread
## over nulls rather than a confidence interval. Table 21 records the arm in the
## analysis's own vocabulary, so map it onto the labels used in the figure.
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

## Two rows of two. tag_levels restores the panel letters, which the shared
## legend would otherwise suppress along with the per-panel subtitles.
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
