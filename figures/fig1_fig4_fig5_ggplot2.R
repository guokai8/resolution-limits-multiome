#!/usr/bin/env Rscript
# Figures 1, 4 and 5. Same house style as fig2_fig3_ggplot2.R.
#
# Figure 1  The resource and the design.
#   a  Donor-by-chip layout; the 27 donors that appear on two chips sit above
#      the rule, and those pairs are the process replicates the paper rests on.
#   b  Three quality dimensions, replicated against non-replicated donors, each
#      divided by its own median so the three share one axis.
#   c  The three inference layers and the statistical unit each one is defined
#      on. This is the schematic that sets up the whole paper.
#
# Figure 4  What peak-gene inference returns at an attainable nucleus number.
#   a  The equalised substrate: 150 nuclei for each of 25 neuronal subtypes.
#   b  How many nuclei a detected peak is seen in -- the sparsity that drives
#      the result.
#   c  The observed promoter-enrichment odds ratio, significantly below one.
#   (The detection-power panel moved to Supplementary Figure 6.)
#
# Figure 5  The diagnostic's reference range.
#   a  Odds ratios across external link sets, with this study's operating point.
#   b  The same link set scored over three analysis windows.
#   c  Feature linkages against library size: no plateau.
#
# Usage: Rscript fig1_fig4_fig5_ggplot2.R [results_dir] [out_dir] [meta_dir]
#
# Unlike Figures 2, 3, 6 and 7 these panels need the primary cohort's donor
# metadata, which is not redistributed here.

args <- commandArgs(trailingOnly = TRUE)
RES  <- path.expand(if (length(args) >= 1) args[1] else "~/Desktop/P5_VulnerableEpigenome/results")
OUT  <- path.expand(if (length(args) >= 2) args[2] else "~/Desktop/ResearchD/figures")
META <- path.expand(if (length(args) >= 3) args[3] else
  "~/Desktop/ResearchD/Data/Other_Datasets/Multiome_Dataset/files")

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                  value = TRUE)[1])), "fig_common.R"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) read_csv(file.path(RES, f), show_col_types = FALSE)

## ======================================================================
## Figure 1
## ======================================================================
## quote = "" and comment.char = "" because donor metadata fields contain
## apostrophes and hash characters that would otherwise truncate rows.
meta <- read.delim(file.path(META, "Multiome_Dataset_Metadata.txt"),
                   colClasses = "character", quote = "", comment.char = "")
## A chip is the 10x run; wells within a chip share the loading and the library
## prep, so "same chip" and "different chip" are what separate the two kinds of
## replicate the paper distinguishes.
meta$Chip <- sub("Well[0-9]+$", "", meta$Well)

## Order donors by how many chips they appear on, then by total nuclei, so the
## replicated donors form a solid block at the top of the heatmap.
cc <- meta %>% count(ID, Chip, name = "n")
nchip <- cc %>% count(ID, name = "n_chip")
ord <- nchip %>%
  left_join(cc %>% group_by(ID) %>% summarise(tot = sum(n), .groups = "drop"), by = "ID") %>%
  arrange(desc(n_chip), desc(tot))
n_rep <- sum(ord$n_chip >= 2)

## Expand to the full donor x chip grid so that absent combinations are drawn
## as empty cells rather than omitted.
grid <- expand.grid(ID = ord$ID, Chip = sort(unique(meta$Chip)),
                    stringsAsFactors = FALSE) %>%
  left_join(cc, by = c("ID", "Chip")) %>%
  mutate(n = tidyr::replace_na(n, 0),
         ID = factor(ID, levels = rev(ord$ID)),
         ## Sort chips numerically: a plain sort puts Chip10 before Chip2.
         Chip = factor(Chip, levels = paste0("Chip", sort(as.integer(
           sub("Chip", "", unique(meta$Chip)))))))

p1a <- ggplot(grid, aes(Chip, ID, fill = log10(n + 1))) +
  geom_tile() +
  ## The rule separates replicated from non-replicated donors. 27 donors are on
  ## two chips; 26 of those clear the 50-nucleus floor on both and are the pairs
  ## every downstream estimate uses.
  geom_hline(yintercept = nrow(ord) - n_rep + 0.5, colour = C_SEA, linewidth = .6) +
  annotate("text", x = 6.5, y = nrow(ord) - n_rep + 1.6, size = 2.3, colour = C_SEA,
           lineheight = .95, vjust = 0,
           label = paste0("Donor-matched process replicates:\n", n_rep,
                          " donors on two chips (26 with \u226550 nuclei on both)")) +
  scale_fill_gradient(low = "#F2F7FC", high = C_ULM) +
  scale_x_discrete(labels = function(x) sub("Chip", "", x)) +
  labs(x = "Chip", y = "Donor (ordered)") +
  theme_pub() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 6.5))

mito <- meta %>%
  mutate(pmito = suppressWarnings(as.numeric(percent_mito))) %>%
  group_by(ID) %>% summarise(pmito = median(pmito, na.rm = TRUE), .groups = "drop")
don <- rd("p5_fig1_donor_table.csv") %>%
  left_join(nchip, by = "ID") %>% left_join(mito, by = "ID") %>%
  mutate(rep = n_chip >= 2)

## Panel b asks one question: are the replicated donors a biased subset? The
## answer is what matters, not the breadth of the QC panel, so the figure keeps
## the three dimensions that bear on it and drops mitochondrial fraction.
## Each dimension is divided by its own median so all three share one axis.
dims <- c(n_nuclei = "Nuclei", median_nCount_RNA = "RNA\ncomplexity",
          median_nCount_ATAC = "ATAC\ncomplexity")
qual <- lapply(names(dims), function(v) {
  x <- suppressWarnings(as.numeric(don[[v]]))
  data.frame(dim = dims[[v]], rep = don$rep, val = x / median(x, na.rm = TRUE))
}) %>% bind_rows() %>% filter(is.finite(val), !is.na(rep)) %>%
  mutate(dim = factor(dim, levels = unname(dims)))

## Two-sided Wilcoxon per dimension. The panel reports the range of the three
## P values rather than three separate numbers.
pv <- sapply(names(dims), function(v) {
  x <- suppressWarnings(as.numeric(don[[v]]))
  suppressWarnings(wilcox.test(x[don$rep], x[!don$rep])$p.value)
})

set.seed(0)   # the jitter is cosmetic, but fixing it keeps the panel reproducible
p1b <- ggplot(qual, aes(dim, val, colour = rep)) +
  geom_point(position = position_jitterdodge(jitter.width = .22, dodge.width = .55),
             size = .75, alpha = .75, stroke = 0) +
  scale_colour_manual(values = c(`TRUE` = C_ULM, `FALSE` = C_GREY)) +
  annotate("text", x = .55, y = Inf, hjust = 0, vjust = 1.25, size = 2.3,
           colour = "grey25", lineheight = 1,
           label = sprintf("blue = replicated\ngrey = not\nP = %.2f-%.2f",
                           min(pv), max(pv))) +
  ## Headroom, or the tallest point in the Nuclei column lands on the legend.
  scale_y_continuous(expand = expansion(mult = c(.05, .30))) +
  labs(x = NULL, y = "Value / median") +
  theme_pub()

## ---- 1c  The statistical-unit schematic ------------------------------------
## The message: one set of nuclei supports three kinds of claim, each defined on
## a different statistical unit, and therefore each with its own resolution.
## Drawn by hand on a blank canvas rather than built from data.
box <- data.frame(
  x    = c(1, 2.5, 4),
  lab  = c("Composition", "Expression", "Regulation"),
  unit = c("donor\nproportion", "donor\npseudobulk", "nucleus,\nfeature\ncovariance"),
  col  = c(L_COMP, L_EXPR, L_REG))

p1d <- ggplot(box) +
  ## The shared source at the top.
  annotate("label", x = 2.5, y = 4.42, label = "Single-nucleus multiome data",
           size = 2.3, family = BASE_FAMILY, fill = "grey95",
           label.size = 0, label.padding = unit(2.6, "pt"), colour = "grey15") +
  ## The branch: one stem, one crossbar, three arrows down.
  annotate("segment", x = 2.5, xend = 2.5, y = 4.18, yend = 3.88,
           colour = "grey55", linewidth = .4) +
  annotate("segment", x = 1, xend = 4, y = 3.88, yend = 3.88,
           colour = "grey55", linewidth = .4) +
  annotate("segment", x = box$x, xend = box$x, y = 3.88, yend = 3.52,
           colour = "grey55", linewidth = .4,
           arrow = arrow(length = unit(3, "pt"), type = "closed")) +
  ## The three claims, in the layer colours used throughout.
  geom_label(aes(x, 3.26, label = lab, fill = I(col)), colour = "white",
             size = 2.15, fontface = "bold", family = BASE_FAMILY,
             label.size = 0, label.padding = unit(2.2, "pt")) +
  ## Down to the statistical unit each claim is defined on.
  annotate("segment", x = box$x, xend = box$x, y = 3.0, yend = 2.62,
           colour = "grey55", linewidth = .4,
           arrow = arrow(length = unit(3, "pt"), type = "closed")) +
  geom_text(aes(x, 2.22, label = unit, colour = I(col)),
            size = 2.05, family = BASE_FAMILY, lineheight = .95) +
  ## And the conclusion the rest of the paper measures.
  annotate("segment", x = .5, xend = 4.5, y = 1.62, yend = 1.62,
           colour = "grey75", linewidth = .35) +
  annotate("text", x = 2.5, y = 1.3, size = 2.35, colour = "grey15",
           family = BASE_FAMILY, lineheight = 1.05,
           label = "Different statistical units\n\u2192 different resolution limits") +
  scale_x_continuous(limits = c(0.3, 4.7)) +
  scale_y_continuous(limits = c(0.95, 4.68)) +
  theme_blank()


fig1 <- p1a + p1b + p1d +
  plot_layout(widths = c(1.1, 0.9, 1.15)) +
  plot_annotation(tag_levels = "a")
save_fig(fig1, OUT, "Fig1_resource", height = 68)

## ======================================================================
## Figure 4
## ======================================================================
cells <- rd("A2_cells_seed0.csv")
fe    <- rd("p5_claimA_features.csv")

nsub <- cells %>% count(celltype, name = "n") %>% arrange(desc(n)) %>%
  mutate(i = row_number())

## Panel a: the substrate. Every bar is 150 by construction; the panel exists to
## show that the comparison is equalised, so the QC statements sit above the
## bars rather than inside them -- white text on the gaps between bars is
## unreadable at print size.
p4a <- ggplot(nsub, aes(i, n)) +
  geom_col(fill = C_ULM, width = .85) +
  geom_hline(yintercept = 150, linetype = "22", linewidth = .4, colour = "grey20") +
  annotate("text", x = nrow(nsub) / 2, y = 157, size = 2.05, colour = "grey25",
           family = BASE_FAMILY, label = "150 nuclei per subtype") +
  scale_y_continuous(limits = c(0, 228), expand = expansion(mult = c(0, .02))) +
  ann(x = .5, y = 226, hjust = 0, size = 2.0, colour = "grey30",
      lab = paste("Test: equalised nucleus design",
                  "(reference link sets at full depth)",
                  "depth ratio 1.0000, no nucleus below target", sep = "\n")) +
  labs(y = "Nuclei", x = paste0(nrow(nsub), " neuronal subtypes")) +
  theme_pub() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## Panel b: in how many of a subtype's 150 nuclei is a detected peak actually
## observed? Points are subtype medians, bars the interquartile range. This is
## the sparsity that makes single-nucleus correlation unreliable here.
pk <- rd("panel_b_nuclei_per_peak_by_subtype.csv") %>%
  arrange(median_nuclei_per_peak) %>% mutate(i = row_number())
med_pk <- median(pk$median_nuclei_per_peak)

p4b <- ggplot(pk, aes(i, median_nuclei_per_peak)) +
  geom_linerange(aes(ymin = q25, ymax = q75), colour = C_LITE, linewidth = 1.1) +
  geom_point(colour = C_ULM, size = 1.2) +
  geom_hline(yintercept = med_pk, colour = C_SEA, linewidth = .55) +
  annotate("text", x = 1, y = 7.2, hjust = 0, size = 2.1, colour = C_SEA,
           lineheight = 1,
           label = sprintf("median %.0f of 150 nuclei\nin every subtype", med_pk)) +
  annotate("text", x = nrow(pk), y = 1.15, hjust = 1, size = 2.0, colour = "grey35",
           label = sprintf("%.0f%% of peaks seen in one nucleus",
                           100 * median(pk$frac_peaks_in_1_nucleus))) +
  scale_y_continuous(limits = c(0.8, 8.2), breaks = c(1, 2, 4, 6, 8)) +
  ann(x = 1, y = 8.1, lab = "Sparse detection", size = 2.5, colour = "grey15") +
  labs(x = "Neuronal subtypes", y = "Nuclei per detected peak") +
  theme_pub() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## Panel c: the measured odds ratio. Written in as a literal because it is the
## single headline number of the layer, quoted identically in the abstract, the
## Results and Supplementary Table 11; it is produced by p5_08.
## Depletion, not absence: the links recovered carry less positional structure
## than the background, which is what sparsity predicts.
orv <- data.frame(or = 0.663, lo = 0.568, hi = 0.768)
p4c <- ggplot(orv) +
  geom_errorbar(aes(1, ymin = lo, ymax = hi), width = .06, linewidth = .5,
                colour = C_SEA) +
  geom_point(aes(1, or), colour = C_SEA, size = 2.6) +
  geom_hline(yintercept = 1, linetype = "22", linewidth = .45, colour = "grey20") +
  ## Two lines, not one: on one line the label runs into the panel border.
  annotate("text", x = 1.22, y = 0.663, hjust = .5, size = 2.3, colour = C_SEA,
           family = BASE_FAMILY, lineheight = 1,
           label = "0.663\n[0.568, 0.768]") +
  annotate("text", x = .72, y = 0.04, hjust = 0, size = 2.0, colour = "grey30",
           family = BASE_FAMILY, lineheight = 1,
           label = "reduced recovery,\nnot absence of regulation") +
  coord_flip(xlim = c(.6, 1.45), ylim = c(0, 1.3)) +
  labs(x = NULL, y = "Promoter-enrichment odds ratio") +
  theme_pub() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())


fig4 <- p4a + p4b + p4c +
  plot_layout(widths = c(1, 1.2, .95)) + plot_annotation(tag_levels = "a")
save_fig(fig4, OUT, "Fig4_regulatory_floor", height = 68)

## ======================================================================
## Figure 5
## ======================================================================
## The four external cellranger-arc link sets come straight from p5_15; this
## study's row is the Layer 3 main analysis and is written in.
##
## The window must match the main analysis. An earlier run at +/-1 Mb used a
## window that was not internally consistent and is not comparable with 0.663;
## it is kept under superseded_inconsistent_window/ so the record is complete,
## and this panel reads window_500000/ only.
W500 <- file.path(dirname(RES), basename(RES), "window_500000")
if (!dir.exists(W500)) W500 <- file.path("data/derived_results", "window_500000")
ext <- list.files(W500, "^p5_promoterOR_.*\\.csv$", full.names = TRUE) %>%
  lapply(read_csv, show_col_types = FALSE) %>% bind_rows() %>%
  transmute(lab = sub("_sorted", "", label), or = OR, lo, hi)
stopifnot(nrow(ext) == 4)   # fail loudly if the window directory is incomplete
or_sets <- bind_rows(
  ext %>% arrange(or),
  data.frame(lab = "This study (150 nuclei)", or = 0.663, lo = 0.568, hi = 0.768)) %>%
  mutate(lab = factor(lab, levels = rev(lab)), dep = or < 1)

p5a <- ggplot(or_sets, aes(or, lab, colour = dep)) +
  ## The span of the four reference link sets.
  annotate("rect", xmin = 2.26, xmax = 3.56, ymin = -Inf, ymax = Inf,
           fill = C_ULM, alpha = .07) +
  ## This study's operating point, highlighted as its own row.
  ## x is a log10 scale, where -Inf and Inf become NaN and the whole rectangle
  ## is dropped without any error, so the bounds have to be finite.
  annotate("rect", xmin = .451, xmax = 8.95, ymin = .42, ymax = 1.5,
           fill = C_SEA, alpha = .10) +
  geom_vline(xintercept = 1, linetype = "22", linewidth = .45, colour = "grey20") +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = .9) +
  geom_point(size = 1.5) +
  geom_text(aes(x = hi * 1.13, label = sprintf("%.2f", or)), hjust = 0, size = 2.2) +
  ann(x = .465, y = 1.42, lab = "operating point", size = 2.1, colour = C_SEA,
      hjust = 0, vjust = .5) +
  ann(x = 2.84, y = 5.42, lab = "reference range", size = 2.1, colour = C_ULM,
      hjust = .5, vjust = .5) +
  scale_colour_manual(values = c(`TRUE` = C_SEA, `FALSE` = C_ULM)) +
  scale_x_log10(limits = c(.45, 9), breaks = c(.5, 1, 2, 5),
                labels = c("0.5", "1", "2", "5")) +
  scale_y_discrete(expand = expansion(add = c(.6, .9))) +
  labs(y = NULL,
       x = "Promoter-enrichment odds ratio (log scale)") +
  theme_pub() +
  theme(axis.text.y = element_text(size = 6.5))

## Panel b: one link set (human_brain_3k) scored at three windows, recomputed by
## p5_15 --window so that the background and the detected set are tightened
## together (see analysis/p5_19_recompute_refs_at_window.sh). Widening the
## window raises the ratio almost twofold, which is why the window has to be
## reported with the number.
win <- data.frame(w = factor(c("\u00b1500 kb", "\u00b11 Mb", "\u00b11.7 Mb"),
                             levels = c("\u00b1500 kb", "\u00b11 Mb", "\u00b11.7 Mb")),
                  or = c(3.558, 4.054, 6.593))
p5b <- ggplot(win, aes(w, or, fill = w)) +
  geom_col(width = .55) +
  geom_text(aes(label = sprintf("%.2f", or)), vjust = -0.5, size = 2.6,
            colour = "grey15") +
  scale_fill_manual(values = c(C_ULM, C_MID, C_LITE)) +
  scale_y_continuous(limits = c(0, 7.6), expand = expansion(mult = c(0, .02))) +
  scale_x_discrete(labels = c("\u00b1500 kb" = "\u00b1500\nkb", "\u00b11 Mb" = "\u00b11\nMb",
                              "\u00b11.7 Mb" = "\u00b11.7\nMb")) +
  labs(x = NULL, y = "Odds ratio, same link set") +
  theme_pub() +
  theme(axis.text.x = element_text(size = 6.8))

## Panel c: link count against nuclei x fragments per nucleus across 28 external
## multiome libraries. b is fitted elsewhere and fixed here; the intercept is
## the median offset at that slope, which is robust to the handful of libraries
## far off the trend. A positive exponent with no plateau means these libraries
## are nowhere near saturating linkage discovery.
lib <- rd("seaad_L3_libraries.csv") %>%
  filter(!is.na(n_cells), !is.na(link), !is.na(atac_frag)) %>%
  mutate(x = n_cells * atac_frag)
b_sat <- 0.421
a_sat <- exp(median(log(lib$link) - b_sat * log(lib$x)))
sat <- data.frame(x = range(lib$x)) %>% mutate(y = a_sat * x^b_sat)

p5c <- ggplot(lib, aes(x, link)) +
  geom_point(colour = C_SEA, size = 1.1, alpha = .85, stroke = 0) +
  geom_line(data = sat, aes(x, y), colour = "grey15", linewidth = .7) +
  annotate("text", x = min(lib$x) * 1.15, y = max(lib$link) * .95, hjust = 0,
           size = 2.3, colour = "grey20",
           label = sprintf("b = +%.2f  (no plateau)", b_sat)) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "bl", size = .25, short = unit(1, "pt"),
                      mid = unit(1.6, "pt"), long = unit(2.4, "pt")) +
  labs(x = "Nuclei \u00d7 fragments per nucleus", y = "Feature linkages detected") +
  theme_pub()

fig5 <- p5a + p5b + p5c +
  plot_layout(widths = c(1.25, .7, 1)) + plot_annotation(tag_levels = "a")
save_fig(fig5, OUT, "Fig5_diagnostic", height = 62)

## Counts behind each figure, so a run that silently lost rows shows up here.
cat(sprintf("Fig1: %d donors, %d chips, %d replicated; P = %.2f-%.2f\n",
            nrow(ord), n_distinct(meta$Chip), n_rep, min(pv), max(pv)))
cat(sprintf("Fig4: %d subtypes x %d nuclei; redundancy median %.1f, means %.3f-%.3f\n",
            nrow(nsub), unique(nsub$n)[1], median(fe$redundancy),
            min(fe$redundancy_mean), max(fe$redundancy_mean)))
cat(sprintf("Fig5: %d link sets, %d SEA-AD libraries\n", nrow(or_sets), nrow(lib)))
cat("-> ", OUT, "\n")
