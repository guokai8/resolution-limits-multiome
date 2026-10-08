# Supplementary Figures 1-6.
#
# Usage: Rscript supp_figs_ggplot2.R <derived_results_dir> <out_dir>
#
# These are the sensitivity analyses and per-stratum breakdowns that the main
# figures summarise. All six read only the deposited tables, so they regenerate
# from this repository alone.
args <- commandArgs(trailingOnly = TRUE)
DR  <- if (length(args) >= 1) args[1] else "data/derived_results"
OUT <- if (length(args) >= 2) args[2] else "figures"
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
      value = TRUE)[1])), "fig_common.R"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## These figures are single panels, so they are printed to a device directly
## rather than going through ggsave and save_fig. cairo_pdf, not the base pdf()
## device: the base device cannot embed Arial and fails with "invalid font type".
save2 <- function(pl, name, w, h) {
  agg_png(file.path(OUT, paste0(name, ".png")), width = w, height = h,
          units = "mm", res = 450); print(pl); invisible(dev.off())
  cairo_pdf(file.path(OUT, paste0(name, ".pdf")), width = w/25.4, height = h/25.4)
  print(pl); invisible(dev.off())
}

## ---- SF1  The nucleus ladder, split by cell type ---------------------------
## Figure 7a and Table 3 report the ladder pooled over cell types. This splits
## it, to show the pooled curve is not hiding heterogeneity between types.
## mh_repeats() pools the seeds as repeated measures; see fig_common.R for why
## that matters.
pri <- read_csv(file.path(DR, "primary_L3_nucleus_ladder.csv"), show_col_types = FALSE) |>
  mh_repeats(c("celltype", "n")) |> mutate(cohort = "Primary (motor cortex)")
ext <- read_csv(file.path(DR, "nabec_L3_nucleus_ladder.csv"), show_col_types = FALSE) |>
  mh_repeats(c("celltype", "n")) |> mutate(cohort = "External (prefrontal cortex)")
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
  labs(x = "Nuclei per cell type", y = "Promoter-enrichment odds ratio") +
  theme_pub() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 6.5),
        strip.background = element_blank(), strip.text = element_text(size = 8.5))
save2(sf1, "SuppFig1_ladder_by_celltype", 180, 78)

## ---- SF2  Raising ATAC depth alone -----------------------------------------
## Does more ATAC depth rescue the equalised 150-nucleus design? It does not:
## the ratio is flat and below one across an eightfold range.
##
## The pooled values come straight from p5_17's B8 table rather than being
## recomputed here. The external replication table encodes the seed in the file
## name and has no seed column, so recomputing in place would silently treat the
## seeds as independent strata and understate the intervals.
dd <- read_csv(file.path(DR, "p5_17_B8_depth_series.csv"), show_col_types = FALSE) |>
  transmute(fold, OR = OR_seeds_as_repeats,
            lo = lo_seeds_as_repeats, hi = hi_seeds_as_repeats)
sf2 <- ggplot(dd, aes(fold, OR)) +
  geom_hline(yintercept = 1, linetype = "22", linewidth = .35, colour = "grey45") +
  geom_line(colour = C_ULM, linewidth = .45) +
  geom_linerange(aes(ymin = lo, ymax = hi), colour = C_ULM, linewidth = .4) +
  geom_point(size = 1.6, colour = C_ULM) +
  scale_x_log10(breaks = c(1, 2, 4, 8), labels = c("1x", "2x", "4x", "8x")) +
  ## The intervals run from 0.265 to 1.405 once the seeds are pooled as repeats.
  ## The limits have to cover that, or ggplot drops the ends of the bars with no
  ## warning and the figure quietly understates the uncertainty.
  scale_y_log10(limits = c(.24, 1.6), breaks = c(.25, .5, 1, 1.5)) +
  ann(x = 1.05, y = 1.5, size = 2.4, colour = "grey15", hjust = 0,
      lab = "ATAC depth alone does not\nrecover regulatory structure") +
  labs(x = "ATAC depth, relative to the equalised design\n(5,564 fragments per nucleus)",
       y = "Promoter-enrichment odds ratio") +
  theme_pub()
save2(sf2, "SuppFig2_depth_series", 110, 72)

## ---- SF3  The background and multiple-testing convention -------------------
## The odds ratio depends on two conventions that have to be reported with it:
## which pairs form the background, and whether peaks that cannot be called at
## the depth analysed enter the Benjamini-Hochberg denominator.
##
## An earlier version compared the two levels of `convention`, which was wrong:
## those two levels differ in ATAC depth by 5 to 22 fold (unmatched 27,973 to
## 120,798 against p5_08 matched 5,564), so that comparison measured depth, not
## convention. The convention contrast lives *within* a row at fixed depth --
## OR, which counts non-callable peaks in the background and the denominator,
## against OR_varpeaks, which excludes them. Keep only the fixed-depth group.
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
  labs(x = NULL, y = "Promoter-enrichment odds ratio\n(same nuclei, same depth)") +
  theme_pub() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 6.2), legend.key.height = unit(13, "pt"))
save2(sf3, "SuppFig3_convention_sensitivity", 120, 82)

## ---- SF4  External replication across eleven strata ------------------------
## Per-stratum odds ratios, so these are read straight from the deposited table
## without pooling across seeds. The `da[0-9]+` files are an earlier external
## run whose pipeline was not aligned with p5_08; it is excluded here and the
## mismatch is documented in protocol/.
ext_rep <- read_csv(file.path(DR, "nabec_L3_external_replication.csv"),
                    show_col_types = FALSE)
fo <- ext_rep |> filter(n == 150, !grepl("da[0-9]+", file)) |>
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
  ## Headroom above the top stratum for the band's label.
  scale_x_discrete(expand = expansion(add = c(.6, 1.75))) +
  ## The shaded band is the primary cohort's own interval, 0.663 [0.568, 0.768].
  ## Say so inside the panel: an unlabelled band is unreadable on its own.
  ann(x = nrow(fo) + .85, y = .663, lab = "primary cohort,\n150 nuclei (95% CI)",
      size = 2.0, colour = "grey30", hjust = .5, vjust = .5) +
  labs(x = NULL, y = "Promoter-enrichment odds ratio at 150 nuclei") +
  theme_pub() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 6.5))
save2(sf4, "SuppFig4_external_replication_forest", 110, 90)

## ---- SF5  190 composition contrasts against the multinomial floor ----------
## Each contrast's donor-level difference divided by the multinomial floor of
## its own design. Left of 1 means the difference is smaller than pure sampling
## noise. Colour records which test calls it significant: pooling cells ignores
## the donor as the statistical unit and calls many contrasts that the
## donor-level test does not.
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
  ## Headroom so the tallest bar is not clipped by the panel border.
  scale_y_continuous(expand = expansion(mult = c(0, .14))) +
  labs(x = "Donor-level difference / multinomial expectation",
       y = "Contrasts", fill = NULL) +
  ## Legend top-left, which is the empty corner here; top-right covers the peak.
  theme_pub() +
  theme(legend.position = c(.02, .98), legend.justification = c(0, 1),
        legend.text = element_text(size = 6.5))
save2(sf5, "SuppFig5_floor_reassessment", 110, 74)

cat("Supplementary Figures 1-6 written\n")

## ---- SF6  Detection power for the equalised design -------------------------
## Moved out of Figure 4 as panel d. It answers the obvious objection to the
## Layer 3 result: was the link set simply too small to show enrichment?
##
## p0 is the promoter-proximal rate in the background, n_test the tested pairs,
## n_obs the links actually recovered. For a true odds ratio `orr`, p1 is the
## implied rate among detected links; se combines the binomial error of the
## detected set with that of the background; z is the resulting test statistic
## as a function of how many links were recovered.
p0     <- 0.0095
n_test <- 3.03e6
n_obs  <- 33227
pw <- lapply(c(3.0, 2.0, 1.5), function(orr) {
  nl <- 10^seq(2, 5, length.out = 160)
  p1 <- orr * p0 / (1 - p0 + orr * p0)
  se <- sqrt(p1 * (1 - p1) / nl + p0 * (1 - p0) / n_test)
  data.frame(nl = nl, z = (p1 - p0) / se, or = factor(orr, levels = c(3, 2, 1.5)))
}) %>% bind_rows()
## Put each curve's label at a chosen height rather than at its end, so the
## three labels do not collide in the upper right.
ztar <- c(`3` = 9.0, `2` = 6.2, `1.5` = 3.6)
tip <- pw %>% group_by(or) %>%
  slice_min(abs(z - ztar[as.character(or)]), n = 1) %>% ungroup()

sf6 <- ggplot(pw, aes(nl, z, colour = or)) +
  geom_line(linewidth = .8) +
  geom_hline(yintercept = 1.96, linetype = "22", linewidth = .4, colour = "grey20") +
  geom_vline(xintercept = n_obs, colour = C_SEA, linewidth = .8) +
  geom_text(data = tip, aes(nl, z, label = paste("OR", or), colour = or),
            hjust = 1.15, vjust = -0.35, size = 2.0, show.legend = FALSE) +
  annotate("text", x = 120, y = 2.4, hjust = 0, size = 2.0, colour = "grey30",
           label = "P = 0.05") +
  annotate("text", x = n_obs * 1.25, y = 1.0, hjust = 0, size = 2.0, colour = C_SEA,
           lineheight = 1, label = "33,227\nlinks") +
  scale_colour_manual(values = c(`3` = C_LITE, `2` = C_MID, `1.5` = C_ULM)) +
  scale_x_log10(breaks = 10^(2:5), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(expand = expansion(mult = c(0, .02))) +
  ## All three curves run far above the panel; clip rather than compress, since
  ## the question is only where they cross P = 0.05.
  coord_cartesian(ylim = c(0, 11)) +
  labs(x = "Detected links", y = "Detection z-score") +
  theme_pub()

save2(sf6, "SuppFig6_detection_power", 95, 72)
