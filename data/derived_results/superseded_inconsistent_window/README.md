# Superseded: promoter-enrichment odds ratios computed with an inconsistent window

These four files are the original reference link-set odds ratios. They bounded the
**background** at ±1 Mb but did not bound the **detected link set**, so 1,770 links beyond
1 Mb (in `human_brain_3k`) entered the numerator while the denominator stopped at 1 Mb.
They are retained because they are what the first version of the manuscript reported and
because Supplementary Note 3 discusses the error.

**Do not cite these.** The current values are in `../window_500000/` (the window of the
primary analysis, and the one the manuscript quotes) and `../window_1000000/`. Regenerate
with `code/analysis/p5_18_download_arc_refs.sh` then `p5_19_recompute_refs_at_window.sh`.
