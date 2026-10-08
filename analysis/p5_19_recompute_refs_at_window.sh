#!/usr/bin/env bash
# =============================================================================
# Step 19: recompute the four reference link sets at a stated window.
# =============================================================================
# Why this exists. The main analysis (p5_08_claimA_features.py) uses a 500 kb
# window, but the four reference sets were originally computed at 1 Mb. Placing
# 0.663 next to 2.45-4.01 therefore compared numbers from different windows,
# and the odds ratio depends on the window by nearly twofold. This script
# recomputes the references at whatever window is asked for, so the comparison
# is like for like.
#
# The superseded 1 Mb outputs are kept under
# data/derived_results/superseded_inconsistent_window/ rather than deleted, so
# the record of what was originally reported stays intact.
#
# Usage:
#   bash code/analysis/p5_19_recompute_refs_at_window.sh 500000
#   bash code/analysis/p5_19_recompute_refs_at_window.sh 1000000   # reproduce the old values
# =============================================================================
set -euo pipefail
WIN="${1:-500000}"
ARC="${ARC:-refs/arc}"
WORK="${WORK:-refs/arc/extracted}"
TSS="${TSS:-data/derived_results/p5_tss_gencode_v32.csv}"
OUT="${OUT:-data/derived_results/window_${WIN}}"
mkdir -p "$WORK" "$OUT"

for tgz in "$ARC"/*_analysis.tar.gz; do
  s=$(basename "$tgz" _analysis.tar.gz)
  d="$WORK/$s"
  if [[ ! -d "$d" ]]; then
    echo "  [tar ] $s"
    mkdir -p "$d"
    tar xzf "$tgz" -C "$d" \
        'analysis/feature_linkage/feature_linkage.bedpe' \
        'analysis/clustering/atac/graphclust/differential_accessibility.csv' 
  fi

  bedpe=$(find "$d" -name 'feature_linkage.bedpe' | head -1)
  [[ -n "$bedpe" ]] || { echo "  [FATAL] $s 里找不到 feature_linkage.bedpe"; exit 1; }

  # The complete peak set is the background the odds ratio is computed against.
  # 10x does not ship it directly; as the supplementary methods note, it is
  # recovered from the Feature ID column of the differential accessibility
  # output, which lists every peak regardless of significance.
  peaks="$d/peaks_all.bed"
  if [[ ! -s "$peaks" ]]; then
    # Any clustering under analysis/clustering/<atac|gex>/*/ carries the full
    # peak list. Prefer ATAC graphclust; fall back to whichever exists.
    da="$d/analysis/clustering/atac/graphclust/differential_accessibility.csv"
    [[ -s "$da" ]] || da=$(find "$d" -name 'differential_accessibility.csv' | head -1)
    [[ -n "$da" && -s "$da" ]] || { echo "  [FATAL] $s 里找不到差异可及性输出"; exit 1; }
    echo "    peak 集来源: ${da#$d/}"
    python3 - "$da" "$peaks" <<'PY'
import csv, re, sys
src, dst = sys.argv[1], sys.argv[2]
seen, n = set(), 0
with open(src, newline="") as fh, open(dst, "w") as out:
    for row in csv.DictReader(fh):
        fid = row.get("Feature ID") or row.get("Feature Id") or next(iter(row.values()))
        m = re.match(r"^(chr[^:]+)[:_](\d+)[-_](\d+)$", str(fid).strip())
        if not m or fid in seen:
            continue
        seen.add(fid); n += 1
        out.write(f"{m.group(1)}\t{m.group(2)}\t{m.group(3)}\n")
print(f"    peak 集: {n:,} 个 -> {dst}")
PY
  fi

  echo "  [run ] $s  窗口 ±${WIN}"
  python3 code/analysis/p5_15_promoter_OR_external.py \
      --bedpe "$bedpe" --peaks "$peaks" --tss "$TSS" \
      --label "$s" --out "$OUT" --window "$WIN"
done

echo
echo "汇总："
python3 - "$OUT" <<'PY'
import csv, glob, os, sys
rows = []
for f in sorted(glob.glob(os.path.join(sys.argv[1], "p5_promoterOR_*.csv"))):
    rows.append(list(csv.DictReader(open(f)))[0])
print(f"{'label':<32}{'window':>10}{'n_link':>10}{'prox':>8}{'OR':>8}{'95% CI':>20}")
for r in rows:
    print(f"{r['label']:<32}{int(r['window']):>10,}{int(r['n_link']):>10,}"
          f"{int(r['prox_link']):>8,}{float(r['OR']):>8.3f}"
          f"{'  ['+format(float(r['lo']),'.3f')+', '+format(float(r['hi']),'.3f')+']':>20}")
PY
