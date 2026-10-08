#!/usr/bin/env bash
# =============================================================================
# Step 18: download the analysis bundles for four 10x Cell Ranger ARC 2.0.0
# demonstration datasets.
# =============================================================================
# These bundles are the source of the reference link sets the paper compares
# against. Both pieces live inside them: feature_linkage.bedpe, and the complete
# peak set, which has to be recovered from the per-cluster differential
# accessibility output. 10x does not publish feature_linkage.bedpe on its own --
# that path returns 403 -- so the whole bundle is the only route.
#
# Usage: bash code/analysis/p5_18_download_arc_refs.sh [target dir]
# Default target: refs/arc, which is gitignored. About 1.31 GB in total.
# =============================================================================
set -euo pipefail
OUT="${1:-refs/arc}"
BASE="https://cf.10xgenomics.com/samples/cell-arc/2.0.0"
SAMPLES=(human_brain_3k pbmc_granulocyte_sorted_3k pbmc_granulocyte_sorted_10k lymph_node_lymphoma_14k)
mkdir -p "$OUT"
for s in "${SAMPLES[@]}"; do
  f="$OUT/${s}_analysis.tar.gz"
  if [[ -s "$f" ]]; then echo "  [skip] $f 已存在"; continue; fi
  echo "  [get ] ${s}_analysis.tar.gz"
  curl -fL --retry 3 --retry-delay 5 -C - -o "$f.part" "$BASE/${s}/${s}_analysis.tar.gz"
  mv "$f.part" "$f"
done
echo "完成："; ls -lh "$OUT"
