#!/usr/bin/env bash
# =============================================================================
# P5 步骤 18 · 下载四个 10x Cell Ranger ARC 2.0.0 演示数据的 analysis 包
# =============================================================================
# 参考 link 集（feature_linkage.bedpe）与完整 peak 集（每簇差异可及性输出）都在
# 这个包里。10x 不单独发布 feature_linkage.bedpe（该路径返回 403）。
#
# 用法：bash code/analysis/p5_18_download_arc_refs.sh [目标目录]
# 默认目标：refs/arc（已在 .gitignore 中）
# 合计约 1.31 GB。
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
