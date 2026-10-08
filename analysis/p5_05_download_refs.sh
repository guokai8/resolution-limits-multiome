#!/usr/bin/env bash
# =============================================================================
# Step 0.7: download the public reference data. Run this on your own machine.
# =============================================================================
# Usage:
#     bash code/p5_05_download_refs.sh            # only what Layer 3 needs (~50 MB)
#     bash code/p5_05_download_refs.sh --all      # plus GWAS and LDSC (~12 GB)
#     REF=/path/to/refs bash code/p5_05_download_refs.sh
#
# Then verify:
#     python3 code/p5_06_build_tss.py --ref refs --data ~/Desktop/ResearchD --out results
#
# Downloads resume (curl -C -) and existing files are skipped, so the script is
# safe to re-run after an interruption. Each file lands as .part and is renamed
# only on success, so a truncated download is never mistaken for a complete one.
# =============================================================================
set -euo pipefail

REF="${REF:-refs}"
ALL=0
[[ "${1:-}" == "--all" ]] && ALL=1
mkdir -p "$REF"
cd "$REF"

get() {  # get <url> <outfile> -- skip if present, resume if partial
  local url="$1" out="$2"
  if [[ -s "$out" ]]; then echo "  [skip] $out 已存在"; return 0; fi
  echo "  [get ] $out"
  curl -fL --retry 3 --retry-delay 5 -C - -o "$out.part" "$url"
  mv "$out.part" "$out"
}

echo "=============================================================="
echo "1/3 · GENCODE 注释（A2 的唯一阻塞项）"
echo "=============================================================="
# The annotation version must match the one the data were aligned against.
# The Ulm multiome used 10x GRCh38-2020-A, which is GENCODE v32 / Ensembl 98.
# Note the gencode.v30.gene_meta.tsv.gz that ships with some pipelines carries
# only gene name and ID -- no coordinates -- so the full GTF is required to get
# TSS positions at all.
get "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz" \
    "gencode.v32.annotation.gtf.gz"

echo
echo "备用镜像（若上面的 EBI 站点慢，Ctrl-C 后手动用这个）："
echo "  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz"
echo "  ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/"

if [[ $ALL -eq 1 ]]; then
  echo
  echo "=============================================================="
  echo "2/3 · GWAS 汇总统计（命题 C 用，约 2 GB）"
  echo "=============================================================="
  # ALS: van Rheenen 2021, 27,205 cases / 110,881 controls, European ancestry
  get "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027164/GCST90027164_buildGRCh37.tsv.gz" \
      "ALS_vanRheenen2021_GCST90027164_GRCh37.tsv.gz"

  echo
  echo "  以下三个作为**阴性/特异性对照**，请到 GWAS Catalog 页面按需取："
  echo "    AD  (Bellenguez 2022)  GCST90027158"
  echo "    PD  (Nalls 2019)       GCST009325"
  echo "    身高 (阴性对照)          GCST90245848"
  echo "  下载页：https://www.ebi.ac.uk/gwas/studies/<ACCESSION>"

  echo
  echo "=============================================================="
  echo "3/3 · S-LDSC 参考（约 10 GB）"
  echo "=============================================================="
  get "https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_baselineLD_v2.2_ldscores.tgz" \
      "1000G_Phase3_baselineLD_v2.2_ldscores.tgz"
  get "https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_plinkfiles.tgz" \
      "1000G_Phase3_plinkfiles.tgz"
  get "https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/weights_hm3_no_hla.tgz" \
      "weights_hm3_no_hla.tgz"
  get "https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_frq.tgz" \
      "1000G_Phase3_frq.tgz"
  for f in *.tgz; do
    [[ -d "${f%.tgz}" ]] || { echo "  [tar ] $f"; tar xzf "$f"; }
  done
else
  echo
  echo "（跳过 GWAS/LDSC。需要时加 --all）"
fi

echo
echo "=============================================================="
echo "完成。下一步跑校验："
echo "  python3 code/p5_06_build_tss.py --ref $REF --data ~/Desktop/ResearchD --out results"
echo "=============================================================="
ls -lh
