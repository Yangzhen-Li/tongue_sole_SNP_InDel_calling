#!/usr/bin/env bash
set -euo pipefail

# Purpose: perform final sample/genotype/site-level SNP filtering for high-coverage WGS data.
# Main output prefix is kept as: sole_G675_snp_final.

PROJECT_DIR="/path/to/tongue_sole_WGS_project"
VCF_DIR="${PROJECT_DIR}/vcf"

INPUT_VCF="${VCF_DIR}/sole675_growth_snp_hardfiltered.vcf.gz"
OUT_PREFIX="${VCF_DIR}/sole_G675_snp_final"
LOG_FILE="${VCF_DIR}/vcftools_filter_new.log"

VCFTOOLS="vcftools"
BGZIP="bgzip"
TABIX="tabix"

need_bin() {
  command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing command in PATH: $1" >&2; exit 1; }
}

need_bin "${VCFTOOLS}"
need_bin "${BGZIP}"
need_bin "${TABIX}"
[[ -f "${INPUT_VCF}" ]] || { echo "ERROR: missing input VCF: ${INPUT_VCF}" >&2; exit 1; }

"${VCFTOOLS}" \
  --gzvcf "${INPUT_VCF}" \
  --remove-filtered-all \
  --remove-indels \
  --min-alleles 2 \
  --max-alleles 2 \
  --minDP 5 \
  --maxDP 200 \
  --minGQ 20 \
  --max-missing 0.90 \
  --recode \
  --recode-INFO-all \
  --out "${OUT_PREFIX}" \
  > "${LOG_FILE}" 2>&1

"${BGZIP}" -f "${OUT_PREFIX}.recode.vcf"
"${TABIX}" -f -p vcf "${OUT_PREFIX}.recode.vcf.gz"

echo "[DONE] WGS SNP filtering finished: ${OUT_PREFIX}.recode.vcf.gz"
