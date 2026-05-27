#!/usr/bin/env bash
set -euo pipefail

# Purpose: perform final SNP filtering for the 20K SNP-array-based target sequencing batch.
# Main output prefix is kept as: sole_834_20k_snp_final.

PROJECT_DIR="/path/to/tongue_sole_20k_target_project"
VCF_DIR="${PROJECT_DIR}/vcf"

INPUT_VCF="${VCF_DIR}/sole_834_20k_snp_hardfiltered.vcf.gz"
OUT_PREFIX="${VCF_DIR}/sole_834_20k_snp_final"
LOG_FILE="${VCF_DIR}/vcftools_filter.log"

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
  --minGQ 20 \
  --minDP 8 \
  --maxDP 500 \
  --min-meanDP 20 \
  --max-meanDP 400 \
  --max-missing 0.95 \
  --recode \
  --recode-INFO-all \
  --out "${OUT_PREFIX}" \
  > "${LOG_FILE}" 2>&1

"${BGZIP}" -f "${OUT_PREFIX}.recode.vcf"
"${TABIX}" -f -p vcf "${OUT_PREFIX}.recode.vcf.gz"

echo "[DONE] 20K target-sequencing SNP filtering finished: ${OUT_PREFIX}.recode.vcf.gz"
