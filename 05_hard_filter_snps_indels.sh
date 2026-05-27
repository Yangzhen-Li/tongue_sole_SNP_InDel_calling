#!/usr/bin/env bash
set -euo pipefail

# Purpose: apply GATK hard filters to raw SNP and InDel VCFs.
# Output SNP file name is kept as: sole675_growth_snp_hardfiltered.vcf.gz.

PROJECT_DIR="/path/to/tongue_sole_WGS_project"
VCF_DIR="${PROJECT_DIR}/vcf"

RAW_SNP_VCF="${VCF_DIR}/sole675_growth_rawsnp.vcf.gz"
RAW_INDEL_VCF="${VCF_DIR}/sole675_growth_rawindel.vcf.gz"
FILTERED_SNP_VCF="${VCF_DIR}/sole675_growth_snp_hardfiltered.vcf.gz"
FILTERED_INDEL_VCF="${VCF_DIR}/sole675_growth_indel_hardfiltered.vcf.gz"
LOG_DIR="${VCF_DIR}/logs"

GATK="gatk"
TABIX="tabix"

mkdir -p "${LOG_DIR}"

need_bin() {
  command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing command in PATH: $1" >&2; exit 1; }
}

need_bin "${GATK}"
need_bin "${TABIX}"
[[ -f "${RAW_SNP_VCF}" ]] || { echo "ERROR: missing raw SNP VCF: ${RAW_SNP_VCF}" >&2; exit 1; }
[[ -f "${RAW_INDEL_VCF}" ]] || { echo "ERROR: missing raw InDel VCF: ${RAW_INDEL_VCF}" >&2; exit 1; }

"${GATK}" VariantFiltration \
  -V "${RAW_SNP_VCF}" \
  --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
  --filter-name "SNP_HARD_FILTER" \
  -O "${FILTERED_SNP_VCF}" \
  > "${LOG_DIR}/snp_hard_filter.log" 2>&1

"${GATK}" VariantFiltration \
  -V "${RAW_INDEL_VCF}" \
  --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0" \
  --filter-name "INDEL_HARD_FILTER" \
  -O "${FILTERED_INDEL_VCF}" \
  > "${LOG_DIR}/indel_hard_filter.log" 2>&1

[[ -f "${FILTERED_SNP_VCF}.tbi" ]] || "${TABIX}" -f -p vcf "${FILTERED_SNP_VCF}"
[[ -f "${FILTERED_INDEL_VCF}.tbi" ]] || "${TABIX}" -f -p vcf "${FILTERED_INDEL_VCF}"

echo "[DONE] Hard filtering finished:"
echo "  SNPs : ${FILTERED_SNP_VCF}"
echo "  InDels: ${FILTERED_INDEL_VCF}"
