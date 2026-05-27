#!/usr/bin/env bash
set -euo pipefail

# Purpose: split the joint-called VCF into raw SNP and raw InDel VCFs.
# Output SNP file name is kept as: sole675_growth_rawsnp.vcf.gz.

PROJECT_DIR="/path/to/tongue_sole_WGS_project"
VCF_DIR="${PROJECT_DIR}/vcf"

INPUT_VCF="${VCF_DIR}/sole675_growth_all.vcf.gz"
RAW_SNP_VCF="${VCF_DIR}/sole675_growth_rawsnp.vcf.gz"
RAW_INDEL_VCF="${VCF_DIR}/sole675_growth_rawindel.vcf.gz"

GATK="gatk"
TABIX="tabix"

need_bin() {
  command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing command in PATH: $1" >&2; exit 1; }
}

need_bin "${GATK}"
need_bin "${TABIX}"
[[ -f "${INPUT_VCF}" ]] || { echo "ERROR: missing input VCF: ${INPUT_VCF}" >&2; exit 1; }

"${GATK}" SelectVariants \
  -V "${INPUT_VCF}" \
  -select-type SNP \
  -O "${RAW_SNP_VCF}"

"${GATK}" SelectVariants \
  -V "${INPUT_VCF}" \
  -select-type INDEL \
  -O "${RAW_INDEL_VCF}"

[[ -f "${RAW_SNP_VCF}.tbi" ]] || "${TABIX}" -f -p vcf "${RAW_SNP_VCF}"
[[ -f "${RAW_INDEL_VCF}.tbi" ]] || "${TABIX}" -f -p vcf "${RAW_INDEL_VCF}"

echo "[DONE] Variant selection finished:"
echo "  SNPs : ${RAW_SNP_VCF}"
echo "  InDels: ${RAW_INDEL_VCF}"
