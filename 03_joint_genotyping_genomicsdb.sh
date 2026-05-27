#!/usr/bin/env bash
set -euo pipefail

# Purpose: jointly genotype single-sample gVCFs using GenomicsDBImport and GenotypeGVCFs.
# Main output name is kept as: sole675_growth_all.vcf.gz.

PROJECT_DIR="/path/to/tongue_sole_WGS_project"
REF="/path/to/SoleRef.genome.fa"

GVCF_DIR="${PROJECT_DIR}/gvcf"
GVCF_GLOB="*.g.vcf.gz"
WORK_DIR="${PROJECT_DIR}/genotype_work"
OUT_DIR="${PROJECT_DIR}/vcf"
FINAL_VCF="${OUT_DIR}/sole675_growth_all.vcf.gz"

TOTAL_THREADS=196
PARALLEL_JOBS=8
GC_THREADS=4
IMPORT_READER_THREADS=6
JAVA_MEM_IMPORT="24g"
JAVA_MEM_GENO="16g"

CONTIG_FILE=""
INCLUDE_REGEX=".*"
EXCLUDE_REGEX="^$"
CONTIG_LIMIT=0
DO_CONSOLIDATE=false

GATK="gatk"
SAMTOOLS="samtools"
BCFTOOLS="bcftools"
TABIX="tabix"

DB_ROOT="${WORK_DIR}/genomicsdb_ws"
PER_CONTIG_DIR="${WORK_DIR}/per_contig_vcfs"
LOG_DIR="${WORK_DIR}/logs"
TMP_DIR="${WORK_DIR}/tmp"
CONTIG_LIST="${WORK_DIR}/contigs.list"
SAMPLE_MAP="${WORK_DIR}/sample_map.txt"

mkdir -p "${WORK_DIR}" "${DB_ROOT}" "${PER_CONTIG_DIR}" "${OUT_DIR}" "${LOG_DIR}" "${TMP_DIR}"

need_bin() {
  command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing command in PATH: $1" >&2; exit 1; }
}

wait_for_slots() {
  while [[ "$(jobs -r -p | wc -l)" -ge "${PARALLEL_JOBS}" ]]; do
    wait -n || exit 1
  done
}

trap 'rc=$?; if [[ $rc -ne 0 ]]; then echo "ERROR: script failed. Check logs: ${LOG_DIR}" >&2; fi; jobs -pr | xargs -r kill 2>/dev/null || true; exit $rc' EXIT

need_bin "${GATK}"
need_bin "${SAMTOOLS}"
need_bin "${BCFTOOLS}"
need_bin "${TABIX}"
[[ -f "${REF}" ]] || { echo "ERROR: missing reference genome: ${REF}" >&2; exit 1; }

[[ -f "${REF}.fai" ]] || "${SAMTOOLS}" faidx "${REF}"
DICT="${REF%.*}.dict"
[[ -f "${DICT}" ]] || "${GATK}" CreateSequenceDictionary -R "${REF}" -O "${DICT}"

if [[ -n "${CONTIG_FILE}" && -s "${CONTIG_FILE}" ]]; then
  cp -f "${CONTIG_FILE}" "${CONTIG_LIST}"
else
  awk '{print $1}' "${REF}.fai" \
    | awk -v inc="${INCLUDE_REGEX}" -v exc="${EXCLUDE_REGEX}" '$0 ~ inc && (exc == "" || $0 !~ exc) {print}' \
    > "${CONTIG_LIST}"

  if [[ "${CONTIG_LIMIT}" -gt 0 ]]; then
    head -n "${CONTIG_LIMIT}" "${CONTIG_LIST}" > "${CONTIG_LIST}.tmp"
    mv -f "${CONTIG_LIST}.tmp" "${CONTIG_LIST}"
  fi
fi

[[ -s "${CONTIG_LIST}" ]] || { echo "ERROR: contig list is empty." >&2; exit 1; }

echo "[INFO] Number of contigs: $(wc -l < "${CONTIG_LIST}")"

shopt -s nullglob
GVCFS=( "${GVCF_DIR}"/${GVCF_GLOB} )
(( ${#GVCFS[@]} > 0 )) || { echo "ERROR: no gVCFs found: ${GVCF_DIR}/${GVCF_GLOB}" >&2; exit 1; }

for gvcf in "${GVCFS[@]}"; do
  [[ -f "${gvcf}.tbi" || -f "${gvcf}.csi" ]] || "${TABIX}" -p vcf "${gvcf}"
done

: > "${SAMPLE_MAP}"
for gvcf in "${GVCFS[@]}"; do
  mapfile -t samples < <("${BCFTOOLS}" query -l "${gvcf}")
  if [[ "${#samples[@]}" -ne 1 ]]; then
    echo "ERROR: expected one sample in ${gvcf}, found ${#samples[@]}" >&2
    exit 1
  fi
  printf "%s\t%s\n" "${samples[0]}" "${gvcf}" >> "${SAMPLE_MAP}"
done

if [[ "$(cut -f1 "${SAMPLE_MAP}" | sort | uniq -d | wc -l)" -gt 0 ]]; then
  echo "ERROR: duplicate sample names found in gVCFs." >&2
  cut -f1 "${SAMPLE_MAP}" | sort | uniq -d >&2
  exit 1
fi

echo "[STEP 1] GenomicsDBImport"
while IFS= read -r contig || [[ -n "${contig}" ]]; do
  [[ -z "${contig}" ]] && continue
  wait_for_slots

  workspace="${DB_ROOT}/ws_${contig}"
  done_marker="${workspace}/.IMPORT_DONE"
  log_file="${LOG_DIR}/${contig}_import.log"

  if [[ -f "${done_marker}" ]]; then
    echo "[SKIP] Import already done: ${contig}"
    continue
  fi

  if [[ -d "${workspace}" && ! -f "${done_marker}" ]]; then
    echo "[WARN] Removing incomplete GenomicsDB workspace: ${workspace}"
    rm -rf "${workspace}"
  fi

  (
    echo "[RUN] Import: ${contig}"
    "${GATK}" --java-options "-Xmx${JAVA_MEM_IMPORT} -Djava.io.tmpdir=${TMP_DIR} -XX:ParallelGCThreads=${GC_THREADS} -XX:+ExitOnOutOfMemoryError" \
      GenomicsDBImport \
      --sample-name-map "${SAMPLE_MAP}" \
      --genomicsdb-workspace-path "${workspace}" \
      --genomicsdb-shared-posixfs-optimizations true \
      --consolidate false \
      --batch-size 50 \
      --reader-threads "${IMPORT_READER_THREADS}" \
      -L "${contig}" \
      > "${log_file}" 2>&1
    touch "${done_marker}"
  ) &
done < "${CONTIG_LIST}"
wait

if [[ "${DO_CONSOLIDATE}" == "true" ]]; then
  echo "[STEP 2] Optional GenomicsDB workspace consolidation"
  while IFS= read -r contig || [[ -n "${contig}" ]]; do
    [[ -z "${contig}" ]] && continue
    workspace="${DB_ROOT}/ws_${contig}"
    log_file="${LOG_DIR}/${contig}_consolidate.log"
    [[ -d "${workspace}" ]] || continue
    "${GATK}" --java-options "-Xmx8g -Djava.io.tmpdir=${TMP_DIR}" \
      GenomicsDBImport \
      --genomicsdb-workspace-path "${workspace}" \
      --genomicsdb-update-workspace-path "${workspace}" \
      --consolidate true \
      -L "${contig}" \
      > "${log_file}" 2>&1 || true
  done < "${CONTIG_LIST}"
fi

echo "[STEP 3] GenotypeGVCFs"
while IFS= read -r contig || [[ -n "${contig}" ]]; do
  [[ -z "${contig}" ]] && continue
  wait_for_slots

  workspace="${DB_ROOT}/ws_${contig}"
  import_done="${workspace}/.IMPORT_DONE"
  out_vcf="${PER_CONTIG_DIR}/${contig}.vcf.gz"
  done_marker="${out_vcf}.DONE"
  log_file="${LOG_DIR}/${contig}_genotype.log"

  [[ -f "${import_done}" ]] || { echo "ERROR: missing import marker for ${contig}" >&2; exit 1; }

  if [[ -f "${done_marker}" && -f "${out_vcf}" && ( -f "${out_vcf}.tbi" || -f "${out_vcf}.csi" ) ]]; then
    echo "[SKIP] Genotyping already done: ${contig}"
    continue
  fi

  rm -f "${out_vcf}" "${out_vcf}.tbi" "${out_vcf}.csi" "${done_marker}"

  (
    echo "[RUN] Genotype: ${contig}"
    "${GATK}" --java-options "-Xmx${JAVA_MEM_GENO} -Djava.io.tmpdir=${TMP_DIR} -XX:ParallelGCThreads=${GC_THREADS} -XX:+ExitOnOutOfMemoryError" \
      GenotypeGVCFs \
      -R "${REF}" \
      -V "gendb://${workspace}" \
      -L "${contig}" \
      --only-output-calls-starting-in-intervals true \
      -O "${out_vcf}" \
      > "${log_file}" 2>&1
    "${TABIX}" -f -p vcf "${out_vcf}"
    touch "${done_marker}"
  ) &
done < "${CONTIG_LIST}"
wait

echo "[STEP 4] Concatenate per-contig VCFs"
CHUNK_LIST="${WORK_DIR}/chunks.list"
: > "${CHUNK_LIST}"
while IFS= read -r contig || [[ -n "${contig}" ]]; do
  [[ -z "${contig}" ]] && continue
  echo "${PER_CONTIG_DIR}/${contig}.vcf.gz" >> "${CHUNK_LIST}"
done < "${CONTIG_LIST}"

missing=0
while IFS= read -r vcf || [[ -n "${vcf}" ]]; do
  [[ -f "${vcf}" ]] || { echo "ERROR: missing per-contig VCF: ${vcf}" >&2; missing=1; }
done < "${CHUNK_LIST}"
(( missing == 0 )) || exit 1

if [[ -f "${FINAL_VCF}" && ( -f "${FINAL_VCF}.tbi" || -f "${FINAL_VCF}.csi" ) ]]; then
  echo "[SKIP] Final VCF already exists and is indexed: ${FINAL_VCF}"
else
  tmp_final="${FINAL_VCF}.tmp"
  rm -f "${tmp_final}" "${tmp_final}.tbi" "${tmp_final}.csi"
  "${BCFTOOLS}" concat --threads "${TOTAL_THREADS}" -f "${CHUNK_LIST}" -Oz -o "${tmp_final}" -a
  "${BCFTOOLS}" index -t -f "${tmp_final}"
  mv -f "${tmp_final}" "${FINAL_VCF}"
  mv -f "${tmp_final}.tbi" "${FINAL_VCF}.tbi" 2>/dev/null || true
  mv -f "${tmp_final}.csi" "${FINAL_VCF}.csi" 2>/dev/null || true
fi

echo "[DONE] Joint genotyping finished. Final VCF: ${FINAL_VCF}"
