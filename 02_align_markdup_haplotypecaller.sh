#!/usr/bin/env bash
set -euo pipefail

# Purpose: align clean FASTQ files to the reference genome, mark duplicates,
# and generate one compressed single-sample gVCF per sample using GATK HaplotypeCaller.
# Input FASTQ names are expected as: SAMPLE_clean_1.fq.gz and SAMPLE_clean_2.fq.gz.
# Main outputs: SAMPLE.sorted.bam, SAMPLE.marked.bam, SAMPLE.g.vcf.gz.

PROJECT_DIR="/path/to/tongue_sole_WGS_project"
REF="/path/to/SoleRef.genome.fa"
PICARD_JAR="/path/to/picard.jar"

CLEAN_DIR="${PROJECT_DIR}/cleandata"
BAM_DIR="${PROJECT_DIR}/bam"
MARKED_DIR="${BAM_DIR}/marked_duplicates"
GVCF_DIR="${PROJECT_DIR}/gvcf"
TMP_DIR="${PROJECT_DIR}/tmp"
ID_LIST="${PROJECT_DIR}/rawdata/id_list"

BWA_THREADS=6
SORT_THREADS=6
HC_THREADS=4
MAX_ALIGN_JOBS=30
MAX_MARKDUP_JOBS=30
MAX_HC_JOBS=20
JAVA_MEM_MARKDUP="16g"
JAVA_MEM_HC="16g"

BWA="bwa"
SAMTOOLS="samtools"
GATK="gatk"
JAVA="java"

mkdir -p "${BAM_DIR}" "${MARKED_DIR}" "${GVCF_DIR}" "${TMP_DIR}"

need_bin() {
  command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing command in PATH: $1" >&2; exit 1; }
}

wait_for_jobs() {
  local max_jobs="$1"
  while [[ "$(jobs -pr | wc -l)" -ge "${max_jobs}" ]]; do
    wait -n || exit 1
  done
}

need_bin "${BWA}"
need_bin "${SAMTOOLS}"
need_bin "${GATK}"
need_bin "${JAVA}"
[[ -f "${REF}" ]] || { echo "ERROR: missing reference genome: ${REF}" >&2; exit 1; }
[[ -f "${PICARD_JAR}" ]] || { echo "ERROR: missing Picard jar: ${PICARD_JAR}" >&2; exit 1; }
[[ -s "${ID_LIST}" ]] || { echo "ERROR: missing or empty ID list: ${ID_LIST}" >&2; exit 1; }

[[ -f "${REF}.fai" ]] || "${SAMTOOLS}" faidx "${REF}"
DICT="${REF%.*}.dict"
[[ -f "${DICT}" ]] || "${GATK}" CreateSequenceDictionary -R "${REF}" -O "${DICT}"

echo "[STEP 1] BWA-MEM alignment and BAM sorting"
while IFS= read -r sample_id || [[ -n "${sample_id}" ]]; do
  [[ -z "${sample_id}" ]] && continue

  clean_file_1="${CLEAN_DIR}/${sample_id}_clean_1.fq.gz"
  clean_file_2="${CLEAN_DIR}/${sample_id}_clean_2.fq.gz"
  sorted_bam="${BAM_DIR}/${sample_id}.sorted.bam"
  log_file="${BAM_DIR}/${sample_id}.bwa_sort.log"

  if [[ ! -f "${clean_file_1}" || ! -f "${clean_file_2}" ]]; then
    echo "[WARN] Missing clean FASTQ files for ${sample_id}; skipped." >&2
    continue
  fi

  if [[ -f "${sorted_bam}" && -f "${sorted_bam}.bai" ]]; then
    echo "[SKIP] Sorted BAM exists for ${sample_id}"
    continue
  fi

  (
    echo "[RUN] Align: ${sample_id}"
    "${BWA}" mem -M -t "${BWA_THREADS}" \
      -R "@RG\tID:${sample_id}\tSM:${sample_id}\tPL:ILLUMINA" \
      "${REF}" "${clean_file_1}" "${clean_file_2}" \
      2> "${log_file}" | \
    "${SAMTOOLS}" sort -@ "${SORT_THREADS}" -o "${sorted_bam}" - \
      >> "${log_file}" 2>&1
    "${SAMTOOLS}" index -@ "${SORT_THREADS}" "${sorted_bam}"
  ) &

  wait_for_jobs "${MAX_ALIGN_JOBS}"
done < "${ID_LIST}"
wait

echo "[STEP 2] Mark duplicates with Picard"
while IFS= read -r sample_id || [[ -n "${sample_id}" ]]; do
  [[ -z "${sample_id}" ]] && continue

  sorted_bam="${BAM_DIR}/${sample_id}.sorted.bam"
  marked_bam="${MARKED_DIR}/${sample_id}.marked.bam"
  metrics_file="${MARKED_DIR}/${sample_id}.markedmetrics.txt"
  log_file="${MARKED_DIR}/${sample_id}.markdup.log"

  if [[ ! -f "${sorted_bam}" ]]; then
    echo "[WARN] Missing sorted BAM for ${sample_id}; skipped." >&2
    continue
  fi

  if [[ -f "${marked_bam}" && -f "${marked_bam}.bai" ]]; then
    echo "[SKIP] Marked BAM exists for ${sample_id}"
    continue
  fi

  (
    echo "[RUN] MarkDuplicates: ${sample_id}"
    "${JAVA}" -Xmx"${JAVA_MEM_MARKDUP}" -jar "${PICARD_JAR}" MarkDuplicates \
      I="${sorted_bam}" \
      O="${marked_bam}" \
      M="${metrics_file}" \
      CREATE_INDEX=true \
      REMOVE_DUPLICATES=false \
      TMP_DIR="${TMP_DIR}" \
      > "${log_file}" 2>&1
  ) &

  wait_for_jobs "${MAX_MARKDUP_JOBS}"
done < "${ID_LIST}"
wait

echo "[STEP 3] GATK HaplotypeCaller in GVCF mode"
while IFS= read -r sample_id || [[ -n "${sample_id}" ]]; do
  [[ -z "${sample_id}" ]] && continue

  marked_bam="${MARKED_DIR}/${sample_id}.marked.bam"
  gvcf_file="${GVCF_DIR}/${sample_id}.g.vcf.gz"
  log_file="${GVCF_DIR}/${sample_id}.haplotypecaller.log"

  if [[ ! -f "${marked_bam}" ]]; then
    echo "[WARN] Missing marked BAM for ${sample_id}; skipped." >&2
    continue
  fi

  if [[ -f "${gvcf_file}" && -f "${gvcf_file}.tbi" ]]; then
    echo "[SKIP] gVCF exists for ${sample_id}"
    continue
  fi

  (
    echo "[RUN] HaplotypeCaller: ${sample_id}"
    "${GATK}" --java-options "-Xmx${JAVA_MEM_HC} -Djava.io.tmpdir=${TMP_DIR} -XX:+ExitOnOutOfMemoryError" \
      HaplotypeCaller \
      -R "${REF}" \
      -I "${marked_bam}" \
      -ERC GVCF \
      --native-pair-hmm-threads "${HC_THREADS}" \
      -O "${gvcf_file}" \
      > "${log_file}" 2>&1
  ) &

  wait_for_jobs "${MAX_HC_JOBS}"
done < "${ID_LIST}"
wait

echo "[DONE] Alignment, duplicate marking, and gVCF calling finished."
