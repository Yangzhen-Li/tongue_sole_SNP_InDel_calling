#!/usr/bin/env bash
set -euo pipefail

# Purpose: perform paired-end FASTQ quality control using fastp.
# Input FASTQ names are expected as: SAMPLE_raw_1.fq.gz and SAMPLE_raw_2.fq.gz.
# Output FASTQ names are kept as: SAMPLE_clean_1.fq.gz and SAMPLE_clean_2.fq.gz.

PROJECT_DIR="/path/to/tongue_sole_WGS_project"
RAW_DIR="${PROJECT_DIR}/rawdata"
CLEAN_DIR="${PROJECT_DIR}/cleandata"
ID_LIST="${RAW_DIR}/id_list"

MAX_JOBS=40
THREADS_PER_JOB=4
FASTP="fastp"

mkdir -p "${CLEAN_DIR}"

need_bin() {
  command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing command in PATH: $1" >&2; exit 1; }
}

wait_for_jobs() {
  while [[ "$(jobs -pr | wc -l)" -ge "${MAX_JOBS}" ]]; do
    wait -n || exit 1
  done
}

need_bin "${FASTP}"
[[ -s "${ID_LIST}" ]] || { echo "ERROR: missing or empty ID list: ${ID_LIST}" >&2; exit 1; }

while IFS= read -r sample_id || [[ -n "${sample_id}" ]]; do
  [[ -z "${sample_id}" ]] && continue

  raw_file_1="${RAW_DIR}/${sample_id}_raw_1.fq.gz"
  raw_file_2="${RAW_DIR}/${sample_id}_raw_2.fq.gz"
  clean_file_1="${CLEAN_DIR}/${sample_id}_clean_1.fq.gz"
  clean_file_2="${CLEAN_DIR}/${sample_id}_clean_2.fq.gz"
  json_report="${CLEAN_DIR}/${sample_id}.json"
  html_report="${CLEAN_DIR}/${sample_id}.html"
  log_file="${CLEAN_DIR}/${sample_id}_fastp.log"

  if [[ ! -f "${raw_file_1}" || ! -f "${raw_file_2}" ]]; then
    echo "[WARN] Missing raw FASTQ files for ${sample_id}; skipped." >&2
    continue
  fi

  (
    echo "[RUN] fastp: ${sample_id}"
    "${FASTP}" \
      -i "${raw_file_1}" \
      -I "${raw_file_2}" \
      -o "${clean_file_1}" \
      -O "${clean_file_2}" \
      -j "${json_report}" \
      -h "${html_report}" \
      --thread "${THREADS_PER_JOB}" \
      > "${log_file}" 2>&1
  ) &

  wait_for_jobs
done < "${ID_LIST}"

wait

echo "[DONE] fastp quality control finished. Output directory: ${CLEAN_DIR}"
