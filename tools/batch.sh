#!/usr/bin/env bash

# -----------------------------
# User configuration
# -----------------------------
CapDAT="./build/capsid_analyzer"
PDBIDS=(1cwp) # eg (1cwp 4btg 1stm)
CANONICAL_FOLDS=(2_0) # eg (2_0 2_1 3_0 3_1 5_0)
CYL_RADIUS=() # eg (35 70 105), if empty it will run just for the default value

FORCE_RERUN="false"       # true = rerun even if sentinel/log already exists
CLEAN_RESULTS_DIR="false" # passed through to CapDAT

set -euo pipefail

# -e: stop on command errors
# -u: stop on unset variables
# -o pipefail: stop hiding failures inside pipelines

# ==========================================================
# CapDAT batch geometry-analysis runner
#
# Features:
# - processes multiple PDBIDs
# - ensures data/ exists
# - downloads missing VIPERdb .vdb.gz inputs with wget
# - decompresses with gunzip
# - runs all requested folds independently
# - writes per-run logs under results/[PDBID]/[fold_name]/output.log
# - keeps a batch summary CSV
# - supports skip-if-results-exist and force-rerun behavior
# - uses temporary download files for safer downloads
# ==========================================================

DATA_DIR="data"
RESULTS_DIR="results"
VIPERDB_BASE_URL="https://viperdb.org/resources/OLIGOMERS"

# Behavior flags
DOWNLOAD_TRIES=3
DOWNLOAD_TIMEOUT=30

# Batch summary outputs
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
SUMMARY_DIR="${RESULTS_DIR}/batch_runs"
SUMMARY_CSV="${SUMMARY_DIR}/batch_summary_${TIMESTAMP}.csv"

# -----------------------------
# Logging helpers
# -----------------------------
log_info() {
  echo "[INFO] $*"
}

log_warn() {
  echo "[WARN] $*" >&2
}

log_error() {
  echo "[ERROR] $*" >&2
}

# -----------------------------
# Dependency checks
# -----------------------------
require_command() {
  local cmd="$1"
  command -v "${cmd}" >/dev/null 2>&1 || {
    log_error "Required command not found: ${cmd}"
    exit 1
  }
}

check_dependencies() {
  require_command wget
  require_command gunzip
  require_command date

  if [ ! -x "${CapDAT}" ]; then
    log_error "CapDAT executable not found or not executable: ${CapDAT}"
    exit 1
  fi
}

# -----------------------------
# Directory helpers
# -----------------------------
ensure_dir() {
  local dir="$1"
  if [ ! -d "${dir}" ]; then
    log_info "Creating directory: ${dir}"
    mkdir -p "${dir}"
  fi
}

ensure_data_dir() {
  ensure_dir "${DATA_DIR}"
}

ensure_summary_dir() {
  ensure_dir "${SUMMARY_DIR}"
}

# -----------------------------
# Input file preparation
# -----------------------------
normalize_pdbid() {
  local raw="$1"
  local norm
  norm="$(printf '%s' "${raw}" | tr '[:upper:]' '[:lower:]')"

  if [ -z "${norm}" ]; then
    log_error "Encountered empty PDBID after normalization"
    return 1
  fi

  printf '%s\n' "${norm}"
}

download_vdb_gz() {
  local pdbid="$1"
  local remote_url="$2"
  local final_gz="$3"
  local tmp_gz="${final_gz}.tmp"

  rm -f "${tmp_gz}"

  log_info "Downloading missing input for PDBID=${pdbid}"
  log_info "Remote URL: ${remote_url}"

  if wget --no-check-certificate --tries="${DOWNLOAD_TRIES}" --timeout="${DOWNLOAD_TIMEOUT}" -O "${tmp_gz}" "${remote_url}"; then
    mv -f "${tmp_gz}" "${final_gz}"
    log_info "Download completed: ${final_gz}"
    return 0
  else
    log_error "Download failed for PDBID=${pdbid}"
    rm -f "${tmp_gz}" "${final_gz}"
    return 1
  fi
}

ensure_vdb_file() {
  local pdbid="$1"
  local vdb_file="${DATA_DIR}/${pdbid}_full.vdb"
  local gz_file="${vdb_file}.gz"
  local remote_url="${VIPERDB_BASE_URL}/${pdbid}_full.vdb.gz"

  if [ -f "${vdb_file}" ]; then
    if [ ! -s "${vdb_file}" ]; then
      log_warn "Existing input file is empty; removing and re-downloading: ${vdb_file}"
      rm -f "${vdb_file}"
    else
      log_info "Found local input file: ${vdb_file}"
      printf '%s\n' "${vdb_file}"
      return 0
    fi
  fi

  if ! download_vdb_gz "${pdbid}" "${remote_url}" "${gz_file}"; then
    return 1
  fi

  log_info "Uncompressing archive: ${gz_file}"
  if gunzip -f "${gz_file}"; then
    log_info "Decompressed successfully: ${vdb_file}"
  else
    log_error "gunzip failed for archive: ${gz_file}"
    rm -f "${gz_file}" "${vdb_file}"
    return 1
  fi

  if [ ! -f "${vdb_file}" ]; then
    log_error "Expected decompressed file not found: ${vdb_file}"
    return 1
  fi

  if [ ! -s "${vdb_file}" ]; then
    log_error "Decompressed file is empty: ${vdb_file}"
    rm -f "${vdb_file}"
    return 1
  fi

  printf '%s\n' "${vdb_file}"
}

# -----------------------------
# Result bookkeeping
# -----------------------------
init_summary_csv() {
  cat > "${SUMMARY_CSV}" <<EOF
timestamp,pdbid,fold_name,status,exit_code,log_file,note
EOF
}

append_summary_row() {
  local timestamp="$1"
  local pdbid="$2"
  local fold_name="$3"
  local status="$4"
  local exit_code="$5"
  local log_file="$6"
  local note="$7"

  printf '%s,%s,%s,%s,%s,%s,%s\n' \
    "${timestamp}" \
    "${pdbid}" \
    "${fold_name}" \
    "${status}" \
    "${exit_code}" \
    "${log_file}" \
    "${note}" >> "${SUMMARY_CSV}"
}

# -----------------------------
# Run policy helpers
# -----------------------------
should_skip_run() {
  local log_file="$1"

  if [ "${FORCE_RERUN}" = "true" ]; then
    return 1
  fi

  if [ -f "${log_file}" ] && [ -s "${log_file}" ]; then
    return 0
  fi

  return 1
}

run_one_fold() {
  local pdbid="$1"
  local fold_name="$2"
  local cyl_radius="${3:-}"
  local run_timestamp
  run_timestamp="$(date +%Y-%m-%dT%H:%M:%S)"

  local run_label="${fold_name}"
  local log_dir="${RESULTS_DIR}/${pdbid}/${fold_name}"
  local capdat_cmd=()

  if [ -n "${cyl_radius}" ]; then
    run_label="${fold_name}_cyl_radius_${cyl_radius}"
    log_dir="${RESULTS_DIR}/${pdbid}/${fold_name}/cyl_radius_${cyl_radius}"
  fi

  local log_file="${log_dir}/output.log"

  ensure_dir "${log_dir}"

  if should_skip_run "${log_file}"; then
    log_warn "Skipping existing run for PDBID=${pdbid}, run=${run_label}"
    log_warn "Existing log file found: ${log_file}"
    append_summary_row "${run_timestamp}" "${pdbid}" "${run_label}" "skipped" "0" "${log_file}" "existing_log_detected"
    return 0
  fi

  log_info "=================================================="
  log_info "Running geometry-analysis"
  log_info "PDBID      : ${pdbid}"
  log_info "fold_name  : ${fold_name}"
  if [ -n "${cyl_radius}" ]; then
    log_info "cyl radius : ${cyl_radius}"
  else
    log_info "cyl radius : (default)"
  fi
  log_info "log file   : ${log_file}"
  log_info "force rerun: ${FORCE_RERUN}"
  log_info "clean dir  : ${CLEAN_RESULTS_DIR}"
  log_info "=================================================="

  capdat_cmd=(
    "${CapDAT}" -i "${pdbid}"
    -l "${log_file}"
    --geometry-analysis
    --geometry_fold_name "${fold_name}"
    --clean-results-dir "${CLEAN_RESULTS_DIR}"
  )

  if [ -n "${cyl_radius}" ]; then
    capdat_cmd+=(--geometry_cylinder_radius "${cyl_radius}")
  fi

  set +e
  "${capdat_cmd[@]}"
  local status=$?
  set -e

  if [ ${status} -ne 0 ]; then
    log_error "Run failed for PDBID=${pdbid}, run=${run_label}, exit code=${status}"
    append_summary_row "${run_timestamp}" "${pdbid}" "${run_label}" "fail" "${status}" "${log_file}" "capdat_run_failed"
    return ${status}
  fi

  log_info "Run completed for PDBID=${pdbid}, run=${run_label}"
  append_summary_row "${run_timestamp}" "${pdbid}" "${run_label}" "ok" "0" "${log_file}" "completed"
  return 0
}

print_final_summary() {
  local total_runs="$1"
  local ok_runs="$2"
  local failed_runs="$3"
  local skipped_runs="$4"
  local download_failed_pdbids="$5"

  echo
  echo "=================================================="
  echo "Batch process finished"
  echo "=================================================="
  echo "Total requested runs : ${total_runs}"
  echo "Successful runs      : ${ok_runs}"
  echo "Failed runs          : ${failed_runs}"
  echo "Skipped runs         : ${skipped_runs}"
  echo "PDBIDs skipped early : ${download_failed_pdbids}"
  echo "Summary CSV          : ${SUMMARY_CSV}"
  echo "=================================================="
}

# -----------------------------
# Main
# -----------------------------
main() {
  check_dependencies
  ensure_data_dir
  ensure_summary_dir
  init_summary_csv

  log_info "=================================================="
  log_info "CapDAT batch geometry-analysis"
  log_info "Executable        : ${CapDAT}"
  log_info "Data dir          : ${DATA_DIR}"
  log_info "Results dir       : ${RESULTS_DIR}"
  log_info "PDBIDs            : ${PDBIDS[*]}"
  log_info "Canonical folds   : ${CANONICAL_FOLDS[*]}"
  if [ "${#CYL_RADIUS[@]}" -gt 0 ]; then
    log_info "Cylinder radii    : ${CYL_RADIUS[*]}"
  else
    log_info "Cylinder radii    : (default from capsid_analyzer)"
  fi
  log_info "FORCE_RERUN       : ${FORCE_RERUN}"
  log_info "CLEAN_RESULTS_DIR : ${CLEAN_RESULTS_DIR}"
  log_info "Summary CSV       : ${SUMMARY_CSV}"
  log_info "=================================================="

  local total_runs=0
  local ok_runs=0
  local failed_runs=0
  local skipped_runs=0
  local download_failed_pdbids=0

  for raw_pdbid in "${PDBIDS[@]}"; do
    local_pdbid="$(normalize_pdbid "${raw_pdbid}")"

    echo
    echo "##################################################"
    echo "[INFO] Processing PDBID=${local_pdbid}"
    echo "##################################################"

    if ! ensure_vdb_file "${local_pdbid}" >/dev/null; then
      log_error "Could not prepare input file for PDBID=${local_pdbid}"
      log_warn "Skipping all folds for PDBID=${local_pdbid}"
      download_failed_pdbids=$((download_failed_pdbids + 1))
      continue
    fi

    for fold_name in "${CANONICAL_FOLDS[@]}"; do
      if [ "${#CYL_RADIUS[@]}" -gt 0 ]; then
        for cyl_radius in "${CYL_RADIUS[@]}"; do
          total_runs=$((total_runs + 1))

          if run_one_fold "${local_pdbid}" "${fold_name}" "${cyl_radius}"; then
            last_status_row="$(tail -n 1 "${SUMMARY_CSV}")"
            case "${last_status_row}" in
              *",skipped,"*)
                skipped_runs=$((skipped_runs + 1))
                ;;
              *",ok,"*)
                ok_runs=$((ok_runs + 1))
                ;;
              *)
                # defensive fallback
                ok_runs=$((ok_runs + 1))
                ;;
            esac
          else
            failed_runs=$((failed_runs + 1))
          fi
        done
      else
        total_runs=$((total_runs + 1))

        if run_one_fold "${local_pdbid}" "${fold_name}"; then
          last_status_row="$(tail -n 1 "${SUMMARY_CSV}")"
          case "${last_status_row}" in
            *",skipped,"*)
              skipped_runs=$((skipped_runs + 1))
              ;;
            *",ok,"*)
              ok_runs=$((ok_runs + 1))
              ;;
            *)
              # defensive fallback
              ok_runs=$((ok_runs + 1))
              ;;
          esac
        else
          failed_runs=$((failed_runs + 1))
        fi
      fi
    done
  done

  print_final_summary "${total_runs}" "${ok_runs}" "${failed_runs}" "${skipped_runs}" "${download_failed_pdbids}"

  if [ "${failed_runs}" -gt 0 ] || [ "${download_failed_pdbids}" -gt 0 ]; then
    exit 1
  fi

  exit 0
}

main "$@"
