#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/.." && pwd)"
matlab_bin="${MATLAB_BIN:-matlab}"
timeout_seconds="${TIMEOUT_SECONDS:-420}"
mark_timeout_seconds="${MARK_TIMEOUT_SECONDS:-120}"

cases=(
  "example1 0.003"
  "example1 0.004"
  "example1 0.005"
  "example1 0.006"
  "example1 0.007"
  "example2 0.003"
  "example2 0.004"
  "example2 0.005"
  "example2 0.006"
  "example2 0.007"
  "example2 0.008"
)
methods=(newton qn1 qn2)

run_method() {
  local example_id="$1"
  local gamma="$2"
  local method="$3"
  local log_file
  local exit_code
  log_file="$(mktemp)"

  echo "RUN ${example_id} ${gamma} ${method}"
  exit_code=0
  if ! timeout "${timeout_seconds}" "${matlab_bin}" -batch \
      "cd('${project_root}'); addpath(fullfile(pwd,'scripts')); run_selected_tolerance_mesh_method('${example_id}', ${gamma}, '${method}');" \
      >"${log_file}" 2>&1; then
    exit_code=$?
  fi

  if grep -Eq '^(saved|skip) ' "${log_file}"; then
    grep -E '^(saved|skip) ' "${log_file}" | tail -n 1
    rm -f "${log_file}"
    return 0
  fi

  if grep -Eq 'Out of memory|smallbin double linked list corrupted' "${log_file}"; then
    echo "EXIT_INSTABILITY ${example_id} ${gamma} ${method}"
  fi

  if [[ "${exit_code}" -eq 124 ]]; then
    echo "MARK_TIMEOUT ${example_id} ${gamma} ${method}"
    timeout "${mark_timeout_seconds}" "${matlab_bin}" -batch \
      "cd('${project_root}'); addpath(fullfile(pwd,'scripts')); run_selected_tolerance_mesh_method('${example_id}', ${gamma}, '${method}', struct('mark_timeout', true));" \
      >/dev/null 2>&1 || true
    rm -f "${log_file}"
    return 0
  fi

  cat "${log_file}"
  rm -f "${log_file}"
  return "${exit_code}"
}

cd "${project_root}"
for case_spec in "${cases[@]}"; do
  read -r example_id gamma <<<"${case_spec}"
  for method in "${methods[@]}"; do
    run_method "${example_id}" "${gamma}" "${method}"
  done
done
