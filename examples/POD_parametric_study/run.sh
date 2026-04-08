#!/bin/bash
set -euo pipefail

if [ -z "${MAIN_DIR:-}" ]; then
    MAIN_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
fi

source "${MAIN_DIR}/scripts/pod_run_helpers.sh"

CASE_DIR=$(pwd)

resolve_neko_exe() {
    local example_dir
    local example_parent

    if [ -n "${NEKO_BIN:-}" ]; then
        printf '%s\n' "${NEKO_BIN}"
        return 0
    fi

    if [ -x "${CASE_DIR}/neko" ]; then
        printf '%s\n' "${CASE_DIR}/neko"
        return 0
    fi

    example_dir="${CASE_DIR}"
    if [[ "${CASE_DIR}" == "${MAIN_DIR}/logs/"* ]]; then
        example_dir="${CASE_DIR/${MAIN_DIR}\/logs/${MAIN_DIR}\/examples}"
    fi
    example_parent="$(dirname "${example_dir}")"

    if [ -x "${example_parent}/neko" ]; then
        printf '%s\n' "${example_parent}/neko"
        return 0
    fi

    echo "Error: no Neko executable found." >&2
    return 1
}

prepare_case_dir() {
    if [ -x "${CASE_DIR}/prepare.sh" ]; then
        ./prepare.sh
    fi
}

run_case() {
    local case_path=$1
    local py_script=$2
    local neko_exe=$3
    local case_name
    local base
    local log_file

    case_name="$(basename "${case_path}")"
    base="${case_name%.case}"

    if pod_case_is_pod "${case_path}"; then
        pod_ensure_adios2_python "${MAIN_DIR}"
        log_file=${LOG_FILE:-"mpmd_${base}.log"}
    else
        log_file=${LOG_FILE:-"neko_${base}.log"}
    fi

    pod_print_runtime_env

    if pod_case_is_pod "${case_path}"; then
        pod_launch_shared_mpi "${case_path}" "${py_script}" "${neko_exe}" \
            "${PY_RANKS}" "${NEKO_RANKS}" "${log_file}"
    else
        pod_launch_neko_only "${case_path}" "${neko_exe}" \
            "${NEKO_RANKS}" "${log_file}"
    fi
}

PY_SCRIPT=${2:-"${MAIN_DIR}/scripts/python/pod_state_recover.py"}
NEKO_RANKS=${NEKO_RANKS:-5}
PY_RANKS=${PY_RANKS:-3}

IS_PARAMETRIC_ROOT=false
if [ -x "${CASE_DIR}/generate_cases.sh" ] && \
   [ -f "${CASE_DIR}/checkpoint.case" ] && \
   [ -f "${CASE_DIR}/pod.case" ]; then
    IS_PARAMETRIC_ROOT=true
fi

NEKO_EXE=$(resolve_neko_exe)

if [ -x "${CASE_DIR}/neko" ] && [ -z "${NEKO_BIN:-}" ] && \
   [ "${CASE_DIR}/neko" -ot "${MAIN_DIR}/sources/problem/problem.f90" ]; then
    echo "Warning: ./neko is older than sources/problem/problem.f90." >&2
fi

if { [ "${IS_PARAMETRIC_ROOT}" = false ] || [ $# -ge 1 ]; } && \
   [ -x "${CASE_DIR}/prepare.sh" ]; then
    prepare_case_dir
fi

if [ "${IS_PARAMETRIC_ROOT}" = true ] && [ $# -lt 1 ]; then
    study_dirs=()

    ./generate_cases.sh

    if [ -d "${CASE_DIR}/checkpoint" ]; then
        study_dirs+=("${CASE_DIR}/checkpoint")
    fi

    while IFS= read -r study_dir; do
        study_dirs+=("${study_dir}")
    done < <(
        find "${CASE_DIR}" -mindepth 1 -maxdepth 1 -type d \
            -name 'POD_nm*_is*' | sort
    )

    if [ ${#study_dirs[@]} -eq 0 ]; then
        echo "Error: no generated study cases found in ${CASE_DIR}" >&2
        exit 1
    fi

    for study_dir in "${study_dirs[@]}"; do
        echo "========================================================="
        echo "Running study case: $(basename "${study_dir}")"
        echo "========================================================="
        (
            cd "${study_dir}"
            NEKO_BIN="${NEKO_EXE}" ./run.sh case.case "${PY_SCRIPT}"
        )
    done
    exit 0
fi

if [ $# -ge 1 ]; then
    CASE_PATH="$1"
else
    shopt -s nullglob
    case_files=( *.case )
    shopt -u nullglob

    if [ ${#case_files[@]} -eq 0 ]; then
        echo "Error: no .case files found in ${CASE_DIR}" >&2
        exit 1
    fi
    if [ ${#case_files[@]} -gt 1 ]; then
        echo "Error: multiple .case files found in ${CASE_DIR}" >&2
        exit 1
    fi
    CASE_PATH="${case_files[0]}"
fi

run_case "${CASE_PATH}" "${PY_SCRIPT}" "${NEKO_EXE}"
