#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
ROOT_DIR="${SCRIPT_DIR}"
while [ ! -d "${ROOT_DIR}/sources" ] && [ "${ROOT_DIR}" != "/" ]; do
    ROOT_DIR="$(dirname "${ROOT_DIR}")"
done
if [ ! -d "${ROOT_DIR}/sources" ]; then
    echo "Error: could not locate repo root from ${SCRIPT_DIR}" >&2
    exit 1
fi

cd "${SCRIPT_DIR}"

resolve_adios2_root() {
    local candidate

    for candidate in \
        "${ADIOS2_PATH:-}" \
        "${ADIOS2_DIR:-}" \
        "${ROOT_DIR}/external/adios2" \
        "${MAIN_DIR:-}/external/adios2"; do
        if [ -n "${candidate}" ] && \
           [ -x "${candidate}/bin/adios2-config" ]; then
            printf '%s\n' "${candidate}"
            return 0
        fi
    done

    if command -v adios2-config >/dev/null 2>&1; then
        dirname "$(dirname "$(realpath "$(command -v adios2-config)")")"
        return 0
    fi

    return 1
}

ensure_adios2_python() {
    local adios2_root
    local pyver

    if python3 -c 'import adios2.bindings' >/dev/null 2>&1; then
        return 0
    fi

    adios2_root=$(resolve_adios2_root) || {
        echo "Error: could not locate ADIOS2 for the POD Python driver." >&2
        return 1
    }

    pyver=$(python3 -c \
        'import sys; print(f"{sys.version_info.major}.'\
'{sys.version_info.minor}")')

    export ADIOS2_PATH="${adios2_root}"
    export PYTHONPATH="${ADIOS2_PATH}/lib/python${pyver}/site-packages"\
"${PYTHONPATH:+:${PYTHONPATH}}"
    export LD_LIBRARY_PATH="${ADIOS2_PATH}/lib:${ADIOS2_PATH}/lib64"\
"${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

    python3 -c 'import adios2.bindings' >/dev/null 2>&1 || {
        echo "Error: adios2.bindings is still not importable." >&2
        return 1
    }

    echo "Using ADIOS2_PATH=${ADIOS2_PATH}"
}

resolve_neko_exe() {
    local example_dir
    local example_parent

    if [ -n "${NEKO_BIN:-}" ]; then
        printf '%s\n' "${NEKO_BIN}"
        return 0
    fi

    if [ -x "${SCRIPT_DIR}/neko" ]; then
        printf '%s\n' "${SCRIPT_DIR}/neko"
        return 0
    fi

    example_dir="${SCRIPT_DIR}"
    if [[ "${SCRIPT_DIR}" == "${ROOT_DIR}/logs/"* ]]; then
        example_dir="${SCRIPT_DIR/${ROOT_DIR}\/logs/${ROOT_DIR}\/examples}"
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
    if [ -x "${SCRIPT_DIR}/prepare.sh" ]; then
        (cd "${SCRIPT_DIR}" && ./prepare.sh)
    fi
}

run_case() {
    local case_path=$1
    local py_script=$2
    local neko_exe=$3
    local case_name
    local base
    local log_file
    local mpirun_cmd

    case_name="$(basename "${case_path}")"
    base="${case_name%.case}"

    if grep -q '"type"[[:space:]]*:[[:space:]]*"pod"' "${case_path}"; then
        ensure_adios2_python
        log_file=${LOG_FILE:-"mpmd_${base}.log"}
        mpirun_cmd=(
            mpirun
            --tag-output
            -n "${PY_RANKS}"
            /usr/bin/env
            NEKO_COMM_ID=1
            NEKO_CTRL_PEER_ROOT="${PY_RANKS}"
            python3 "${py_script}" "${case_path}"
            :
            -n "${NEKO_RANKS}"
            /usr/bin/env
            NEKO_COMM_ID=0
            NEKO_CTRL_PEER_ROOT=0
            "${neko_exe}" "${case_path}"
        )

        echo "Launching shared MPI job:"
    else
        log_file=${LOG_FILE:-"neko_${base}.log"}
        mpirun_cmd=(
            mpirun
            --tag-output
            -n "${NEKO_RANKS}"
            "${neko_exe}" "${case_path}"
        )

        echo "Launching MPI job:"
    fi

    echo "Using Python:      $(which python3)"
    echo "Using mpirun:      $(which mpirun)"
    echo "CONDA_PREFIX:      ${CONDA_PREFIX:-<unset>}"
    echo "PYTHONPATH (start):${PYTHONPATH:-<unset>}"
    echo "LD_LIBRARY_PATH:   ${LD_LIBRARY_PATH:-<unset>}"
    echo "---------------------------------------------------------"

    printf '  %q' "${mpirun_cmd[@]}"
    printf '\n'
    echo "Output:            ${log_file}"
    "${mpirun_cmd[@]}" > "${log_file}" 2>&1
}

PY_SCRIPT=${2:-"${ROOT_DIR}/sources/state_recovery/POD_state_recover"\
"/pod_state_recover.py"}
NEKO_RANKS=${NEKO_RANKS:-5}
PY_RANKS=${PY_RANKS:-3}

IS_PARAMETRIC_ROOT=false
if [ -x "${SCRIPT_DIR}/generate_cases.sh" ] && \
   [ -f "${SCRIPT_DIR}/checkpoint.case" ] && \
   [ -f "${SCRIPT_DIR}/pod.case" ]; then
    IS_PARAMETRIC_ROOT=true
fi

NEKO_EXE=$(resolve_neko_exe)

if [ -x "${SCRIPT_DIR}/neko" ] && [ -z "${NEKO_BIN:-}" ] && \
   [ "${SCRIPT_DIR}/neko" -ot \
   "${ROOT_DIR}/sources/problem/problem.f90" ]; then
    echo "Warning: ./neko is older than sources/problem/problem.f90." >&2
fi

if { [ "${IS_PARAMETRIC_ROOT}" = false ] || [ $# -ge 1 ]; } && \
   [ -x "${SCRIPT_DIR}/prepare.sh" ]; then
    prepare_case_dir
fi

if [ "${IS_PARAMETRIC_ROOT}" = true ] && [ $# -lt 1 ]; then
    study_dirs=()

    ./generate_cases.sh

    if [ -d "${SCRIPT_DIR}/checkpoint" ]; then
        study_dirs+=("${SCRIPT_DIR}/checkpoint")
    fi

    while IFS= read -r study_dir; do
        study_dirs+=("${study_dir}")
    done < <(
        find "${SCRIPT_DIR}" -mindepth 1 -maxdepth 1 -type d \
            -name 'POD_nm*_is*' | sort
    )

    if [ ${#study_dirs[@]} -eq 0 ]; then
        echo "Error: no generated study cases found in ${SCRIPT_DIR}" >&2
        exit 1
    fi

    for study_dir in "${study_dirs[@]}"; do
        echo "========================================================="
        echo "Running study case: $(basename "${study_dir}")"
        echo "========================================================="
        NEKO_BIN="${NEKO_EXE}" \
            "${study_dir}/run.sh" "${study_dir}/case.case" "${PY_SCRIPT}"
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
        echo "Error: no .case files found in ${SCRIPT_DIR}" >&2
        exit 1
    fi
    if [ ${#case_files[@]} -gt 1 ]; then
        echo "Error: multiple .case files found in ${SCRIPT_DIR}" >&2
        exit 1
    fi
    CASE_PATH="${case_files[0]}"
fi

run_case "${CASE_PATH}" "${PY_SCRIPT}" "${NEKO_EXE}"
