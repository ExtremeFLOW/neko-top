#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd -- "${SCRIPT_DIR}/../.." && pwd)

cd "${SCRIPT_DIR}"

resolve_case_file() {
    local requested=${1:-}
    local case_files

    if [ -n "${requested}" ]; then
        printf '%s\n' "${requested}"
        return 0
    fi

    shopt -s nullglob
    case_files=( *.case )
    shopt -u nullglob

    if [ ${#case_files[@]} -eq 1 ]; then
        printf '%s\n' "${case_files[0]}"
        return 0
    fi

    if [ ${#case_files[@]} -eq 0 ]; then
        echo "Error: no .case file found in ${SCRIPT_DIR}" >&2
    else
        echo "Error: multiple .case files found in ${SCRIPT_DIR}" >&2
    fi
    exit 1
}

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

CASE_FILE=$(resolve_case_file "${1:-}")
PY_SCRIPT=${2:-"${ROOT_DIR}/sources/state_recovery/POD_state_recover"\
"/pod_state_recover.py"}

NEKO_RANKS=${NEKO_RANKS:-10}
PY_RANKS=${PY_RANKS:-4}
LOG_FILE=${LOG_FILE:-mpmd.log}

ensure_adios2_python

echo "Using Python:      $(which python3)"
echo "Using mpirun:      $(which mpirun)"
echo "CONDA_PREFIX:      ${CONDA_PREFIX:-<unset>}"
echo "PYTHONPATH (start):${PYTHONPATH:-<unset>}"
echo "LD_LIBRARY_PATH:   ${LD_LIBRARY_PATH:-<unset>}"
echo "---------------------------------------------------------"

mpirun_cmd=(
    mpirun
    --tag-output
    -n "${PY_RANKS}"
    /usr/bin/env
    NEKO_COMM_ID=1
    NEKO_CTRL_PEER_ROOT="${PY_RANKS}"
    python3 "${PY_SCRIPT}" "${CASE_FILE}"
    :
    -n "${NEKO_RANKS}"
    /usr/bin/env
    NEKO_COMM_ID=0
    NEKO_CTRL_PEER_ROOT=0
    ./neko "${CASE_FILE}"
)

echo "Launching shared MPI job:"
printf '  %q' "${mpirun_cmd[@]}"
printf '\n'
echo "Combined output:   ${LOG_FILE}"

"${mpirun_cmd[@]}" > "${LOG_FILE}" 2>&1
