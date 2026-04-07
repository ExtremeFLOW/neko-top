#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd -- "${SCRIPT_DIR}/../.." && pwd)

cd "${SCRIPT_DIR}"

if [ $# -ge 1 ]; then
    CASE_FILE="$1"
else
    case_files=( *.case )
    if [ ${#case_files[@]} -eq 1 ]; then
        CASE_FILE="${case_files[0]}"
    elif [ ${#case_files[@]} -eq 0 ]; then
        echo "Error: no .case file found in $(pwd)" >&2
        exit 1
    else
        echo "Error: multiple .case files found in $(pwd): ${case_files[*]}" >&2
        exit 1
    fi
fi

PY_SCRIPT=${2:-"${ROOT_DIR}/sources/state_recovery/POD_state_recover"\
"/pod_state_recover.py"}

NEKO_RANKS=${NEKO_RANKS:-10}
PY_RANKS=${PY_RANKS:-4}
LOG_FILE=${LOG_FILE:-mpmd.log}

# Resolve ADIOS2_PATH
if [ -z "${ADIOS2_PATH:-}" ]; then
    if [ -n "${ADIOS2_DIR:-}" ]; then
        ADIOS2_PATH="$ADIOS2_DIR"
    elif [ -d "${ROOT_DIR}/external/adios2" ]; then
        ADIOS2_PATH="${ROOT_DIR}/external/adios2"
    elif [ -n "${MAIN_DIR:-}" ] && [ -d "$MAIN_DIR/external/adios2" ]; then
        ADIOS2_PATH="$MAIN_DIR/external/adios2"
    fi
fi

if [ -n "${ADIOS2_PATH:-}" ]; then
    pyver=$(python3 -c \
        'import sys; print("{}.{}".format(sys.version_info.major, '\
'sys.version_info.minor))')
    export PYTHONPATH="${ADIOS2_PATH}/lib/python${pyver}/site-packages:"\
"${PYTHONPATH:-}"
    export LD_LIBRARY_PATH="${ADIOS2_PATH}/lib:${ADIOS2_PATH}/lib64:"\
"${LD_LIBRARY_PATH:-}"
    echo "Using ADIOS2_PATH=${ADIOS2_PATH}"
else
    echo "Warning: ADIOS2_PATH is not set; in-situ Python may fail to " \
         "import ADIOS2." >&2
fi

echo "Using Python:      $(which python3)"
echo "Using mpirun:      $(which mpirun)"
echo "CONDA_PREFIX:      ${CONDA_PREFIX}"
echo "PYTHONPATH (start):${PYTHONPATH:-<unset>}"
echo "LD_LIBRARY_PATH:   ${LD_LIBRARY_PATH:-<unset>}"
echo "---------------------------------------------------------"

mpirun_cmd=(
    mpirun
    --tag-output
    -n "${PY_RANKS}"
    python3 "${PY_SCRIPT}" "${CASE_FILE}"
    :
    -n "${NEKO_RANKS}"
    ./neko "${CASE_FILE}"
)

echo "Launching shared MPI job:"
printf '  %q' "${mpirun_cmd[@]}"
printf '\n'
echo "Combined output:   ${LOG_FILE}"

"${mpirun_cmd[@]}" > "${LOG_FILE}" 2>&1
