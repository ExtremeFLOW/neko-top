#!/bin/bash
set -euo pipefail

if [ -z "${MAIN_DIR:-}" ]; then
    MAIN_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
fi

source "${MAIN_DIR}/scripts/mpmd_run_helpers.sh"

CASE_DIR=$(pwd)
CASE_FILE=$(mpmd_resolve_case_file "${CASE_DIR}" "${1:-}")
PY_SCRIPT=${2:-"${MAIN_DIR}/scripts/python/pod_state_recover.py"}

NEKO_RANKS=${NEKO_RANKS:-6}
PY_RANKS=${PY_RANKS:-6}

if [ -x "${CASE_DIR}/prepare.sh" ] && [ ! -f "${CASE_DIR}/box.nmsh" ]; then
    ./prepare.sh
fi

if [ -n "${NEKO_BIN:-}" ]; then
    NEKO_EXE="${NEKO_BIN}"
elif [ -x "./neko" ]; then
    NEKO_EXE="./neko"
elif [ -x "../neko" ]; then
    NEKO_EXE="../neko"
else
    echo "Error: no Neko executable found." >&2
    exit 1
fi

CASE_NAME="$(basename "${CASE_FILE}")"
CASE_BASE="${CASE_NAME%.case}"
LOG_FILE=${LOG_FILE:-"mpmd_${CASE_BASE}.log"}

mpmd_ensure_adios2_python "${MAIN_DIR}"

mpmd_print_runtime_env
mpmd_launch_shared "${CASE_FILE}" "${PY_SCRIPT}" "${NEKO_EXE}" \
    "${PY_RANKS}" "${NEKO_RANKS}" "${LOG_FILE}"
