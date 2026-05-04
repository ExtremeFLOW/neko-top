#!/bin/bash
set -euo pipefail

if [ -z "${MAIN_DIR:-}" ]; then
    MAIN_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
fi

source "${MAIN_DIR}/scripts/mpmd_run_helpers.sh"

CASE_DIR=$(pwd)
CASE_FILE=$(mpmd_resolve_case_file "${CASE_DIR}" "${1:-}")
PY_SCRIPT=${2:-"${MAIN_DIR}/scripts/python/pod_state_recover.py"}

NEKO_RANKS=${NEKO_RANKS:-10}
PY_RANKS=${PY_RANKS:-4}
LOG_FILE=${LOG_FILE:-pySEMTools.log}

mpmd_ensure_adios2_python "${MAIN_DIR}"

mpmd_print_runtime_env
mpmd_launch_shared "${CASE_FILE}" "${PY_SCRIPT}" "./neko" \
    "${PY_RANKS}" "${NEKO_RANKS}" "${LOG_FILE}"
