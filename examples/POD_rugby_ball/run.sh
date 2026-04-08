#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd -- "${SCRIPT_DIR}/../.." && pwd)

cd "${SCRIPT_DIR}"

source "${ROOT_DIR}/scripts/pod_run_helpers.sh"

CASE_FILE=$(pod_resolve_case_file "${SCRIPT_DIR}" "${1:-}")
PY_SCRIPT=${2:-"${ROOT_DIR}/scripts/python/pod_state_recover.py"}

NEKO_RANKS=${NEKO_RANKS:-10}
PY_RANKS=${PY_RANKS:-4}
LOG_FILE=${LOG_FILE:-mpmd.log}

pod_ensure_adios2_python "${ROOT_DIR}"

pod_print_runtime_env
pod_launch_shared_mpi "${CASE_FILE}" "${PY_SCRIPT}" "./neko" \
    "${PY_RANKS}" "${NEKO_RANKS}" "${LOG_FILE}"
