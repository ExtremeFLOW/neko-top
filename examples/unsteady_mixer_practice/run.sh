#!/bin/bash
set -euo pipefail

if [ -z "${MAIN_DIR:-}" ]; then
    MAIN_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
fi

source "${MAIN_DIR}/scripts/pod_run_helpers.sh"

CASE_DIR=$(pwd)
CASE_FILE=$(pod_resolve_case_file "${CASE_DIR}" "${1:-}")
PY_SCRIPT=${2:-"${MAIN_DIR}/scripts/python/pod_state_recover.py"}

NEKO_RANKS=${NEKO_RANKS:-1}
PY_RANKS=${PY_RANKS:-4}
LOG_FILE=${LOG_FILE:-mpmd.log}

if [ -n "${NEKO_BIN:-}" ]; then
    NEKO_EXE="${NEKO_BIN}"
elif [ -x "${CASE_DIR}/neko" ]; then
    NEKO_EXE="${CASE_DIR}/neko"
elif [ -x "${MAIN_DIR}/examples/unsteady_mixer/neko" ]; then
    NEKO_EXE="${MAIN_DIR}/examples/unsteady_mixer/neko"
else
    echo "Error: no Neko executable found for unsteady_mixer_practice." >&2
    echo "Build the examples or set NEKO_BIN explicitly." >&2
    exit 1
fi

pod_prepare_python_runtime "${MAIN_DIR}"
pod_print_runtime_env
pod_launch_shared_mpi "${CASE_FILE}" "${PY_SCRIPT}" "${NEKO_EXE}" \
    "${PY_RANKS}" "${NEKO_RANKS}" "${LOG_FILE}"
