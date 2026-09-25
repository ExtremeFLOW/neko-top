#!/bin/bash
set -euo pipefail

if [ -z "${MAIN_DIR:-}" ]; then
    MAIN_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
fi

CASE_FILE=${CASE_FILE:-low_Re.case}
PYTHON_SCRIPT=${PYTHON_SCRIPT:-"${MAIN_DIR}/scripts/python/pod_state_recover.py"}
NEKO_EXE=${NEKO_BIN:-./neko}
NEKO_RANKS=${NEKO_RANKS:-8}
PY_RANKS=${PY_RANKS:-48}
TOTAL_RANKS=$((NEKO_RANKS + PY_RANKS))
LOG_FILE=${LOG_FILE:-mpmd.log}

if [ -z "${PYTHON_BIN:-}" ]; then
    PYTHON_BIN=$(command -v python3 || command -v python || true)
fi

if [ ! -f "${CASE_FILE}" ]; then
    echo "Error: case file not found: ${CASE_FILE}" >&2
    exit 1
fi

if [ ! -f "${PYTHON_SCRIPT}" ]; then
    echo "Error: Python driver not found: ${PYTHON_SCRIPT}" >&2
    exit 1
fi

if [ -z "${PYTHON_BIN:-}" ]; then
    echo "Error: could not find python3 or python in PATH." >&2
    exit 1
fi

if [ ! -x "${NEKO_EXE}" ]; then
    echo "Error: Neko executable not found or not executable: ${NEKO_EXE}" >&2
    echo "Build examples/low_Re or set NEKO_BIN explicitly." >&2
    exit 1
fi

if [ "${NEKO_RANKS}" -lt 1 ]; then
    echo "Error: NEKO_RANKS must be at least 1." >&2
    exit 1
fi

if [ "${PY_RANKS}" -lt 1 ]; then
    echo "Error: PY_RANKS must be at least 1." >&2
    exit 1
fi

if ! command -v srun >/dev/null 2>&1; then
    echo "Error: srun not found in PATH." >&2
    exit 1
fi

source "${MAIN_DIR}/scripts/mpmd_run_helpers.sh"

mpmd_ensure_adios2_python "${MAIN_DIR}"
mpmd_print_runtime_env

cat <<'EOF' > select_gpu
#!/bin/bash

export ROCR_VISIBLE_DEVICES=${SLURM_LOCALID:-0}
export NEKO_GS_COMM=${NEKO_GS_COMM:-MPI}
sleep "${NEKO_STARTUP_DELAY:-20}"
exec "$@"
EOF

chmod +x ./select_gpu

rm -f ./mpmd.conf

for ((rank=0; rank<NEKO_RANKS; rank++)); do
    echo "${rank} /usr/bin/env NEKO_COMM_ID=0 NEKO_CTRL_PEER_ROOT=${NEKO_RANKS} ./select_gpu ${NEKO_EXE} ${CASE_FILE}" >> mpmd.conf
done

for ((rank=NEKO_RANKS; rank<NEKO_RANKS + PY_RANKS; rank++)); do
    echo "${rank} /usr/bin/env NEKO_COMM_ID=1 NEKO_CTRL_PEER_ROOT=0 ${PYTHON_BIN} ${PYTHON_SCRIPT} ${CASE_FILE}" >> mpmd.conf
done

echo "Launching ${NEKO_RANKS} Neko GPU ranks and ${PY_RANKS} Python ranks"
echo "Case:   ${CASE_FILE}"
echo "Neko:   ${NEKO_EXE}"
echo "Python: ${PYTHON_BIN} ${PYTHON_SCRIPT}"
echo "Output: ${LOG_FILE}"

if ! srun --unbuffered --multi-prog mpmd.conf > "${LOG_FILE}" 2>&1; then
    rm -f ./select_gpu ./mpmd.conf
    echo "Error: shared MPMD launch failed. See ${LOG_FILE}." >&2
    exit 1
fi

rm -f ./select_gpu ./mpmd.conf
