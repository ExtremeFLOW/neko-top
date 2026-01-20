#!/bin/bash
set -euo pipefail

CASE_FILE=${1:-POD_rugby_ball.case}
PY_SCRIPT=${2:-pod_state_recover.py}

NEKO_RANKS=${NEKO_RANKS:-10}
PY_RANKS=${PY_RANKS:-4}
PY_STARTUP_DELAY=${PY_STARTUP_DELAY:-3}

if [ -x ./prepare.sh ]; then
    ./prepare.sh
fi

# conda
if [ -f "/scratch/nobis/miniconda3/etc/profile.d/conda.sh" ]; then
    . "/scratch/nobis/miniconda3/etc/profile.d/conda.sh"
else
    eval "$(/scratch/nobis/miniconda3/bin/conda shell.bash hook)"
fi

conda activate /scratch/nobis/POSTDOC/PYSEMTOOLS/miniconda3/.conda/envs/pySEMTools

# Resolve ADIOS2_PATH
if [ -z "${ADIOS2_PATH:-}" ]; then
    if [ -n "${ADIOS2_DIR:-}" ]; then
        ADIOS2_PATH="$ADIOS2_DIR"
    elif [ -n "${MAIN_DIR:-}" ] && [ -d "$MAIN_DIR/external/adios2" ]; then
        ADIOS2_PATH="$MAIN_DIR/external/adios2"
    fi
fi

if [ -n "${ADIOS2_PATH:-}" ]; then
    pyver=$(python3 -c 'import sys; print("{}.{}".format(sys.version_info.major, sys.version_info.minor))')
    export PYTHONPATH="${ADIOS2_PATH}/lib/python${pyver}/site-packages:${PYTHONPATH:-}"
    export LD_LIBRARY_PATH="${ADIOS2_PATH}/lib:${ADIOS2_PATH}/lib64:${LD_LIBRARY_PATH:-}"
    echo "Using ADIOS2_PATH=${ADIOS2_PATH}"
else
    echo "Warning: ADIOS2_PATH is not set; in-situ Python may fail to import ADIOS2." >&2
fi

echo "Using Python:      $(which python3)"
echo "Using mpirun:      $(which mpirun)"
echo "CONDA_PREFIX:      ${CONDA_PREFIX}"
echo "PYTHONPATH (start):${PYTHONPATH:-<unset>}"
echo "LD_LIBRARY_PATH:   ${LD_LIBRARY_PATH:-<unset>}"
echo "---------------------------------------------------------"



# Start python first
mpirun --tag-output -n "${PY_RANKS}" python3 "${PY_SCRIPT}" "${CASE_FILE}" > python.log 2>&1 &
PY_PID=$!

sleep "${PY_STARTUP_DELAY}"

# Run neko
mpirun --tag-output -n "${NEKO_RANKS}" ./neko "${CASE_FILE}" > neko.log 2>&1 

# If you want to cleanly stop python when neko exits:
# kill "${PY_PID}" 2>/dev/null || true
# wait "${PY_PID}" 2>/dev/null || true
