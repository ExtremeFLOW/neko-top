#!/bin/bash -l

# In this file make changes to the SBATCH variables to control the LUMI
# hpc settings for the unsteady mixer POD practice run.

# =============================================================================
# Define the SBATCH options here.

# -- Technical options

# Queue name
#SBATCH --partition=dev-g

# Ask for one GPU-backed Neko rank and four CPU-only Python ranks.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --gpus-per-node=1
#SBATCH --mem=480GB

# Time specifications (dd-hh:mm:ss)
#SBATCH --time 00:10:00

# -- Notification options

# Set the email to receive to and when to receive it
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nobis@kth.se

# -- Mandatory options, change with great care.

# Definitions of output files.
#SBATCH --output output.log
#SBATCH --error error.log
#SBATCH --account=project_465002526

# ============================================================================ #
# Determine if the script is run on the HPC or locally

set -e

if [[ -z "$SLURM_JOB_NAME" && (($# > 0)) ]]; then
    example=$1
elif [ ! -z "$SLURM_JOB_NAME" ]; then
    example=$SLURM_JOB_NAME
else
    printf "ERROR: No example supplied" >&2
    exit 1
fi

ml craype-accel-amd-gfx90a rocm

source functions.sh

export MPICH_GPU_SUPPORT_ENABLED=1
export ATP_ENABLED=true

if [ -f "${MAIN_DIR}/build/pod_runtime.env" ]; then
    source "${MAIN_DIR}/build/pod_runtime.env"
else
    echo "Error: missing ${MAIN_DIR}/build/pod_runtime.env" >&2
    echo "Run ./setup.sh -e after activating the target Python environment." >&2
    exit 1
fi

# Set these explicitly if you want a different case, Python executable, or
# recovery script than the defaults in examples/unsteady_mixer_practice/run.sh.
export CASE_FILE="${CASE_FILE:-unsteady_mixer_practice.case}"
export PYTHON_BIN="${PYTHON_BIN:-$(command -v python3 || command -v python)}"
export PYTHON_SCRIPT="${PYTHON_SCRIPT:-${MAIN_DIR}/scripts/python/pod_state_recover.py}"
export NEKO_RANKS="${NEKO_RANKS:-1}"
export PY_RANKS="${PY_RANKS:-4}"
export NEKO_STARTUP_DELAY="${NEKO_STARTUP_DELAY:-20}"

run $example
