#!/bin/bash -l

# In this file make changes to the SBATCH variables to control the LUMI
# hpc settings for the med_mixer POD run.

# =============================================================================
# Define the SBATCH options here.

# -- Technical options

# Queue name
#SBATCH --partition=small-g

# Ask for a full-node layout on LUMI-G: eight GPU-backed Neko ranks and
# forty-eight CPU-only Python ranks, matching the 56 usable CPU cores/node.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --gpus-per-node=8
#SBATCH --mem=480GB

# Time specifications (dd-hh:mm:ss)
#SBATCH --time 23:30:00

# -- Notification options

# Set the email to receive to and when to receive it
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nobis@kth.se

# -- Mandatory options, change with great care.

# Definitions of output files.
#SBATCH --open-mode=append
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

cat <<'EOF' > select_gpu
#!/bin/bash

export ROCR_VISIBLE_DEVICES=${SLURM_LOCALID:-0}
exec "$@"
EOF

chmod +x ./select_gpu
trap 'rm -f ./select_gpu' EXIT

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
# recovery script than the defaults in examples/med_mixer/run.sh.
export CASE_FILE="${CASE_FILE:-unsteady_mixer_practice.case}"
export PYTHON_BIN="${PYTHON_BIN:-$(command -v python3 || command -v python)}"
export PYTHON_SCRIPT="${PYTHON_SCRIPT:-${MAIN_DIR}/scripts/python/pod_state_recover.py}"
export NEKO_RANKS="${NEKO_RANKS:-8}"
export PY_RANKS="${PY_RANKS:-48}"
export NEKO_STARTUP_DELAY="${NEKO_STARTUP_DELAY:-20}"

run $example
