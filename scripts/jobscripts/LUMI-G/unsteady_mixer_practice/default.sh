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

export NEKO_RANKS="${NEKO_RANKS:-1}"
export PY_RANKS="${PY_RANKS:-4}"
export NEKO_RANKS_PER_NODE="${NEKO_RANKS_PER_NODE:-1}"
export PY_RANKS_PER_NODE="${PY_RANKS_PER_NODE:-4}"

export GPU_CORES="${GPU_CORES:-49}"
export CPU_CORES="${CPU_CORES:-2,3,4,5}"

run $example
