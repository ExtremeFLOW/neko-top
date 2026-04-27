#!/bin/bash

# In this file make changes to the SBATCH variables to control the LUMI
# hpc settings for POD state-recovery runs that need GPU Neko ranks and
# CPU-only Python ranks in the same MPI_COMM_WORLD.

# =============================================================================
# Define the SBATCH options here.

# -- Technical options

# Queue name
#SBATCH --partition=standard-g

# Ask for a uniform 8 GPU + 8 CPU-rank layout per node.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --gpus-per-node=8

# Time specifications (dd-hh:mm:ss)
#SBATCH --time 01-00:00:00

# -- Notification options

# Set the email to receive to and when to receive it
#SBATCH --mail-type=END

# -- Mandatory options, change with great care.

# Definitions of output files.
#SBATCH --output output.log
#SBATCH --error error.log

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

parse_first_int() {
    local value=${1:-}

    if [[ "${value}" =~ ^([0-9]+) ]]; then
        printf '%s\n' "${BASH_REMATCH[1]}"
        return 0
    fi

    return 1
}

# ============================================================================ #
# Export the per-node layout consumed by scripts/pod_run_helpers.sh

source functions.sh

export MPICH_GPU_SUPPORT_ENABLED=1
export GPU_CORES="${GPU_CORES:-49,57,17,25,1,9,33,41}"

# Available CPU cores on a LUMI-G node, excluding the GPU-near cores used
# by the 8 Neko ranks above.
export CPU_CORES="${CPU_CORES:-2,3,4,5,6,7,10,11,12,13,14,15,18,19,20,21,22,23,26,27,28,29,30,31,34,35,36,37,38,39,42,43,44,45,46,47,50,51,52,53,54,55,58,59,60,61,62,63,65,66,67,68,69,70,71,73,74,75,76,77,78,79,81,82,83,84,85,86,87,89,90,91,92,93,94,95,97,98,99,100,101,102,103,105,106,107,108,109,110,111,113,114,115,116,117,118,119,121,122,123,124,125,126,127}"

IFS=',' read -r -a gpu_cores_arr <<< "${GPU_CORES}"
IFS=',' read -r -a cpu_cores_arr <<< "${CPU_CORES}"

tasks_per_node=$(parse_first_int "${SLURM_NTASKS_PER_NODE:-${SLURM_TASKS_PER_NODE:-16}}") || {
    echo "ERROR: Could not determine tasks per node from Slurm." >&2
    exit 1
}

job_nodes=$(parse_first_int "${SLURM_JOB_NUM_NODES:-${SLURM_NNODES:-1}}") || {
    echo "ERROR: Could not determine the number of Slurm nodes." >&2
    exit 1
}

export NEKO_RANKS_PER_NODE="${NEKO_RANKS_PER_NODE:-${#gpu_cores_arr[@]}}"
export PY_RANKS_PER_NODE="${PY_RANKS_PER_NODE:-$((tasks_per_node - NEKO_RANKS_PER_NODE))}"

if [ "${PY_RANKS_PER_NODE}" -lt 1 ]; then
    echo "ERROR: Need at least one CPU-side Python rank per node." >&2
    exit 1
fi

if [ "${#cpu_cores_arr[@]}" -lt "${PY_RANKS_PER_NODE}" ]; then
    echo "ERROR: Requested ${PY_RANKS_PER_NODE} Python ranks per node but only ${#cpu_cores_arr[@]} CPU cores were listed." >&2
    exit 1
fi

export NEKO_RANKS="${NEKO_RANKS:-$((job_nodes * NEKO_RANKS_PER_NODE))}"
export PY_RANKS="${PY_RANKS:-$((job_nodes * PY_RANKS_PER_NODE))}"
export POD_NEKO_STARTUP_DELAY="${POD_NEKO_STARTUP_DELAY:-20}"

run $example

# ==============================   End of File   ==============================
