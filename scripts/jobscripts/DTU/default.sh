#!/bin/bash

# In this file make changes to the BSUB variables to control the DTU hpc
# settings.
#
# In addition, all modules should be loaded and python virtualenv should be
# setup if python is used in either testing or visualisation.
# Modules and python setups can be done in a seperate file and supplied through
# the FILES variable in submit.sh. This will ensure a uniform setup.

# =============================================================================
# Define the BSUB options here.

# --  Technical Options

# Queue name
#BSUB -q "gpua100"

# Ask for n cores placed on R host.
#BSUB -n 4
#BSUB -R "span[ptile=4]"
#BSUB -gpu "num=1:mode=exclusive_process"

# Memory specifications. Amount we need and when to kill the
# program using too much memory.

#BSUB -R "rusage[mem=5GB]"
#BSUB -M 5GB

# Time specifications (hh:mm)
#BSUB -W 02:00

# -- Notification options

# Set the email to recieve to and when to recieve it
#BSUB -Ne    # Send notification at completion

# -- Mandatory options, change with great care.

# Definitions of output files.
#BSUB -oo output.log
#BSUB -eo error.err

# ============================================================================ #
# Determine if the script is run on the HPC or locally

set -e

if [[ ! -f "functions.sh" ]]; then
    exit 0
elif [[ -z "$LSB_JOBNAME" && (($# > 0)) ]]; then
    example=$1
elif [ ! -z "$LSB_JOBNAME" ]; then
    example=$LSB_JOBNAME
else
    printf "ERROR: No example supplied" >&2
    exit 1
fi

# LSF 10 exposes the GPUs allocated to the job through LSB_GPU_ID. Map its
# numeric IDs to UUIDs because CUDA and nvidia-smi may use different ordering.
if [[ -n "${LSB_GPU_ID:-}" && -z "${CUDA_VISIBLE_DEVICES+x}" ]] && \
        command -v nvidia-smi >/dev/null 2>&1; then
    read -r -a gpu_ids <<< "${LSB_GPU_ID//,/ }"
    gpu_uuids=()

    for gpu_id in "${gpu_ids[@]}"; do
        if gpu_uuid=$(nvidia-smi --id="$gpu_id" --query-gpu=uuid \
                --format=csv,noheader,nounits 2>/dev/null); then
            gpu_uuids+=("$gpu_uuid")
        fi
    done

    if ((${#gpu_uuids[@]} == ${#gpu_ids[@]})); then
        export CUDA_VISIBLE_DEVICES="$(IFS=,; printf '%s' "${gpu_uuids[*]}")"
    fi

    unset gpu_id gpu_ids gpu_uuid gpu_uuids
fi

source functions.sh
run $example

# ==============================   End of File   ==============================
