#!/bin/bash

# In this file make changes to the SBATCH variables to control the LUMI
# hpc settings.
#
# In addition, all modules should be loaded and python virtualenv should be
# setup if python is used in either testing or visualisation.
# Modules and python setups can be done in a separate file and supplied through
# the FILES variable in submit.sh. This will ensure a uniform setup.

# =============================================================================
# Define the SBATCH options here.

# --  Technical Options

# Queue name
#SBATCH --partition=standard

# Ask for n cores placed on R host.
#SBATCH --nodes=1
#SBATCH --tasks-per-node=128

# Time specifications (dd-hh:mm:ss)
#SBATCH --time 00-00:05:00

# -- Notification options

# Set the email to receive to and when to receive it
#SBATCH --mail-type=END    # Send notification at completion

# -- Mandatory options, change with great care.

# Definitions of output files.
#SBATCH --output output.log
#SBATCH --error error.log

# ============================================================================ #
# Determine if the script is run on the HPC or locally

set -e

if [[ ! -f "functions.sh" ]]; then
    exit 0
elif [[ -z "$SLURM_JOB_NAME" && (($# > 0)) ]]; then
    example=$1
elif [ ! -z "$SLURM_JOB_NAME" ]; then
    example=$SLURM_JOB_NAME
else
    printf "ERROR: No example supplied" >&2
    exit 1
fi

# ============================================================================ #
# Select which GPU to map to which core
source functions.sh
run $example

# ==============================   End of File   ==============================
