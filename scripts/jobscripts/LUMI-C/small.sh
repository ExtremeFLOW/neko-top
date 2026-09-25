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
#SBATCH --partition=small

# Ask for n cores placed on R host.
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

# Time specifications (dd-hh:mm:ss)
#SBATCH --time 00-00:10:00

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

if [[ -z "$SLURM_JOB_NAME" && (($# > 0)) ]]; then
    example=$1
elif [ ! -z "$SLURM_JOB_NAME" ]; then
    example=$SLURM_JOB_NAME
else
    printf "ERROR: No example supplied" >&2
    exit 1
fi

# Assign the number of threads to use for OpenMP parallel regions
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# ============================================================================ #
# Select which GPU to map to which core
source functions.sh
run $example

# ==============================   End of File   ==============================
