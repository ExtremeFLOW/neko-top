#!/bin/bash

# In this file make changes to the SBATCH variables to control the Marenostrum
# hpc settings.
#
# In addition, all modules should be loaded and python virtualenv should be
# setup if python is used in either testing or visualisation.
# Modules and python setups can be done in a seperate file and supplied through
# the FILES variable in submit.sh. This will ensure a uniform setup.

# =============================================================================
# Define the SBATCH options here.

# --  Technical Options

# Queue name
#SBATCH --qos=acc_debug

# Ask for n cores placed on R host.
#SBATCH --ntasks 2
#SBATCH --cpus-per-task 20
#SBATCH --gres=gpu:2

# Time specifications (hh:mm)
#SBATCH --time 00-00:01:00 # 10 minutes

# -- Notification options

# Set the email to recieve to and when to recieve it
#SBATCH --mail-type=END    # Send notification at completion

# -- Mandatory options, change with great care.

# Definitions of output files.
#SBATCH --output output.log
#SBATCH --error error.err

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

source functions.sh
run $example

# ==============================   End of File   ==============================
