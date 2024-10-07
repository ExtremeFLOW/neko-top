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
#SBATCH --qos=gp_debug

# Ask for n cores placed on R host.
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1

# Time specifications (hh:mm)
#SBATCH --time 00-00:10:00 # 10 minutes

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

if [[ -z "$LSB_JOBNAME" && (($# > 0)) ]]; then
    example=$1
elif [ ! -z "$LSB_JOBNAME" ]; then
    example=$LSB_JOBNAME
else
    printf "ERROR: No example supplied" >&2
    exit 1
fi

# Load the required modules
if [ ! -z $(which module) ]; then
    module --silent load mpi/4.1.4-gcc-12.2.0-binutils-2.39 openblas/0.3.23 cuda/12.2
fi

source functions.sh
run $example

# ==============================   End of File   ==============================
