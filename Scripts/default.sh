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
#BSUB -q "gpuv100"

# Ask for n cores placed on R host.
#BSUB -n 4
#BSUB -R "span[ptile=4]"
#BSUB -gpu "num=1:mode=exclusive_process"

# Memory specifications. Amount we need and when to kill the
# program using too much memory.

#BSUB -R "rusage[mem=10GB]"
#BSUB -M 10GB

# Time specifications (hh:mm)
#BSUB -W 00:10

# -- Notification options

# Set the email to recieve to and when to recieve it
#BSUB -Ne    # Send notification at completion

# -- Mandatory options, change with great care.

# Definitions of output files.
#BSUB -o output.out
#BSUB -e error.err

# ============================================================================ #
# Determine if the script is run on the HPC or locally

if [[ -z "$LSB_JOBNAME" && (($# > 0)) ]]; then
    example=$1
elif [ ! -z "$LSB_JOBNAME" ]; then
    example=$LSB_JOBNAME
else
    echo "No example supplied" >&2
    exit 1
fi

# Load the required modules
if [ ! -z $(which module) ]; then
    module --silent load mpi/4.1.4-gcc-12.2.0-binutils-2.39 openblas/0.3.23 cuda/12.2
fi

casefile=$(find . -name "*.case")

START=$(date +%s)
mpirun --pernode ./neko $casefile 1>neko.out
END=$(date +%s)
DIFF=$(($END - $START))

echo "It took $DIFF seconds"

if [ ! -s "error.err" ]; then
    rm -fr $RPATH/$example && mkdir -p $RPATH/$example

    # Move the field and boundary files to the results folder
    mkdir -p $RPATH/$example/field $RPATH/$example/bdry
    mv -t $RPATH/$example/field field*
    mv -t $RPATH/$example/bdry bdry*

    # Move all files which are not the error or executable files to the log folder
    find ./ -type f -not -name "error.err" -not -name "neko" -exec mv -t $RPATH/$example {} +

    # Remove all but the log files
    find ./ -type f -not -name "error.err" -not -name "output.out" -delete

    rm -f output.out
    touch output.out
fi

# ==============================   End of File   ==============================
