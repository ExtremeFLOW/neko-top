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

#BSUB -R "rusage[mem=20GB]"
#BSUB -M 20GB

# Time specifications (hh:mm)
#BSUB -W 10:00

# -- Notification options

# Set the email to recieve to and when to recieve it
#BSUB -Ne    # Send notification at completion

# -- Mandatory options, change with great care.

# Definitions of output files.
#BSUB -o output.out
#BSUB -e error.err

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

# ============================================================================ #
# Execute the example

printf "=%.0s" {1..80} && printf "\n"
printf "Running example: %s.\n" $example

# Run the example
printf "=%.0s" {1..80} && printf "\n"
printf "Executing Neko.\n\n"

casefile=$(find . -name "*.case")
if [ -z "$casefile" ]; then
    printf "ERROR: No case file found.\n" >&2
    exit 1
elif [ $(echo "$casefile" | wc -l) -gt 1 ]; then
    printf "ERROR: Multiple case files found.\n" >&2
    exit 1
fi

casename=$(basename -- ${casefile%.*})
printf "See $casename.out for the status output.\n"
{ time $(mpirun --pernode ./neko $casefile 1>$casename.out 2>error.err); } 2>&1

if [ -s "error.err" ]; then
    printf "\nERROR: An error occured during execution. See error.err for details.\n" >&2
    exit 1
else
    printf "\nNeko execution concluded.\n"
fi

# ============================================================================ #
# Move the results to the results folder
results=$RPATH/$example
printf "=%.0s" {1..80} && printf "\n"
printf "Moving files to results folder: \n\t$results\n\n"

# Remove the results folder if it exists and create a new one
mkdir -p $results && rm -fr $results/*

# Move all the nek5000 files to the results folder and compress them.
for nek in $(find ./ -name "*.nek5000"); do
    printf "Archiving:  %s\n" $nek

    base=${nek%.*}
    field=$(ls $base.f*)
    mkdir -p $results/$base
    mv -t $results/$base $nek $field
done
printf "\n"

# Move all files which are not the error or executable files to the log folder
find ./ -type f \
    -not -name "error.err" \
    -not -name "neko" \
    -not -name "output.out" \
    -not -name "*.chkp" \
    -exec mv -t $results {} +

if [ -s "error.err" ]; then
    printf "ERROR: An error occured during execution. See error.err for details.\n"
    exit 1
else
    printf "=%.0s" {1..80} && printf "\n"
    printf "Example concluded.\n"
    printf "=%.0s" {1..80} && printf "\n"
fi

# Remove all but the log files
find ./ -type f -not -name "error.err" -not -name "output.out" -delete

# Clear the output file to indicate successful completion
cp -ft $results output.out
rm -f output.out
touch output.out

# ==============================   End of File   ==============================
