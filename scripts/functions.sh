#!/bin/bash
# ============================================================================ #
# Define helpers for the run script

function cleanup {
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

    # Move all files which are not the error or executable files to the log
    # folder
    find ./ -type f \
        -not -name "error.err" \
        -not -name "neko" \
        -not -name "output.log" \
        -not -name "*.chkp" \
        -exec mv -t $results {} +

    if [ -s "error.err" ]; then
        printf "ERROR: An error occured during execution.\n"
        printf "See error.err for details.\n"
        return 1
    else
        printf "=%.0s" {1..80} && printf "\n"
        printf "Example concluded.\n"
        printf "=%.0s" {1..80} && printf "\n"
    fi

    # Remove all but the log files
    find ./ -type f -not -name "error.err" -not -name "output.log" -delete

    # Clear the output file to indicate successful completion
    cp -ft $results output.log
    rm -f output.log
    touch output.log
}

function prepare {
    set +e

    # ------------------------------------------------------------------------ #
    # Report the Job environment if it exists

    if [ ! -z "$SLURM_JOB_NAME" ]; then
        printf "=%.0s" {1..80} && printf "\n"
        printf "SLURM Job: %s\n" $SLURM_JOB_NAME
        printf "=%.0s" {1..80} && printf "\n"
        printf "Job ID: %s\n" $SLURM_JOB_ID
        printf "Job Node: %s\n" $SLURM_JOB_NODELIST
        printf "Job Partition: %s\n" $SLURM_JOB_PARTITION
        printf "Job Account: %s\n" $SLURM_JOB_ACCOUNT
        printf "Job Time: %s\n" $SLURM_JOB_TIME
        printf "Job Memory: %s\n" $SLURM_JOB_MEMORY
        printf "Job CPUs: %s\n" $SLURM_JOB_CPUS_PER_NODE
        printf "Job GPUs: %s\n" $SLURM_JOB_GPUS
        printf "Job QOS: %s\n" $SLURM_JOB_QOS
        printf "Job Reservation: %s\n" $SLURM_JOB_RESERVATION
        printf "Job Work Directory: %s\n" $SLURM_SUBMIT_DIR
        printf "Job Output Directory: %s\n" $SLURM_SUBMIT_DIR
        printf "Job Output File: %s\n" $SLURM_JOB_NAME
        printf "Job Error File: %s\n" $SLURM_JOB_NAME

    elif [ ! -z "$LSB_JOBNAME" ]; then
        printf "=%.0s" {1..80} && printf "\n"
        printf "LSF10 Job: %s\n" $LSB_JOBNAME
        printf "=%.0s" {1..80} && printf "\n"
        printf "Job ID: %s\n" $LSB_JOBID
        printf "Job CPUs: %s\n" $LSB_DJOB_NUMPROC

    fi

    [ -f $MAIN_DIR/prepare.env ] && source $MAIN_DIR/prepare.env
    if [ ! -z "$(which module 2>>/dev/null)" ]; then
        printf "\nModules:\n"
        module list 2>&1
    fi

    # ------------------------------------------------------------------------ #
    # Ensure the environment is set up

    [ -z "$NEKO_DIR" ] && NEKO_DIR="$MAIN_DIR/external/neko"
    if [ -z "$JSON_FORTRAN_DIR" ]; then
        JSON_FORTRAN_DIR="$MAIN_DIR/external/json-fortran"
    fi

    JSON_FORTRAN=$(find $JSON_FORTRAN_DIR -type d \
        -exec test -f '{}'/libjsonfortran.so \; -print)
    export LD_LIBRARY_PATH="$JSON_FORTRAN:$LD_LIBRARY_PATH"
    export PATH="$NEKO_DIR/bin:$PATH"

    # ------------------------------------------------------------------------ #
    # Run preparation if it exists

    if [ -f "prepare.sh" ]; then
        printf "=%.0s" {1..80} && printf "\n"
        printf "Preparing example.\n\n"
        printf "=%.0s" {1..80} && printf "\n"

        { time ./prepare.sh; } 2>&1

        if [ -s "error.err" ]; then
            printf "\nERROR: An error occured during preparation.\n"
            printf "See error.err for details.\n"
            return 1
        else
            printf "\nPreparation concluded.\n"
        fi
    fi

    # ------------------------------------------------------------------------ #
    # Find the Neko executable

    if [ -f ./neko ]; then
        neko=$(realpath ./neko)
    elif [ ! -z "$(ls *.f90 2>>/dev/null)" ]; then
        printf "=%.0s" {1..80} && printf "\n"
        printf "Building user Neko\n"
        printf "=%.0s" {1..80} && printf "\n"

        $NEKO_DIR/bin/makeneko *.f90
        neko=$(realpath ./neko)
    else
        neko=$(realpath $NEKO_DIR/bin/neko)
    fi

    if [ ! -f "$neko" ]; then
        printf "ERROR: Neko executable not found."
        return 1
    fi
    export neko
}

# ============================================================================ #
# Define the run function

function run {
    set -e

    # ------------------------------------------------------------------------ #
    # Set up the environment and find neko
    prepare

    # ------------------------------------------------------------------------ #
    # Execute the example

    printf "=%.0s" {1..80} && printf "\n"
    printf "Running example: %s.\n" $example

    # Run the example
    printf "=%.0s" {1..80} && printf "\n"
    printf "Executing Neko.\n\n"

    casefile=($(find . -name "*.case"))
    if [[ ${#casefile[@]} -eq 0 ]]; then
        printf "ERROR: No case file found.\n"
        return 1
    elif [[ ${#casefile[@]} -eq 1 ]]; then
        casefile=${casefile[0]}
        logfile=$(basename -- ${casefile%.*}).log
    else
        logfile=$(basename -- $(dirname $(realpath $0))).log
    fi
    printf "See $logfile for the status output.\n"

    if [ -f "run.sh" ]; then
        { time ./run.sh 1>$logfile; } 2>&1
    elif [ ! -z "$SLURM_JOB_NAME" ]; then
        {
            time srun --gpu-bind=single:1 $neko $casefile 1>$logfile
        } 2>&1
    else
        # Look for the number of cores to use
        if [ ! -z "$CUDA_VISIBLE_DEVICES" ]; then
            ncores=$(echo $CUDA_VISIBLE_DEVICES | tr "," "\n" | wc -l)
        elif [ ! -z "$LSB_DJOB_NUMPROC" ]; then
            ncores=$LSB_DJOB_NUMPROC
        fi

        if [ -z "$ncores" ]; then
            ncores="1"
        fi
        { time $(mpirun -n $ncores $neko $casefile 1>$logfile 2>error.err); } 2>&1
    fi

    if [ -s "error.err" ]; then
        printf "\nERROR: An error occured during execution.\n"
        printf "See error.err for details.\n"
        return 1
    else
        printf "\nNeko execution concluded.\n"
    fi

    # ------------------------------------------------------------------------ #
    # Clean up the results
    cleanup
}
