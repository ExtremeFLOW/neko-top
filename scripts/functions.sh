#!/bin/bash

# ============================================================================ #
# Define the run function

function run {
    set +e # Do not exit on error

    # Find the case file and define the log file
    casefile=($(find . -name "*.case" -or -name "*.json"))
    if [[ ${#casefile[@]} -eq 0 ]]; then
        printf >&2 "ERROR: No case file found.\n"
        return 1
    elif [[ ${#casefile[@]} -eq 1 ]]; then
        casefile=${casefile[0]}
        logfile=$(basename -- ${casefile%.*}).log
    else
        logfile=$(basename -- $(dirname $(realpath $0))).log
    fi

    # Neko-TOP logs to its own stream, next to the Neko log file.
    nekotop_logfile=${logfile%.log}_nekotop.log

    # Check for recoverable errors in the error log
    if [ -s error.log ]; then
        grep "ERROR: Optimizer stopped after reaching the maximum runtime" \
            error.log >/dev/null || return 1

        # If the error is recoverable, clear the error log and continue
        echo "" > error.log
    fi

    if [ -f "$logfile" ] || [ -f "$nekotop_logfile" ]; then
        # Move old log files to folder with counter padded to 2 digits

        old_run=run_$(find ./ -maxdepth 1 -type d -name "run_*" | wc -l)
        old_run=$(printf "%s_%02d" "run" $((10#${old_run#run_} + 1)))
        mkdir -p ./$old_run

        find ./ -maxdepth 1 -not -empty -type f -name "*.log" \
            -not -name "output.log" -not -name "error.log" \
            -exec mv -ft ./$old_run {} \;
        cp -ft ./$old_run output.log
        [ -s error.log ] && cp -ft ./$old_run error.log

        # Reset the log files
        printf "Ready" >./output.log
    fi

    # Run the example
    printf "Executing Neko.\n" > ./output.log
    printf "See $logfile for the status output.\n"
    printf "See $nekotop_logfile for the Neko-TOP status output.\n"
    export NEKO_LOG_FILE=$logfile
    export NEKO_TOP_LOG_FILE=$nekotop_logfile

    # ------------------------------------------------------------------------ #
    # Set up the environment and find neko
    prepare 2>error.log || return 1
    if [ -s ./error.log ]; then
        printf "ERROR: An error occured during preparation.\n"
        printf "See error.log for details.\n"
        return 1
    fi
    rm -fr error.log && touch error.log

    # ------------------------------------------------------------------------ #
    # Execute the example

    printf "=%.0s" {1..80} && printf "\n"
    printf "Running example: %s.\n" $example
    printf "=%.0s" {1..80} && printf "\n"

    TIME_START=$(date +%s)
    if [ -f "run.sh" ]; then
        ./run.sh 2>error.log

    elif [[ -n "$SLURM_JOB_NAME" && -n "$CPU_BIND" ]]; then
        srun -u --cpu-bind=${CPU_BIND} $neko $casefile 2>error.log

    elif command -v srun 2>&1 1>/dev/null; then
        srun -u $neko $casefile 2>error.log

    elif [ -n "$(which mpirun 2>/dev/null)" ]; then
        # Look for the number of cores to use
        if [ ! -z "$NPROCS" ]; then
            ncores=$NPROCS
        elif [ ! -z "$CUDA_VISIBLE_DEVICES" ]; then
            ncores=$(echo $CUDA_VISIBLE_DEVICES | tr "," "\n" | wc -l)
        elif [ ! -z "$LSB_DJOB_NUMPROC" ]; then
            ncores=$LSB_DJOB_NUMPROC
        else
            nsockets="$(lscpu | grep "Socket(s)" | awk '{print $2}')"
            ncores="$(lscpu | grep "Core(s) per socket" | awk '{print $4}')"
            ncores=$((nsockets * ncores))
        fi

        if [ -z "$ncores" ]; then
            ncores=1
        fi

        mpirun --tag-output -n $ncores $neko $casefile 2>error.log

        # Remove all lines printed from mpi rank > 0 and remove the mpi tag
        sed -i '/^\[[0-9]*,[1-9]*\]/d' error.log
        sed -i 's/\[1,0\]<stderr>://g' error.log

    else
        $neko $casefile 2>error.log

    fi
    TIME_END=$(date +%s)

    # ------------------------------------------------------------------------ #
    # Check for errors and normal end

    if [ -s ./error.log ]; then
        printf "ERROR: An error occurred during execution.\n"
        printf "See error.log for details.\n"
        return 1
    fi

    printf "\nExample concluded.\n"
    TIME_DIFF=$((TIME_END - TIME_START))
    printf "Execution time: %02d:%02d:%02d\n" \
        $((TIME_DIFF / 3600)) $((TIME_DIFF % 3600 / 60)) $((TIME_DIFF % 60))

    # ------------------------------------------------------------------------ #
    # Clean up the results
    cleanup 2>error.log || return 1

    if [ -s ./error.log ]; then
        printf "ERROR: An error occurred during cleanup.\n"
        printf "See error.log for details.\n"
        return 1
    fi
}

# ============================================================================ #
# Define the prepare function

function prepare {
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

    source $MAIN_DIR/scripts/dependencies.sh
    find_json_fortran
    find_neko

    # ------------------------------------------------------------------------ #
    # Run preparation if it exists

    if [ -f "prepare.sh" ]; then
        printf "=%.0s" {1..80} && printf "\n"
        printf "Running user provided preparation script.\n"
        printf "=%.0s" {1..80} && printf "\n"

        prep_sh=$(realpath ./prepare.sh)
        if [ -f "./select_gpu" ]; then
            prep_sh="./select_gpu $prep_sh"
        fi

        if [[ -n "$SLURM_JOB_NAME" && -n "$CPU_BIND" ]]; then
            srun --nodes=1 --ntasks=1 $prep_sh
            sleep 1 # Make sure SLURM have time to clean up.
        else
            $prep_sh
        fi

    fi

    # ------------------------------------------------------------------------ #
    # Find the Neko executable

    printf "=%.0s" {1..80} && printf "\n"
    printf "Finding Neko executable.\n"
    printf "=%.0s" {1..80} && printf "\n"

    if [ -f ./neko ]; then
        printf "Using user provided Neko executable.\n"
        neko=$(realpath ./neko)

    elif [ ! -z "$(ls *.f90 2>>/dev/null)" ]; then
        printf "Building user Neko based on the following files\n"
        for f in $(ls *.f90); do printf "\t- %s\n" $f; done

        $NEKO_DIR/bin/makeneko *.f90 2>&1 || echo "makeneko failed" >&2
        neko=$(realpath ./neko)

    else
        printf "Using default Neko executable.\n"
        neko=$(realpath $NEKO_DIR/bin/neko)
    fi

    if [ ! -f "$neko" ]; then
        printf >&2 "ERROR: Neko executable not found."
        return 1
    fi

    if [ -f "./select_gpu" ]; then
        neko="./select_gpu $neko"
    fi

    export neko
}

# ============================================================================ #
# Define the cleanup function

function cleanup {

    if [ -s ./error.log ]; then
        printf >&2 "ERROR: Will not clean up due to errors.\n"
        return 1
    fi

    # ------------------------------------------------------------------------ #
    # Remove links to data files and folders

    find ./ -maxdepth 1 -type l -delete

    # ------------------------------------------------------------------------ #
    # Move the results to the results folder
    results=$RPATH/$example
    printf "=%.0s" {1..80} && printf "\n"
    printf "Cleaning up the results.\n"
    printf "=%.0s" {1..80} && printf "\n"
    printf "Files will be moved to results folder: \n\t$results\n\n"

    # Remove the results folder if it exists and create a new one
    [ -d $results ] && rm -fr $results/*
    mkdir -p $results

    # Move all the nek5000 files to the results folder.
    printf "Archiving nek5000 files.\n"
    for nek in $(find ./ -name "*.nek5000"); do
        printf "\t- %s\n" ${nek##*/}

        base=$(basename ${nek%[0-9]*.*})
        directory=$(dirname ${nek%[0-9]*.*})
        pattern=$directory/${base%.*}

        mkdir -p $pattern
        mv -t $pattern $nek ${nek%[0-9]*.*}[0-9]*.f* 2>/dev/null || true
    done
    printf "\n"

    if [ -s ./error.log ]; then
        printf >&2 "ERROR: An error occurred during archiving.\n"
        return 1
    fi

    # Move all files which are not the error or executable files to the log
    # folder
    rsync -a --no-links --remove-source-files \
        --exclude "output.log" \
        --exclude "error.log" \
        --exclude "neko" \
        --exclude "*.smod" \
        ./ $results

    if [ -s ./error.log ]; then
        printf >&2 "ERROR: An error occurred during rsync.\n"
        return 1
    fi

    # Remove all but the log files
    find ./ -type f -not -name "error.log" -not -name "output.log" -delete
    find ./ -type d -empty -delete

    # ------------------------------------------------------------------------ #
    # Last minute checks

    if [ -s ./error.log ]; then
        printf >&2 "ERROR: An error occurred during clean-up.\n"
        return 1
    fi

    # ------------------------------------------------------------------------ #
    # Clear the output file to indicate successful completion
    printf "Clean-up concluded.\n"
    printf "=%.0s" {1..80} && printf "\n"
    cp -ft $results output.log
    printf "" >output.log
}
