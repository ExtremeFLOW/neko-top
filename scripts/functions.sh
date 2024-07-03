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
        printf "ERROR: An error occured during execution. " >&2
        printf "See error.err for details.\n" >&2
        exit 1
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

    # ------------------------------------------------------------------------ #
    # Ensure the environment is set up

    [ -z "$NEKO_DIR" ] && NEKO_DIR="$MAIN_DIR/external/neko"
    if [ -z "$JSON_FORTRAN_DIR" ]; then
        JSON_FORTRAN_DIR="$MAIN_DIR/external/json-fortran"
    fi

    JSON_FORTRAN=$(find $JSON_FORTRAN_DIR -type d \
        -exec test -f '{}'/libjsonfortran.so \; -print)
    export LD_LIBRARY_PATH="$JSON_FORTRAN:$LD_LIBRARY_PATH"

    # Run preparation if it exists
    if [ -f "prepare.sh" ]; then
        printf "=%.0s" {1..80} && printf "\n"
        printf "Preparing example.\n\n"

        ./prepare.sh

        if [ -s "error.err" ]; then
            printf "\nERROR: An error occured during preparation. " >&2
            printf "See error.err for details.\n" >&2
            exit 1
        else
            printf "\nPreparation concluded.\n"
        fi
    fi

    # ------------------------------------------------------------------------ #
    # Find the Neko executable

    if [ -f neko ]; then
        neko=$(realpath ./neko)
    elif [ ! -z "$(ls *.f90 2>/dev/null)" ]; then
        printf "=%.0s" {1..80} && printf "\n"
        printf "Building user Neko\n"
        $NEKO_DIR/bin/makeneko *.f90
        neko=$(realpath ./neko)
    else
        neko=$(realpath $NEKO_DIR/bin/neko)
    fi

    if [ ! -f "$neko" ]; then
        printf "ERROR: Neko executable not found." >&2
        exit 1
    fi
    export neko
}

# ============================================================================ #
# Define the run function

function run {

    function handler() {
        printf "An error occured" >error.err
    }
    trap 'handler' SIGINT

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

    casefile=$(find . -name "*.case")
    if [ -z "$casefile" ]; then
        printf "ERROR: No case file found.\n" >&2
        exit 1
    fi

    if [ -f "run.sh" ]; then
        export PATH="$NEKO_DIR/bin:$PATH"

        logfile=$(basename -- $(dirname $(realpath $0)))

        { time ./run.sh 2>error.err; } 2>&1
    else
        casename=$(basename -- ${casefile%.*})
        printf "See $casename.log for the status output.\n"

        {
            time $(mpirun --pernode $neko $casefile 1>$casename.log 2>error.err)
        } 2>&1
    fi

    if [ -s "error.err" ]; then
        printf "\nERROR: An error occured during execution. " >&2
        printf "See error.err for details.\n" >&2
        exit 1
    else
        printf "\nNeko execution concluded.\n"
    fi

    # ------------------------------------------------------------------------ #
    # Clean up the results
    cleanup

}
