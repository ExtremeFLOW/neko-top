#!/bin/bash

# ============================================================================ #
# Define the help function
function help() {
    echo -e "Execution of Neko cases."

    echo -e "\n\e[4mUsage:\e[0m"
    echo -e "  run.sh [options] [example]"

    echo -e "\n\e[4mDescription:\e[0m"
    echo -e "  This script works as a function which run the examples specified"
    echo -e "  through the command line."
    echo -e ""
    echo -e "  The <EXAMPLE> input refers to the name or pattern of the case"
    echo -e "  file. The examples folder is searched for the specified pattern."
    echo -e "  If multiple matching case files are found, then all of them are"
    echo -e "  run."
    echo -e ""
    echo -e "  See Readme for additional details."
    echo -e ""

    printf "\e[4mOptions:\e[0m\n"
    printf "  -%-1s, --%-10s %-60s\n" "a" "all" "Run all journals available."
    printf "  -%-1s, --%-10s %-60s\n" "c" "clean" "Clean artifacts from previous runs."
    printf "  -%-1s, --%-10s %-60s\n" "d" "delete" "Delete previous runs."
    printf "  -%-1s, --%-10s %-60s\n" "h" "help" "Print help."
    printf "  -%-1s, --%-10s %-60s\n" "n" "neko" "Look for examples in neko."

    printf "\n\e[4mAvailable case files:\e[0m\n"
    for case in $(find $EPATH -name "*.case" 2>/dev/null); do
        printf '  - %-s\n' ${case#$EPATH/}
    done
    exit 0
}
if [ $# -lt 1 ]; then help; fi

# ============================================================================ #
# User defined inputs.
export MAIN_DIR=$(dirname $(realpath $0))
CURRENT_DIR=$(pwd)

# Define all needed folders relative to the project folder. (without trailing /)
export EPATH="$MAIN_DIR/examples"                 # Examples scripts
export RPATH="$MAIN_DIR/results"                  # Result export location
export LPATH="$MAIN_DIR/logs"                     # Logging locations
export SPATH="$MAIN_DIR/scripts/"                 # Scripts folder
export HPATH="$MAIN_DIR/scripts/jobscripts/LSF10" # Submission settings
export DPATH="$MAIN_DIR/data"                     # Meshes

if [ -f "$MAIN_DIR/prepare.sh" ]; then
    source $MAIN_DIR/prepare.sh
fi
[ -z "$NEKO_DIR" ] && export NEKO_DIR="$MAIN_DIR/external/neko"

# End of user inputs
# ============================================================================ #
# Keywords

# Empty input
if [ $# -lt 1 ]; then help; fi

case_files=""
for in in $@; do

    # Opions and flags
    if [[ ${in:0:2} == "--" ]]; then
        case "${in:2}" in
        "all") ALL=true ;;       # Run all examples available
        "clean") CLEAN=true ;;   # Clean logs
        "help") help ;;          # Print help
        "neko") NEKO=true ;;     # Look for example in neko
        "delete") DELETE=true ;; # Delete previous runs
        *) printf '  %-10s %-67s\n' "Invalid option:" "$in" && exit 1 ;;
        esac

    elif [[ ${in:0:1} == "-" ]]; then
        for ((i = 1; i < ${#in}; i++)); do
            case "${in:$i:1}" in
            "a") ALL=true ;;    # Run all examples available
            "c") CLEAN=true ;;  # Clean logs
            "h") help ;;        # Print help
            "n") NEKO=true ;;   # Look for example in neko
            "d") DELETE=true ;; # Delete previous runs
            *) printf '  %-10s %-67s\n' "Invalid option:" "${in:$i:1}" && exit 1 ;;
            esac
        done
    fi
done

if [ $NEKO ]; then
    export EPATH="$NEKO_DIR/examples"
    export RPATH="$RPATH/neko"
    export LPATH="$LPATH/neko"
    export HPATH="$HPATH/neko"
fi

for in in $@; do
    if [ ${in:0:1} == "-" ]; then continue; fi

    # Ignore invalid inputs
    if [[ -z $(ls $EPATH/$in 2>/dev/null) ]]; then
        printf '  %-12s %-67s\n' "Not Found:" "$in"

    # Extract the examples from the input
    elif [[ ! $ALL ]]; then
        for case in $(find $EPATH/$in -name "*.case"); do
            case_files+="${case#$EPATH/} "
        done
    fi
done

if [ $ALL ]; then
    case_files=""
    for case in $(find $EPATH/ -name "*.case"); do
        case_files+="${case#$EPATH/} "
    done
fi

if [ ! "$case_files" ]; then
    exit
fi

# ============================================================================ #
# Handle settings

if [ $DELETE ]; then
    printf 'Do you wish to delete ALL results? [Yes No]\n'
    read -p '> ' yn
    case $yn in
    [Yy]*) echo "Removing..." && rm -fr $RPATH && echo "Results removed" ;;
    *) echo "Results not removed" ;;
    esac
    printf 'Logs have been cleaned.\n'
    rm -fr $LPATH
fi

# ============================================================================ #
# Define functions for running and submitting the examples.

# Function for running the examples
function Run() {
    case=$1
    echo "Launching case on local machine" 1>output.log

    # Run the case and pipe stdout and stderr to the log files
    printf '  %-12s %-s\n' "Started:" "$case"
    ./job_script.sh $case >output.log 2>error.err
}

# Function for submitting the examples
function Submit() {
    case=$1
    echo "Submitting case on HPC" 1>output.log

    export BSUB_QUIET=Y
    bsub -J $case -env "all" <job_script.sh
    printf '  %-12s %-s\n' "Submitted:" "$case"
}
INTERRUPTED=0
function handler() {
    if [ "$MAIN_DIR" != "$(pwd)" ]; then
        printf "Interrupted" >error.err
    fi
    INTERRUPTED=1
}
trap 'handler' SIGINT

# ============================================================================ #
# Run the examples
full_start=$(date +%s.%N)
QUEUE=""

printf "\n\e[4mQueueing case files.\e[0m\n"
for case in $case_files; do
    case_name=$(basename ${case%.case})
    case_dir=$(dirname $case)

    # Define the name of the current exampel, if there are multiple cases in the
    # same folder, we add the case name to the example name.
    example=${case_dir#$EPATH/}
    if [[ ! -f "$EPATH/$case_dir/run.sh" &&
        $(find $EPATH/$case_dir -name "*.case" | wc -l) > 1 ]]; then
        example=$example/$case_name
    fi

    log=$LPATH/$example && mkdir -p $log

    # Setup the log folder
    if [[ -f $log/output.log && "$(head -n 1 $log/output.log)" == "Ready" ]]; then
        continue
    fi

    if [ $CLEAN ]; then
        rm -fr $log/*
    else
        find $log -type f -name "*.log" -delete
    fi

    # Remove old output and error files
    rm -f $log/output.log $log/error.err
    touch $log/output.log $log/error.err

    # Find the setting file for the case recursively
    setting=$HPATH/${case%.*}.sh
    [ ! -f $setting ] && setting=$(dirname $setting)/default.sh
    while [ ! -f $setting ]; do
        setting=$(realpath "$(dirname ${setting%/*})")/default.sh
    done

    # Copy the case files to the log folder
    cp -ft $log $EPATH/$case
    # Copy all data from the case folder to the log folder
    find $EPATH/$case_dir/* -maxdepth 0 \
        -not -name "*.case" -and -not -name "*.nmsh" -exec rsync -r {} $log \;
    # Copy the job script to the log folder
    cp -f $setting $log/job_script.sh
    cp -f $SPATH/functions.sh $log/functions.sh

    # Create symbolic links to the mesh files to avoid copying massive files
    for file in $(find $EPATH/$case_dir -name "*.nmsh" 2>/dev/null); do
        ln -fs $file $log
    done

    if [ -d "$DPATH" ]; then ln -fs $DPATH $log; fi
    if [ -d "$MAIN_DIR/data_local" ]; then ln -fs $MAIN_DIR/data_local $log; fi

    # Indicate that the case is ready to be run
    printf 'Ready' >$log/output.log

    QUEUE="$QUEUE $example"
    if [ -z "$(which bsub)" ]; then
        printf '  %-12s %-s\n' "Queued:" "$example"
    fi
done

# Done with the setup
# ============================================================================ #
# Move to the directory submit or run the code and return

for example in $QUEUE; do

    # Move to the log folder and submit the job
    cd $LPATH/$example
    if [ $INTERRUPTED == 1 ]; then
        printf "Interrupted" >error.err
    elif [ "$(which bsub)" ]; then
        Submit $example
    else
        Run $example
    fi
    cd $CURRENT_DIR
done

if [ ! "$(which bsub)" ]; then
    $MAIN_DIR/status.sh
fi

printf "\n"
# # EOF # #
