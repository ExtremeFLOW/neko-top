#!/bin/bash
set -e # Exit with nonzero exit code if anything fails
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
    printf "  -%-1s, --%-10s %-60s\n" "s" "submit" "Submit the examples to the queue."
    printf "  -%-1s, --%-10s %-60s\n" " " "dry-run" "Dry run the script."

    printf "\n\e[4mAvailable case files:\e[0m\n"
    for case in $(find $EPATH -name "*.case" 2>/dev/null); do
        printf '  - %-s\n' ${case#$EPATH/}
    done
}
if [ $# -lt 1 ]; then help; fi

# ============================================================================ #
# User defined inputs.

# Assign default values to the options
ALL=false
CLEAN=false
NEKO=false
DELETE=false
CLUSTER=""
DRY=false

# List possible options
OPTIONS=all,clean,help,neko,delete,submit:,dry-run
OPT="a,c,h,n,s:,d"

# Parse the inputs for options
PARSED=$(getopt --options=$OPT --longoptions=$OPTIONS --name "$0" -- "$@")
eval set -- "$PARSED"

# Loop through the options and set the variables
while true; do
    case "$1" in
    "-a" | "--all") ALL=true && shift ;;          # Run all examples available
    "-c" | "--clean") CLEAN=true && shift ;;      # Clean logs
    "-h" | "--help") help && exit ;;              # Print help
    "-n" | "--neko") NEKO=true && shift ;;        # Look for example in neko
    "-d" | "--delete") DELETE=true && shift ;;    # Delete previous runs
    "-s" | "--submit") CLUSTER="$2" && shift 2 ;; # Submit to the queue
    "--dry-run") DRY=true && shift ;;             # Dry run

    # End of options
    "--") shift && break ;;
    esac
done

# ============================================================================ #
# Define environment
export MAIN_DIR=$(dirname $(realpath $0))
CURRENT_DIR=$(pwd)

# Execute the preparation script if it exists
if [ -f "$MAIN_DIR/prepare.env" ]; then
    source $MAIN_DIR/prepare.env
fi

# Define all needed folders relative to the project folder. (without trailing /)
export EPATH="$MAIN_DIR/examples"    # Examples scripts
export RPATH="$MAIN_DIR/results"     # Result export location
export LPATH="$MAIN_DIR/logs"        # Logging locations
export SPATH="$MAIN_DIR/scripts/"    # Scripts folder
export DPATH="$MAIN_DIR/data"        # Official data
export DLPATH="$MAIN_DIR/data_local" # Local data

# Define the job script folder
if [ ! -z "$CLUSTER" ]; then
    export HPATH="$MAIN_DIR/scripts/jobscripts/$CLUSTER" # Submission settings
else
    export HPATH="$MAIN_DIR/scripts/jobscripts" # Submission settings
fi

[ -z "$NEKO_DIR" ] && export NEKO_DIR="$MAIN_DIR/external/neko"
export NEKO_DIR=$(realpath $NEKO_DIR)

if [ "$NEKO" == true ]; then
    export EPATH="$NEKO_DIR/examples"
    export RPATH="$RPATH/neko"
    export LPATH="$LPATH/neko"
    export HPATH="$HPATH/neko"
fi

# End of user inputs
# ============================================================================ #
# Find the examples to run

example_list=()
for in in $@; do
    [ "$ALL" == true ] && break

    # Decompose the input into the directory and the base name
    [ $(dirname $in) == "." ] && dir="" || dir=$(dirname $in)
    base=$(basename $in)

    # Extract the examples from the input
    matches=($(find $EPATH/$dir -maxdepth 1 -type d -name "$base"))
    matches+=($(find $EPATH/$dir -maxdepth 1 -type f -name "$base"))
    matches+=($(find $EPATH/$dir -maxdepth 1 -type f -name "$base.case"))

    for match in ${matches[@]}; do
        file_list=()
        if [ -d $match ]; then
            file_list=($(find $match -name "run.sh" 2>/dev/null))
            file_list+=($(find $match -name "*.case" 2>/dev/null))
        fi
        if [ -f $match ]; then
            file_list+=($match)
        fi

        for file in ${file_list[@]}; do
            dir=$(dirname $file)
            if [[ -f $dir/run.sh ]]; then
                example_list+=("${dir#$EPATH/}/run.sh")
            elif [ $(basename $file) == "run.sh" ]; then
                example_list+=("${dir#$EPATH/}")
            else
                example_list+=("${file#$EPATH/}")
            fi
        done

    done
done

if [ "$ALL" == true ]; then
    file_list=($(find $EPATH -name "run.sh" 2>/dev/null))
    file_list+=($(find $EPATH -name "*.case" 2>/dev/null))

    example_list=()
    for file in ${file_list[@]}; do
        dir=$(dirname $file)
        if [[ -f $dir/run.sh ]]; then
            example_list+=("${dir#$EPATH/}/run.sh")
        elif [ $(basename $file) == "run.sh" ]; then
            example_list+=("${dir#$EPATH/}")
        else
            example_list+=("${file#$EPATH/}")
        fi
    done
fi

# Make sure run.sh in parent folders are used if present.
for i in ${!example_list[@]}; do
    example=${example_list[$i]}

    run_file=$(dirname ${example%/run.sh})/run.sh
    while [[ $run_file != "./run.sh" && ! -f $EPATH/$run_file ]]; do
        run_file=$(dirname ${run_file%/run.sh})/run.sh
    done

    if [[ -f $EPATH/$run_file && ${example: -3} == '.sh' ]]; then

        printf >&2 "\e[1;31mInvalid run file:\e[m\n"
        printf >&2 "$EPATH/$example\n"
        printf >&2 "\tNested run files are not allowed.\n"

        unset example_list[$i]
    elif [[ -f $EPATH/$run_file && ${example: -3} != '.sh' ]]; then
        example_list[$i]=$run_file
    fi
done

# Case files may not be nested in example folders
for i in ${!example_list[@]}; do
    example=${example_list[$i]}
    parent=$(dirname ${example%/*.*})
    while [[ $parent != "." &&
        -z "$(find $EPATH/$parent -maxdepth 1 -name '*.case')" ]]; do
        parent=$(dirname $parent)
    done

    if [ $parent != "." ]; then

        printf >&2 "\e[1;31mInvalid example file:\e[m\n"
        printf >&2 "$EPATH/$example\n"
        printf <&2 "\tNested examples are not allowed.\n"
        printf <&2 "\tMove the $example file to the root of example suite\n"
        if [ ${example: -5} == ".case" ]; then
            printf >&2 "\tor create a run.sh file in the parent folder.\n"
        fi

        unset example_list[$i]
    fi
done

# Remove duplicates and check for nested examples
example_list=($(echo "${example_list[@]}" | tr ' ' '\n' | sort -u))

# Check if any examples were found, if not, exit.
if [ -z $example_list ]; then
    printf "No examples found.\n" >&2 && exit
fi

# ============================================================================ #
# Handle settings

if [ "$DELETE" == true ]; then
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
    cd $LPATH/$example
    printf '\t%-12s %-s\n' "Started:" "$1"
    source $SPATH/functions.sh
    run $1 1>output.log 2>error.err
    cd $CURRENT_DIR
}

# Function for submitting the examples
function Submit() {

    # Find the setting file for the case recursively
    setting=$HPATH/${case%.*}.sh
    while [[ ! -f $setting && ! -z "$setting" ]]; do
        setting=$(dirname ${setting%/default.sh})/default.sh
    done
    setting=$(realpath $setting)

    if [ ! -f $setting ]; then
        printf >&2 "\e[1;31mInvalid setting file:\e[m\n"
        printf >&2 "$HPATH/${case%.*}.sh\n"
        printf >&2 "\tNo setting file found for the case.\n"
        exit 1
    fi
    cp -f $setting $log/job_script.sh

    # Run the submission based on which cluster we attempt to use.
    cd $LPATH/$example
    if [ $CLUSTER == "DTU" ]; then
        export BSUB_QUIET=Y
        bsub -J $1 -env "all" <job_script.sh

    elif [ $CLUSTER == "M5" ]; then
        if [ -z "$M5_ACCOUNT" ]; then
            printf >&2 "No account specified for Marenostrum5.\n"
            printf >&2 "Please set the M5_ACCOUNT variable in the environment.\n"
            exit 1
        fi
        sbatch -A $M5_ACCOUNT -J $1 job_script.sh 1>/dev/null 2>error.err

    else
        printf >&2 "No or invalid cluster specified for submission.\n"
        printf >&2 "\t- DTU for the DTU cluster.\n"
        printf >&2 "\t- M5 for the Marenostrum5 cluster.\n"
        exit 1
    fi

    printf '\t%-12s %-s\n' "Submitted:" "$1"
    cd $CURRENT_DIR
}

# Definition of a interrupt handler
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
set +e # Do not exit on error during execution
full_start=$(date +%s.%N)
QUEUE=""

printf "\n\e[4mQueueing case files.\e[0m\n"
for case in ${example_list[@]}; do

    case_name=$(basename ${case%.*})
    case_dir=$(dirname $case)

    # Define the name of the current exampel, if there are multiple cases in the
    # same folder, we add the case name to the example name.
    example=$case_dir
    if [[ ${case: -5} == ".case" &&
        $(find $EPATH/$case_dir -name "*.case" | wc -l) > 1 ]]; then
        example=$example/$case_name
    fi

    log=$LPATH/$example && mkdir -p $log
    [ "$CLEAN" = true ] && rm -fr $log/*

    # Setup the log folder
    if [[ -f "$log/output.log" &&
        "$(head -n 1 $log/output.log)" == "Ready" ]]; then
        rm -f $log/error.err && touch $log/error.err

        [ ! -z "$CLUSTER" ] && printf '\t%-12s %-s\n' "Queued:" "$example"
        QUEUE="$QUEUE $example"
        continue
    fi

    # Remove old output and error files
    find $log -type f -name "*.log" -or -name "error.err" -delete
    touch $log/output.log $log/error.err

    # Copy the case files to the log folder
    if [ ${case: -3} == ".sh" ]; then
        find $EPATH/$case_dir -name "*.case" -exec cp -ft $log {} +
    elif [ ${case: -5} == ".case" ]; then
        cp -ft $log $EPATH/$case
    fi

    # Copy all data from the case folder to the log folder
    find $EPATH/$case_dir/* -maxdepth 0 \
        -not -name "*.case" -and -not -name "*.nmsh" \
        -exec rsync -r {} $log \;

    # Create symbolic links to the mesh files to avoid copying massive files
    for file in $(find $EPATH/$case_dir -name "*.nmsh" 2>/dev/null); do
        ln -fs $file $log
    done

    # Copy the job script to the log folder
    cp -f $SPATH/functions.sh $log/functions.sh

    # Assign links to the data folders
    if [ -d "$DPATH" ]; then ln -fs $DPATH $log; fi
    if [ -d "$DLPATH" ]; then ln -fs $DLPATH $log; fi

    # Indicate that the case is ready to be run
    printf 'Ready' >$log/output.log

    QUEUE="$QUEUE $example"
    [ -z "$CLUSTER" ] && printf '\t%-12s %-s\n' "Queued:" "$example"
done

# Done with the setup
# ============================================================================ #
# Move to the directory submit or run the code and return

# If we are just doing a dry-run, we exit here
if [ "$DRY" == true ]; then
    $MAIN_DIR/status.sh
    exit 0
fi

for example in $QUEUE; do

    # Move to the log folder and submit the job
    if [ $INTERRUPTED == 1 ]; then
        continue
    elif [ ! -z "$CLUSTER" ]; then
        Submit $example
    else
        Run $example
    fi
done

if [ -z "$CLUSTER" ]; then
    $MAIN_DIR/status.sh
fi

printf "\n"
# # EOF # #
