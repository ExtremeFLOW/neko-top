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
    echo -e "  Please note: Cases are in fact json files. We support regular"
    echo -e "  json files, json files hidden under the '.case' filename."
    echo -e ""
    echo -e "  See Readme for additional details."
    echo -e ""

    printf "\e[4mOptions:\e[0m\n"
    printf "  -%-1s, --%-10s %-60s\n" "a" "all" "Run all journals available."
    printf "  -%-1s, --%-10s %-60s\n" "c" "clean" "Clean artifacts from previous runs."
    printf "  -%-1s, --%-10s %-60s\n" "d" "delete" "Delete previously completed runs."
    printf "  -%-1s, --%-10s %-60s\n" "h" "help" "Print help."
    printf "  -%-1s, --%-10s %-60s\n" "n" "neko" "Look for examples in neko."
    printf "  -%-1s, --%-10s %-60s\n" "s" "submit" "Submit the examples to a cluster."
    printf "  -%-1s, --%-10s %-60s\n" " " "dry-run" "Dry run the script."
    printf "  -%-1s, --%-10s %-60s\n" "r" "re-run" "Re-run the examples."
    printf "  -%-1s, --%-10s %-60s\n" "p" "procs" "Number of processors to use."
    printf "  -%-1s, --%-10s %-60s\n" " " "sequential" "Submit the examples sequentially."
    printf "  -%-1s, --%-10s %-60s\n" " " "njobs" "Number of jobs to submit per example."

    printf "\n\e[4mEnvironment:\e[0m\n"
    printf "  -%-1s %-60s\n" "NEKO_DIR" "Path to the Neko installation."
    printf "  -%-1s %-60s\n" "NPROCS" "Number of processors to use."

    printf "\n\e[4mAvailable case files:\e[0m\n"
    for case in $(find $EPATH -name "*.case" 2>/dev/null); do
        printf '  - %-s\n' ${case#$EPATH/}
    done
}
if [ $# -lt 1 ]; then help; fi

# ============================================================================ #
# Define the main directory

export CURRENT_DIR=$(pwd)
export MAIN_DIR=$(dirname $(realpath $0))
export EXTERNAL_DIR="$MAIN_DIR/external"

# ============================================================================ #
# User defined inputs.

# Assign default values to the options
ALL=false
CLEAN=false
NEKO=false
DELETE=false
CLUSTER=""
SEQUENTIAL=false
DRY=false
RERUN=false
N_JOBS=1

# List possible options
OPTIONS=all,clean,help,neko,delete,submit:,dry-run,re-run,sequential,procs:,njobs:
OPT=a,c,h,n,s:,d,r,p:

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
    "-s" | "--submit") CLUSTER="${2^^}" && shift 2 ;; # Submit to the queue
    "-r" | "--re-run") RERUN=true && shift ;;     # Re-run the examples
    "-p" | "--procs") NPROCS="$2" && shift 2 ;;    # Number of processors to use

    # Long option with no short option
    "--dry-run") DRY=true && shift ;;             # Dry run
    "--sequential") SEQUENTIAL=true && shift ;;   # Submit sequentially
    "--njobs") N_JOBS="$2" && shift 2 ;;          # Number of jobs to submit per example

    # End of options
    "--") shift && break ;;
    esac
done

# Fix case of the cluster name
CLUSTER=$(echo $CLUSTER | awk '{print toupper($0)}')

# ============================================================================ #
# Define environment

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
if [ -n "$CLUSTER" ]; then
    export HPATH="$MAIN_DIR/scripts/jobscripts/$CLUSTER" # Submission settings
    if [ ! -d "$HPATH" ]; then
        printf >&2 "\e[1;31mInvalid Cluster:\e[m $CLUSTER\n"
        exit 1
    fi
fi

[ -z "$NEKO_DIR" ] && export NEKO_DIR="$MAIN_DIR/external/neko"
export NEKO_DIR=$(realpath $NEKO_DIR)

if [ "$NEKO" == true ]; then
    export EPATH="$NEKO_DIR/examples"
    export RPATH="$RPATH/neko"
    export LPATH="$LPATH/neko"
fi

# End of user inputs
# ============================================================================ #
# Find the examples to run
set +e # Do not exit on error

example_list=()
for in in $@; do
    [ "$ALL" == true ] && break

    # Decompose the input into the directory and the base name
    [ $(dirname $in) == "." ] && dir="" || dir=$(dirname $in)
    base=$(basename $in)

    # Extract the examples from the input
    matches=($(find $EPATH/$dir -mindepth 1 -maxdepth 1 -type d -name "$base"))
    matches+=($(find $EPATH/$dir -mindepth 1 -maxdepth 1 -type f -name "$base"))
    matches+=($(find $EPATH/$dir -mindepth 1 -maxdepth 1 -type f -name "$base.case"))
    matches+=($(find $EPATH/$dir -mindepth 1 -maxdepth 1 -type f -name "$base.json"))

    for match in ${matches[@]}; do
        file_list=()
        if [ -d $match ]; then
            file_list=($(find $match -name "run.sh" 2>/dev/null))
            file_list+=($(find $match -name "*.case" 2>/dev/null))
            file_list+=($(find $match -name "*.json" 2>/dev/null))
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
    file_list+=($(find $EPATH -name "*.json" 2>/dev/null))

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
    while [ $parent != "." ]; do

        if [[ -n "$(find $EPATH/$parent -maxdepth 1 -name '*.case' -or -name '*.json')" ]]; then

            printf >&2 "\e[1;31mInvalid example file:\e[m\n"
            printf >&2 "$EPATH/$example\n"
            printf >&2 "\tNested examples are not allowed.\n"
            printf >&2 "\tMove the $example file to the root of example suite\n"
            if [[ ${example: -5} == ".case" || ${example: -5} == ".json" ]]; then
                printf >&2 "\tor create a run.sh file in the parent folder.\n"
            fi

            unset example_list[$i]

            parent="."
        else
            parent=$(dirname $parent)
        fi
    done

done

# Remove duplicates and check for nested examples
example_list=($(echo "${example_list[@]}" | tr ' ' '\n' | sort -u))

# If multiple examples with same name and  different file extensions are found
# we stop the execution and print an error message.
for example in ${example_list[@]}; do

    matches=($(
        find $EPATH -wholename "$EPATH/${example%.*}.json" \
            -or -wholename "$EPATH/${example%.*}.case"
    ))

    if [ ${#matches[@]} -gt 1 ]; then
        printf >&2 "\e[1;31mInvalid example file:\e[m ${example%.*}\n"
        printf >&2 "\tMultiple examples with the same name found.\n"
        printf >&2 "\tPlease remove the duplicates.\n"
        for match in ${matches[@]}; do
            printf >&2 "\t- ${match#$EPATH/}\n"
        done
        exit 1
    fi
done

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
    source functions.sh
    run $1 1>output.log 2>error.log
    cd $CURRENT_DIR
}

function ExampleUsesPodStateRecovery() {
    local example_dir=$1
    local case_file

    while IFS= read -r case_file; do
        if grep -q '"type"[[:space:]]*:[[:space:]]*"pod"' "$case_file"; then
            return 0
        fi
    done < <(
        find "$example_dir" -maxdepth 1 -type f \
            \( -name "*.case" -or -name "*.json" \) | sort
    )

    return 1
}

# Function for submitting the examples
function Submit() {

    # Run the submission based on which cluster we attempt to use.
    cd $LPATH/$example

    if [ $CLUSTER == "MN5" ]; then
        if [ -z "$MN5_ACCOUNT" ]; then
            printf >&2 "No account specified for Marenostrum5.\n"
            printf >&2 "Using SLURM environment variables if available\n"
            printf >&2 "Assign the 'MN5_ACCOUNT' environment variable to avoid"
            printf >&2 "this message."
        else
            export SBATCH_ACCOUNT="$MN5_ACCOUNT"
        fi

    elif [[ $CLUSTER == "LUMI-C" || $CLUSTER == "LUMI-G" ]]; then
        if [ -z "$LUMI_ACCOUNT" ]; then
            printf >&2 "No account specified for LUMI.\n"
            printf >&2 "Using SLURM environment variables if available\n"
            printf >&2 "Assign the 'LUMI_ACCOUNT' environment variable to avoid"
            printf >&2 "this message."
        else
            export SBATCH_ACCOUNT="$LUMI_ACCOUNT"
        fi
    fi

    if [ -n "$(which bsub 2>/dev/null)" ]; then
        bsub -J $1 -env "all" <job_script.sh
    elif [ -n "$(which sbatch 2>/dev/null)" ]; then
        if [ "$(squeue -h --name=$1 --me | wc -l)" -gt 0 ]; then
            printf '\t%-12s %-s\n' "In queue:" "$1"
            cd $CURRENT_DIR
            return
        fi

        # Deal with sequential submission and job dependencies
        id=""
        DEP=""
        for i in $(seq 1 $N_JOBS); do
            if [[ -n "$id" && -n "$SEQ_DEP" ]]; then
                DEP="--dependency=afterany:$id:$SEQ_DEP"
            elif [ -n "$id" ]; then
                DEP="--dependency=afterany:$id"
            elif [ -n "$SEQ_DEP" ]; then
                DEP="--dependency=afterany:$SEQ_DEP"
            fi

            id=$(sbatch --parsable -J $1 $DEP job_script.sh)
            if [ "$SEQUENTIAL" == true ]; then
                SEQ_DEP="$SEQ_DEP:$id"
            fi
        done
    else
        printf >&2 "Unknown submission system.\n"
        exit 1
    fi

    printf '\t%-12s %-s\n' "Submitted:" "$1"
    cd $CURRENT_DIR
}

# Definition of a interrupt handler
INTERRUPTED=0
function handler() {
    if [ "$MAIN_DIR" != "$(pwd)" ]; then
        printf "Interrupted" >>error.log
    fi
    INTERRUPTED=1
}
trap 'handler' SIGINT

# ============================================================================ #
# Compile the example executables

if [[ "$NEKO" != true && -d $MAIN_DIR/build ]]; then
    printf "\n\e[4mCompiling the examples.\e[0m\n"
    cmake --build $MAIN_DIR/build --target Examples --parallel

    # Check if the compilation was successful
    if [ $? -ne 0 ]; then
        printf >&2 "\e[1;31mCompilation failed.\e[m\n"
        exit 1
    fi
fi

# ============================================================================ #
# Run the examples
full_start=$(date +%s.%N)
QUEUE=""

printf "\n\e[4mQueueing case files.\e[0m\n"
for case in ${example_list[@]}; do

    case_name=$(basename ${case%.*})
    case_dir=$(dirname $case)

    # Define the name of the current exampel, if there are multiple cases in the
    # same folder, we add the case name to the example name.
    example=$case_dir

    if [[ $(find $EPATH/$case_dir -name "*.case" -or -name "*.json" | wc -l) > 1 ]]; then
        example=$example/$case_name
    fi

    if [ "$RERUN" == false ] && [ -d "$RPATH/$example" ]; then
        printf '\t\e[1;32m%-12s\e[m %-s\n' "Complete:" "$example"
        continue
    fi

    case "$CLUSTER" in
        "MN5" | "LUMI-C" | "LUMI-G")
            if [[ "$(squeue -h --name=$example --me | wc -l)" -gt 0 ]]; then
                printf '\t\e[1;33m%-12s\e[m %s %-s\n' "In queue:" "$example"
                continue
            fi
        ;;
        "") ;;
    esac

    export log=$LPATH/$example
    if [[ "$CLEAN" == true && -d "$log" ]]; then
        rm -fr $log
    fi

    # Setup the log folder
    if [[ -f "$log/output.log" &&
        "$(head -n 1 $log/output.log)" == "Ready" ]]; then
        rm -f $log/error.log && touch $log/error.log

        [ -n "$CLUSTER" ] && printf '\t%-12s %-s\n' "Queued:" "$example"
        QUEUE="$QUEUE $example"
        continue

    elif [[ -s "$log/error.log" ]]; then
        # Move old log files to folder with counter padded to 2 digits
        old_run=run_$(find $log -maxdepth 1 -type d -name "run_*" | wc -l)
        old_run=$(printf "%s_%02d" "run" $((10#${old_run#run_} + 1)))
        mkdir -p $log/$old_run

        find $log -maxdepth 1 -not -empty -type f -name "*.log" \
            -exec mv -ft $log/$old_run {} \;

        touch $log/output.log $log/error.log
        echo "Ready" >$log/output.log

        [ -n "$CLUSTER" ] && printf '\t%-12s %-s\n' "Restarting:" "$example"
        QUEUE="$QUEUE $example"
        continue

    elif [[ -f "$log/output.log" ]]; then
        printf '\t\e[1;33m%-12s\e[m %s %-s\n' "Skipped:" "$example"
        continue
    fi

    mkdir -p $log
    touch $log/output.log $log/error.log

    # Copy the case files to the log folder
    if [ ${case: -6} == "run.sh" ]; then
        find $EPATH/$case_dir \( -name "*.case" -or -name "*.json" \) \
            -exec cp -ft $log {} +
    elif [ ${case: -5} == ".case" ]; then
        cp -ft $log $EPATH/$case
    elif [ ${case: -5} == ".json" ]; then
        cp -ft $log $EPATH/$case
    fi

    # Copy all data from the case folder to the log folder
    find $EPATH/$case_dir/* -maxdepth 0 \
        -not -name "*.case" -and -not -name "*.json" -and -not -name "*.nmsh" \
        -exec rsync -r {} $log \;

    # Create symbolic links to the mesh files to avoid copying massive files
    find $EPATH/$case_dir -name "*.nmsh" -exec ln -fs {} $log \;

    # Copy the job script to the log folder
    cp -f $SPATH/functions.sh $log/functions.sh

    # If we are submitting to a cluster, look for the associated jobscript
    if [ -n "$CLUSTER" ]; then

        setting=$HPATH/${case%.*}.sh
        exact_setting=$setting

        # Find the setting file for the case recursively
        while [[ ! -f $setting && "$(dirname $setting)" != "$HPATH" ]]; do
            setting=$(dirname ${setting%/default.sh})/default.sh
        done

        if [ ! -f "$exact_setting" ] && \
           [ -f "$HPATH/pod/default.sh" ] && \
           ExampleUsesPodStateRecovery "$EPATH/$case_dir"; then
            setting="$HPATH/pod/default.sh"
        fi

        setting=$(realpath $setting)

        if [ ! -f $setting ]; then
            printf >&2 "\e[1;31mInvalid setting file:\e[m\n"
            printf >&2 "$HPATH/${case%.*}.sh\n"
            printf >&2 "\tNo setting file found for the case.\n"
            exit 1
        else
            cp -f $setting $log/job_script.sh
        fi
    fi

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
    exit $?
fi

for example in $QUEUE; do

    # Move to the log folder and submit the job
    if [ $INTERRUPTED == 1 ]; then
        continue
    elif [ -n "$CLUSTER" ]; then
        Submit $example
    else
        Run $example
    fi
done

if [ -z "$CLUSTER" ]; then
    $MAIN_DIR/status.sh
    exit $?
fi

printf "\n"
# # EOF # #
