#!/bin/bash

set -e

# ============================================================================ #
# User defined inputs.

# Define all needed folders relative to the project folder. (without trailing /)
EPATH="Examples" # Examples scripts
RPATH="Results"  # Result export location
LPATH="logs"     # Logging locations
HPATH="Scripts"  # Submition settings

# End of user inputs
# ============================================================================ #
# Help commands
function help() {
    echo "run.sh EXAMPLE [EXAMPLE ...]"
    echo "  This script works as a function which run the examples specified"
    echo "  through the command line."
    echo ""
    echo "  The <EXAMPLE> input refers to the name of the subfolder located in"
    echo "  $EPATH. This folder should contain everything needed to run the "
    echo "  example."
    echo ""
    echo "  See Readme for details."
    echo ""
    echo "Available examples;"
    for example in $(find $EPATH -d 2>/dev/null); do
        if [ -f $example/*.case ]; then
            t_name=${example#$EPATH/}
            printf '  - %-s\n' $t_name
        fi
    done
}

# Variable setups
MAINDIR=$PWD

export EPATH=$MAINDIR/$EPATH
export RPATH=$MAINDIR/$RPATH
export LPATH=$MAINDIR/$LPATH
export HPATH=$MAINDIR/$HPATH

mkdir -p $RPATH $LPATH
# ============================================================================ #
# Keywords

# Empty input
if [ -z "$1" ]; then
    help && exit
fi

examples=""

for in in $@; do

    # Opions and flags
    if [[ ${in::1} == "-" ]]; then
        for ((i = 1; i < ${#in}; i++)); do
            case "${in:$i:1}" in
            "a" | "-all") ALL=true ;;       # Run all examples available
            "c" | "-clean") CLEAN=true ;;   # Clean logs
            "h" | "-help") help && exit ;;  # Print help
            "d" | "-delete") DELETE=true ;; # Delete previous runs
            esac
        done

    # Ignore invalid inputs
    elif [[ -z $(ls $EPATH/$in 2>/dev/null) ]]; then
        printf '  %-10s %-67s\n' "Not Found:" "$in"

    # Extract the examples from the input
    elif [[ ! $ALL ]]; then
        for example in $(find $EPATH/$in -type d); do

            case $(find $example -maxdepth 1 -name "*.case" | wc -l) in
            0) ;;
            1) examples+="${example#$EPATH/} " ;;
            *) printf "Multiple case files in: ${example#$EPATH/}\n" >&2 && exit 1 ;;
            esac

        done
    fi
done

if [ $ALL ]; then
    examples=""
    for example in $(find $EPATH/ -type d); do
        case $(find $example -maxdepth 1 -name "*.case" | wc -l) in
        0) ;;
        1) examples+="${example#$EPATH/} " ;;
        *) printf "Multiple case files in: ${example#$EPATH/}\n" >&2 && exit 1 ;;
        esac
    done
fi

if [ ! "$examples" ]; then
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
    example=$1
    cd $LPATH/$example

    echo "Launching example on local machine" 1>output.out

    # Run the example and pipe stdout and stderr to the log files
    ./job_script.sh $example >output.out 2>error.err &
    printf '  %-13s %-64s\n' "Started:" "${example:0:64}"
    wait

    cd $MAINDIR
}

# Function for submitting the examples
function Submit() {
    example=$1
    cd $LPATH/$example

    touch output.out
    echo "Launching example on HPC" 1>output.out

    export BSUB_QUIET=Y
    bsub -J $example -env "RPATH=$RPATH" <./job_script.sh
    printf '  %-13s %-64s\n' "Submitted:" "${example:0:64}"

    cd $MAINDIR
}

# ============================================================================ #
# Run the examples
full_start=$(date +%s.%N)

printf "\n\e[4mQueueing examples.\e[0m\n"
for example in $examples; do
    if [ ! -f $EPATH/$example/neko ]; then
        printf '  %-13s %-64s\n' "Not Compiled:" "$example"
        continue
    fi

    # Create the folder needed
    log=$LPATH/$example

    # Cleanup old log files
    if [ $CLEAN ]; then
        rm -fr $log
    fi
    mkdir -p $log
    cd $log

    # Remove old output and error files
    rm -f output.out error.err
    touch output.out error.err

    # Determine if we have a HPC file
    setting=$HPATH/${example}.sh
    if [ ! -f $setting ]; then setting=$HPATH/${example%/*}.sh; fi
    if [ ! -f $setting ]; then setting=$HPATH/default.sh; fi

    cp -rt . $EPATH/$example/neko $EPATH/$example/*.case $EPATH/$example/*.nmsh
    cp $setting ./job_script.sh

    # Done with the setup
    # ======================================================================== #
    # Move to the directory submit the code and return

    if [ "$(which bsub)" ]; then
        Submit $example
    else
        Run $example
    fi

done

if [ ! "$(which bsub)" ]; then
    ./status.sh
fi

printf "\n"
# # EOF # #
