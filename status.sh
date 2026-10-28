#!/bin/bash

# ============================================================================ #
# User defined inputs.

# Define all needed folders relative to the Testing folder. (without trailing /)
RPATH="results" # Result export location
LPATH="logs"    # Logging locations

# End of user inputs
# ============================================================================ #
# Help commands
function help() {
    echo "status.sh"
    echo "  This script prints the status of the tests either currently running"
    echo "  or completed. Any completed test which have exited with an error"
    echo "  gets the error printed as well."
}

# Variable setups
MAIN_DIR=$(dirname $(realpath $0))

RPATH=$MAIN_DIR/$RPATH
LPATH=$MAIN_DIR/$LPATH

[ ! -d $LPATH ] && exit 0
EXIT_STATUS=0

# ============================================================================ #
# Keywords

# Print help
for in in $@; do
    if [ "$in" = "-h" ] || [ "$in" = "--help" ]; then
        help && exit
    fi
done

# ============================================================================ #
# Print status

# List all the tests, if there are none we return
tests=($(find -L $LPATH -type d -exec test -f '{}'/output.log \; -print | sort -u))
for ((i = 0; i < ${#tests[@]}; i++)); do tests[$i]="${tests[$i]#$LPATH/}"; done

# Trim tests called `run_*` which are not actual tests
filtered_tests=()
for test in ${tests[@]}; do
    if [[ ! $(basename $test) == run_* ]]; then
        filtered_tests+=("$test")
    fi
done
tests=("${filtered_tests[@]}")

if [ ${#tests[@]} -eq 0 ]; then
    printf "No tests found.\n"
    exit 0
fi

# If we are running in LSF-10 mode, print the running jobs.
if [ $(which bjobs 2>/dev/null) ]; then
    printf "\n\e[4mRunning jobs.\e[m\n"
    bjobs -ro -noheader "time_left:8 job_name"
elif [ $(which squeue 2>/dev/null) ]; then
    printf "\n\e[4mRunning jobs.\e[m\n"
    squeue -rho "%.10t %9L %j" -u $USER
fi

printf "\n\e[4mTest status.\e[m\n"

for test in ${tests[@]}; do
    # Check if there were errors. Print them if there were.
    if [ -s $LPATH/$test/error.log ]; then

        if [ "$(head -n 1 $LPATH/$test/error.log)" = "Interrupted" ]; then
            printf '\t\e[1;31m%-12s\e[m %-s\n' "Interrupted:" "$test"
        else
            printf '\t\e[1;31m%-12s\e[m %-s\n' "Error:" "$test"
            EXIT_STATUS=1
        fi
    elif [[ -s $LPATH/$test/output.log && ! -s $LPATH/$test/error.log ]]; then
        file=($(find $LPATH/$test -type f -name "*.case"))

        # If more than one file exists
        if [[ ${#file[@]} -gt 2 ]]; then
            file+="$LPATH/$test/output.log"
        fi

        if [ "$(head -n 1 $LPATH/$test/output.log)" = "Ready" ]; then
            printf '\t\e[1;33m%-12s\e[m %s %-s\n' "Pending:" "$test"
        elif [ ${#file} -gt 0 ]; then
            for f in ${file[@]}; do
                logfile=${f%.*}.log

                if [ ! -f $logfile ]; then
                    printf '\t\e[1;33m%-12s\e[m %s %-s\n' "Starting:" "$test"
                    continue
                fi

                if [ "$(tail -n 1 ${f%.*}.log | xargs)" == "Normal end." ]; then
                    stat="Complete:"
                else
                    stat="Running:"
                    progress=$(
                        tail -n 100 "${f%.*}.log" |        # Get the last 1000 lines
                            grep '^\s*Step = ' |           # Get all timestamps
                            tail -n 1 |                    # Get the last line
                            sed -e 's/.*\[\(.*\)].*/\1/' | # Get the progress
                            xargs                          # Trim whitespace
                    )
                fi
                printf '\t\e[1;33m%-12s\e[m' "$stat"
                if [[ "$stat" == "Running:" && ! -z "$progress" ]]; then
                    printf ' [%7s ]' "$progress"
                fi

                if [[ $(basename $f) == "output.log" || ${#file[@]} -lt 2 ]]; then
                    printf " %s\n" "$test"
                else
                    printf " %s\n" "$test/$(basename $f)"
                fi
            done
        else
            printf '\t\e[1;33m%-12s\e[m %s %-s\n' "Starting:" "$test"
        fi

    elif [[ -d $RPATH/$test && ! -s $LPATH/$test/output.log && ! -s $LPATH/$test/error.log ]]; then
        printf '\t\e[1;32m%-12s\e[m %-s\n' "Complete:" "$test"
    fi
done

# ============================================================================ #
# Print errors for all unfinished tests

for test in ${tests[@]}; do
    # Check if there were errors. Print them if there were.
    if [ -s $LPATH/$test/error.log ]; then

        if [ "$(head -n 1 $LPATH/$test/error.log)" = "Interrupted" ]; then
            continue
        fi

        # Print the header with example name
        printf '\n\e[4;31m%-s\e[m' "${test:0:79}"
        printf '\e[4;31m%.0s_\e[m' $(seq 1 $((80 - ${#test}))) && printf '\n'

        # Print the error messagge, sitting between *** ERROR and ERROR STOP
        start_line=$(grep -n '\*\*\* ERROR' $LPATH/$test/error.log | head -n 1 | cut -d: -f1)
        end_line=$(grep -n 'STOP' $LPATH/$test/error.log | head -n 1 | cut -d: -f1)
        end_line=$((end_line - 1)) # Remove the line with ERROR STOP

        if [ -z "$start_line" ] || [ -z "$end_line" ]; then
            head -n 20 $LPATH/$test/error.log
        else
            sed -n "${start_line},${end_line}p" $LPATH/$test/error.log
        fi

    fi
done
printf "\n"

# Remove all empty folders in the logs folder
find $LPATH -type d -empty -delete

exit $EXIT_STATUS
# # EOF # #
