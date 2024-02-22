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

if [[ ! -d $RPATH && ! -d $LPATH ]]; then exit 0; fi

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
for test in $(find $LPATH -type d 2>/dev/null); do
    if [ -f $test/output.out ]; then
        tests+="${test#$LPATH/} "
    fi
done

if [ -z "$tests" ]; then exit; fi

# If we are running in LSF-10 mode, print the running jobs.
if [ "$(which bsub)" ]; then
    printf "\n\e[4mRunning jobs.\e[m\n"
    bjobs -ro -noheader "time_left:8 job_name"
fi

printf "\n\e[4mTest status.\e[m\n"

for test in $tests; do
    if [[ -d $RPATH/$test && ! -s $LPATH/$test/output.out && ! -s $LPATH/$test/error.err ]]; then
        printf '  \e[1;32m%-10s\e[m %-s\n' "Complete:" "$test"
        rm -fr $LPATH/$test
    fi
done

for test in $tests; do
    if [[ -f $LPATH/$test/neko && ! -s $LPATH/$test/error.err ]]; then
        file=$(find $LPATH/$test -type f -name "*.case")
        if [ -f "${file%.*}.out" ]; then
            progress=$(
                grep 't = ' "${file%.*}.out" |     # Get all lines with t = in them
                    tail -n 1 |                    # Get the last line
                    sed -e 's/.*\[\(.*\)].*/\1/' | # Get the progress
                    xargs                          # Trim whitespace
            )
            printf '  \e[1;33m%-10s\e[m [ %6s ] %s %-s\n' "Running:" "$progress" "$test"
        else
            printf '  \e[1;33m%-10s\e[m %s %-s\n' "Pending:" "$test"
        fi
    fi
done

for test in $tests; do
    # Check if there were errors. Print them if there were.
    if [ -s $LPATH/$test/error.err ]; then
        printf '  \e[1;31m%-10s\e[m %-s\n' "Error:" "$test"
    fi
done

# ============================================================================ #
# Print errors for all unfinished tests

for test in $tests; do
    # Check if there were errors. Print them if there were.
    if [ -s $LPATH/$test/error.err ]; then
        printf '\n\e[4;31m%-s\e[m' "${test:0:79}"
        printf '\e[4;31m%.0s_\e[m' $(seq 1 $((80 - ${#test}))) && printf '\n'
        head -n 10 $LPATH/$test/error.err | fold -w 80
    fi
done
printf "\n"

# Remove all empty folders in the logs folder
find $LPATH -type d -empty -delete
# # EOF # #
