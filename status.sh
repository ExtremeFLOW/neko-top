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

# ============================================================================ #
# Progress of running examples
#
# Neko writes its log through a buffer, so the last line is often only partly
# written. Only complete step banners are used, i.e. lines like
#   Step =       1234 t =   1.234000000000000E+00             [  10.28% ]
# The logs are read backwards and every search stops at its first match.

# Print the percentage of the last complete step banner of the current solve,
# without the "%". With "adjoint" as second argument only adjoint steps count,
# skipping the forward steps the adjoint re-runs from its checkpoints.
function step_progress() {
    tac "$1" 2>/dev/null | LC_ALL=C awk -v adjoint="$2" '
        /[|] Adjoint / { seen = 1 }
        /Starting (simulation|iteration)/ { exit }
        /Step = +[0-9]+ +t = +[-+0-9.Ee]+ +[[] *[-0-9.]+% []]/ {
            if (adjoint != "" && !seen) next
            sub(/.*[[] */, ""); sub(/%.*/, ""); print; exit
        }'
}

# Print the progress of a running example. A flow solve gives the percentage
# of its last step, e.g. "[ 73.63% ]". An optimization gives the optimizer
# iteration, the phase (F: forward, A: adjoint) and the completed fraction of
# that phase, e.g. "[  12 A  26.37% ]"; this needs the default log level
# (NEKO_LOG_LEVEL >= 1).
# Usage: example_progress <log file> <case file> <output.log>
function example_progress() {
    local log=$1 case_file=$2 out=$3
    local iter marker phase=F pct unsteady=0

    if ! grep -Eq '"optimization"[[:space:]]*:' "$case_file"; then
        pct=$(step_progress "$log")
        [ -n "$pct" ] && printf '[%7s ]' "$pct%"
        return
    fi

    # The iteration is printed to stdout: "Starting iteration N", or
    # "Iteration : N" when restarting from a checkpoint. Iteration 0 is the
    # evaluation of the initial design. Two matches cover a partly written
    # last line.
    iter=$(tac "$out" 2>/dev/null |
        LC_ALL=C grep -a -m 2 -F -e 'Starting iteration' -e 'Iteration      :' |
        LC_ALL=C awk '/(Starting iteration|Iteration +:)[ \t]+[0-9]+[ \t]*$/ {
            print $NF; exit }')

    # The phase is that of the last solve started, "Starting simulation" or
    # "Adjoint Starting simulation"; the forward steps the adjoint re-runs
    # from its checkpoints print no such header. Only adjoint steps print
    # "| Adjoint" solver rows, so one of those ends the search early. With the
    # log in output.log, "Starting iteration" means the iteration's forward
    # solve has not started yet.
    marker=$(tac "$log" 2>/dev/null |
        LC_ALL=C grep -a -m 1 -F -e '| Adjoint ' -e 'Starting simulation' \
            -e 'Starting iteration')
    if [ -z "$marker" ]; then
        # Nothing solved yet, or an optimization without a flow solve.
        [ -n "$iter" ] && printf '[%4d ]' "$iter"
        return
    fi

    if [[ "$marker" == *Adjoint* ]]; then
        phase=A
        pct=$(step_progress "$log" adjoint)
    else
        pct=$(step_progress "$log")
    fi

    # output.log is written at once, a separate log through a buffer: an
    # iteration announced after the last flush has not shown its forward yet.
    if [[ "$phase" == "A" && "$out" -nt "$log" ]]; then
        phase=F pct=""
    fi

    # An unsteady adjoint runs backwards in time and its banner shows the time
    # of the forward solve, so its percentage counts down.
    grep -Eq '"unsteady"[[:space:]]*:[[:space:]]*true' "$case_file" && unsteady=1

    LC_ALL=C awk -v i="${iter:-0}" -v ph="$phase" -v p="$pct" -v u="$unsteady" '
        BEGIN {
            if (p == "") p = 0
            else if (ph == "A" && u) p = 100 - p
            p = p < 0 ? 0 : p > 100 ? 100 : p
            printf "[%4d %s %6.2f%% ]", i, ph, p
        }'
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

                # Without NEKO_LOG_FILE, Neko writes its log to stdout, which
                # ends up in output.log next to the iteration lines.
                if [ ! -f $logfile ] && grep -q -F '(version: ' $LPATH/$test/output.log; then
                    logfile=$LPATH/$test/output.log
                fi

                if [ ! -f $logfile ]; then
                    printf '\t\e[1;33m%-12s\e[m %s %-s\n' "Starting:" "$test"
                    continue
                fi

                if [ "$(tail -n 1 $logfile | xargs)" == "Normal end." ]; then
                    stat="Complete:"
                else
                    stat="Running:"
                    progress=$(example_progress "$logfile" "$f" "$LPATH/$test/output.log")
                fi
                printf '\t\e[1;33m%-12s\e[m' "$stat"
                if [[ "$stat" == "Running:" && ! -z "$progress" ]]; then
                    printf ' %s' "$progress"
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
