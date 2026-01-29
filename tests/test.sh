#!/bin/bash
set -e # Exit with non-zero exit code if anything fails
# ============================================================================ #
# Define the help function

function help() {
    echo -e "test.sh [options] TEST_NAME [TEST_NAME ...]"
    echo -e "  Script to execute tests."
    echo -e ""
    echo -e " TEST_NAME options:"
    echo -e "  unit           Run unit tests."
    echo -e "  sensitivity    Run sensitivity regression tests."
    echo -e "  mma            Run MMA regression tests."
    echo -e ""
    echo -e " Options:"
    echo -e "  -h, --help      Show this help message and exit."
    echo -e "  -a, --all       Run all tests (unit and sensitivity)."
    echo -e ""
    echo -e "  See Readme for additional details."
    exit 0
}

# Handle options
OPT="h,a"
OPTIONS="help,all,procs:"
parser=$(getopt --options=$OPT --longoptions=$OPTIONS --name "$0" -- "$@")

ALL=false
eval set -- "$parser"
while true; do
    case "$1" in
        -h|--help) help && shift ;;
        -a|--all) ALL=true && shift ;;
        --procs) NP=$2 && shift 2 ;;
        --) shift && break ;;
        *) echo "Unexpected option: $1" && help ;;
    esac
done

if [ $# -eq 0 ] && [ "$ALL" == false ]; then help; fi

# ============================================================================ #
# Determine which tests to run

UNIT_TEST=false
SENSITIVITY_TEST=false
MMA_TEST=false
REGRESSION_TEST=false

if [ "$ALL" == true ]; then
    UNIT_TEST=true
    SENSITIVITY_TEST=true
    MMA_TEST=true
    REGRESSION_TEST=true
else
    for arg in "$@"; do
        case "$arg" in
            unit) UNIT_TEST=true ;;
            regression) REGRESSION_TEST=true ;;
            sensitivity) SENSITIVITY_TEST=true ;;
            mma) MMA_TEST=true ;;
            *) echo "Invalid argument: $arg" && help ;;
        esac
    done
fi

# ============================================================================ #
# Set up environment variables and source dependencies

MAIN_DIR=$(realpath $(dirname $0)/../)

[ -f "$MAIN_DIR/prepare.env" ] && source "$MAIN_DIR/prepare.env"
[ -z "$EXTERNAL_DIR" ] && export EXTERNAL_DIR="$MAIN_DIR/external"
source "$MAIN_DIR/scripts/dependencies.sh"
find_neko

set +e # Continue execution even if individual tests fail
# ============================================================================ #
# Run the unit tests

if [ "$UNIT_TEST" == true ]; then
    ctest -L unit -C Debug -O test_report.log --test-dir $MAIN_DIR/build \
        --output-on-failure -j ${NP:-1} --no-tests=ignore

    if [ $? -ne 0 ]; then
        UNIT_SUCCESS=Failure
    else
        UNIT_SUCCESS=Success
    fi
else
    UNIT_SUCCESS=Skipped
fi

# ============================================================================ #
# Run the regression tests

if [ "$REGRESSION_TEST" == true ]; then
    ctest -L regression -C Debug -O test_report.log --test-dir $MAIN_DIR/build \
        --output-on-failure -j ${NP:-1} --no-tests=ignore

    if [ $? -ne 0 ]; then
        REGRESSION_SUCCESS=Failure
    else
        REGRESSION_SUCCESS=Success
    fi
else
    REGRESSION_SUCCESS=Skipped
fi

# ============================================================================ #
# Run the sensitivity regression tests

if [ "$SENSITIVITY_TEST" == true ]; then
    $MAIN_DIR/tests/regression/sensitivity/run.sh -np=$NP

    if [ $? -ne 0 ]; then
        SENSITIVITY_SUCCESS=Failure
    else
        SENSITIVITY_SUCCESS=Success
    fi
else
    SENSITIVITY_SUCCESS=Skipped
fi

# ============================================================================ #
# Run the mma regression tests

if [ "$MMA_TEST" == true ]; then
    $MAIN_DIR/tests/regression/mma/run.sh -np=$NP

    if [ $? -ne 0 ]; then
        MMA_SUCCESS=Failure
    else
        MMA_SUCCESS=Success
    fi
else
    MMA_SUCCESS=Skipped
fi

# ============================================================================ #
# Summary of test results

printf "\nTest Summary:\n"
printf "=========================\n"
printf "%-15s %s\n" "Unit:"        "${UNIT_SUCCESS}"
printf "%-15s %s\n" "Regression:"  "${REGRESSION_SUCCESS}"
printf "%-15s %s\n" "Sensitivity:" "${SENSITIVITY_SUCCESS}"
printf "%-15s %s\n" "MMA:"         "${MMA_SUCCESS}"

# Exit with failure if any test failed
# Exit non-zero if any *_SUCCESS variable is "Failure"
while IFS= read -r name; do
    if [[ "${!name}" == "Failure" ]]; then
        printf "\nERROR: One or more tests failed.\n" >&2
        exit 1
    fi
done < <(compgen -v | grep '_SUCCESS$' || true)

exit 0
