#!/bin/bash

# ============================================================================ #
# Define the help function

function help() {
    echo -e "test.sh [options]"
    echo -e "  Script to execute the tests."
    echo -e ""
    echo -e " Options:"
    echo -e "  -h, --help      Show this help message and exit."
    echo -e ""
    echo -e "  See Readme for additional details."
    exit 0
}

# Handle options
OPT="h,a"
OPTIONS="help,all"
parser=$(getopt --options=$OPT --longoptions=$OPTIONS --name "$0" -- "$@")

ALL=false
eval set -- "$parser"
while true; do
    case "$1" in
        -h|--help) help && shift ;;
        -a|--all) ALL=true && shift ;;
        --) shift && break ;;
        *) echo "Unexpected option: $1" && help ;;
    esac
done

[ "$ALL" == true ] && UNIT_TEST=true || UNIT_TEST=false
[ "$ALL" == true ] && SENSITIVITY_TEST=true || SENSITIVITY_TEST=false

for arg in "$@"; do
    if [ "$arg" == "unit" ]; then
        UNIT_TEST=true
    elif [ "$arg" == "sensitivity" ]; then
        SENSITIVITY_TEST=true
    else
        echo "Invalid argument: $arg"
        help
    fi
done

# ============================================================================ #
# Set up environment variables and source dependencies

MAIN_DIR=$(realpath $(dirname $0)/../)

EXTERNAL_DIR="$MAIN_DIR/external"
source "$MAIN_DIR/scripts/dependencies.sh"
find_hdf5
find_neko

# ============================================================================ #
# Run the unit tests

if [ "$UNIT_TEST" == true ]; then
    ctest -C Debug -O test_report.log --verbose --test-dir $MAIN_DIR/build

    if [ $? -ne 0 ]; then
        echo "Some tests failed. Check test_report.log for details."
        exit 1
    else
        echo "Unit tests passed successfully."
        exit 0
    fi
fi

# ============================================================================ #
# Run the sensitivity regression tests

if [ "$SENSITIVITY_TEST" == true ]; then
    $MAIN_DIR/tests/regression/sensitivity/run.sh

    if [ $? -ne 0 ]; then
        echo "Sensitivity regression tests failed."
        exit 1
    else
        echo "Regression tests passed successfully."
        exit 0
    fi
fi
