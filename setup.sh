#!/bin/bash
set -e # Exit with nonzero exit code if anything fails
# ============================================================================ #
# Print the help message
function help() {
    echo -e "Usage: $0 [options]"
    echo -e "Options:"
    echo -e "\t-h, --help        Show this help message and exit"
    echo -e "\t-t, --test        Run the tests after the installation"
    echo -e "\t-c, --clean       Clean the build directory before compiling"
    echo -e "\t-q, --quiet       Suppress output"
    echo -e "\t-d, --device      Device type to compile for (off, CUDA)"
    echo -e ""
    echo -e "Compilation and setup of Neko-TOP, this script will install all"
    echo -e "the dependencies and compile the Neko-TOP code."
    echo -e ""
    echo -e "Environment Variables:"
    echo -e "\tNEKO_DIR          The directory where Neko is installed"
    echo -e "\tJSON_FORTRAN_DIR  The directory where JSON-Fortran is installed"
    echo -e "\tNEK5000_DIR       The directory where Nek5000 is installed"
    echo -e "\tPFUNIT_DIR        The directory where PFUnit is installed"
    echo -e "\tGSLIB_DIR         The directory where GSLIB is installed"
    echo -e "\tCUDA_DIR          The directory where CUDA is installed"
    echo -e "\tBLAS_DIR          The directory where BLAS is installed"
}

# Assign default values to the options
DEVICE_TYPE="OFF"
CLEAN=false
QUIET=false
TEST=false

# List possible options
OPTIONS=help,test,clean,quiet,device:
OPT=h,t,c,q,d:

# Parse the inputs for options
PARSED=$(getopt --options=$OPT --longoptions=$OPTIONS --name "$0" -- "$@")
eval set -- "$PARSED"

# Loop through the options and set the variables
while true; do
    case "$1" in
    "-h" | "--help") help && exit ;;                  # Print help
    "-t" | "--test") TEST=true && shift ;;            # Build the tests
    "-c" | "--clean") CLEAN=true && shift ;;          # Clean compilation
    "-q" | "--quiet") QUIET=true && shift ;;          # Suppress output
    "-d" | "--device") DEVICE_TYPE="$2" && shift 2 ;; # Device type

    # End of options
    "--") shift && break ;;
    esac
done
export TEST CLEAN QUIET DEVICE_TYPE

# ============================================================================ #
# Set main directories

CURRENT_DIR=$(pwd)
MAIN_DIR=$(dirname $(realpath $0))
EXTERNAL_DIR="$MAIN_DIR/external"

# ============================================================================ #
# Execute the preparation script if it exists and prepare the environment

# Execute the preparation script if it exists
if [ -f "$MAIN_DIR/prepare.env" ]; then
    source $MAIN_DIR/prepare.env
fi
source $MAIN_DIR/scripts/dependencies.sh

# Ensure local dependencies are used if they are not defined by the environment
[ -z "$NEKO_DIR" ] && NEKO_DIR="$EXTERNAL_DIR/neko"
[ -z "$JSON_FORTRAN_DIR" ] && JSON_FORTRAN_DIR="$EXTERNAL_DIR/json-fortran"
[ -z "$NEK5000_DIR" ] && NEK5000_DIR="$EXTERNAL_DIR/Nek5000"
[ -z "$GSLIB_DIR" ] && GSLIB_DIR="$NEK5000_DIR/3rd_party/gslib"
[ -z "$PFUNIT_DIR" ] && PFUNIT_DIR="$EXTERNAL_DIR/pFUnit"

# Define standard compilers if they are not defined as environment variables
if [ -z "$CC" ]; then export CC=$(which gcc); else export CC; fi
if [ -z "$CXX" ]; then export CXX=$(which g++); else export CXX; fi
if [ -z "$FC" ]; then export FC=$(which gfortran); else export FC; fi
if [ -z "$NVCC" ]; then export NVCC=$(which nvcc); else export NVCC; fi

# Everything past this point should be general across all setups.
# ============================================================================ #
# Install dependencies (See scripts/dependencies.sh for details)

check_system_dependencies           # Check for system dependencies.
find_json_fortran $JSON_FORTRAN_DIR # Re-defines the JSON_FORTRAN_DIR variable.
find_gslib $GSLIB_DIR               # Re-defines the GSLIB_DIR variable.
find_pfunit $PFUNIT_DIR             # Re-defines the PFUNIT_DIR variable.
find_neko $NEKO_DIR                 # Re-defines the NEKO_DIR variable.

# Done settng up external dependencies
# ============================================================================ #
# Compile the Neko-TOP and example codes.

# Set the variables for the compilation
VARIABLES=("-DJSON_FORTRAN_DIR=$JSON_FORTRAN")
[ "$TEST" == true ] && VARIABLES+=("-DBUILD_TESTING=ON")
[ "$TEST" == true ] && VARIABLES+=("-DPFUNIT_DIR=$PFUNIT_DIR/cmake")
[ "$DEVICE_TYPE" != "OFF" ] && VARIABLES+=("-DDEVICE_TYPE=$DEVICE_TYPE")

# Clean the build directory if the clean flag is set
[ "$CLEAN" == true ] && rm -rf $MAIN_DIR/build

printf "Compiling the example codes and Neko-TOP\n"
cmake -B $MAIN_DIR/build -S $MAIN_DIR "${VARIABLES[@]}"
cmake --build $MAIN_DIR/build --parallel

# ============================================================================ #
# Print the status of the build

printf "Neko-TOP Installation Complete\n"
printf "=%.0s" {1..80} && printf "\n"
printf "Neko installed to:\n"
printf "\t$NEKO_DIR\n"
printf "Supported features:\n"
printf "\tMPI: YES\n"
printf "\tTests: " && [[ "$TEST" == true ]] && printf "YES\n" || printf "NO\n"
printf "\tDevice: $DEVICE_TYPE\n"
printf "=%.0s" {1..80} && printf "\n"
if [ "$TEST" == true ]; then
    printf "To run the tests, execute the following command:\n"
    printf "\tctest --test-dir $MAIN_DIR/build\n"
fi
