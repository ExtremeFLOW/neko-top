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
    echo -e "\t-d, --device      Device type to compile for (off, CUDA, HIP)"
    echo -e "\t    --doc         Build the documentation"
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
    echo -e "\tHIP_DIR           The directory where HIP is installed"
    echo -e "\tBLAS_DIR          The directory where BLAS is installed"
    echo -e "\tCMAKE_VARIABLES   Additional variables to pass to CMake"
}

# ============================================================================ #
# Set main directories

export CURRENT_DIR=$(pwd)
export MAIN_DIR=$(dirname $(realpath $0))
export EXTERNAL_DIR="$MAIN_DIR/external"

# ============================================================================ #
# Parse the options

# Assign default values to the options
DEVICE_TYPE="NONE"
CLEAN=false
CLEAN_NEKO=false
QUIET=false
TEST=false
DOCS=false

# List possible options
OPTIONS=help,test,clean,clean-neko,quiet,device:,doc
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

    # Purely long settings
    "--doc") DOCS=true && shift ;;              # Build the documentation
    "--clean-neko") CLEAN_NEKO=true && shift ;; # Clean Neko

    # End of options
    "--") shift && break ;;
    esac
done

# Check if the device type has changed
if [ -f "$MAIN_DIR/build/CMakeCache.txt" ]; then
    CURRENT_DEVICE_TYPE="$(grep -oP '(?<=DEVICE_TYPE:STRING=).*' $MAIN_DIR/build/CMakeCache.txt)"
    if [ "$DEVICE_TYPE" != "$CURRENT_DEVICE_TYPE" ]; then
        echo "Device type has changed, cleaning the build directory"
        CLEAN=true
        CLEAN_NEKO=true
    fi
fi

export TEST CLEAN CLEAN_NEKO QUIET DEVICE_TYPE

# ============================================================================ #
# Execute the preparation script if it exists and prepare the environment

printf "=%.0s" {1..80} && printf "\n"
printf "Preparing environment.\n\n"

# Execute the preparation script if it exists
if [ -f "$MAIN_DIR/prepare.env" ]; then
    source $MAIN_DIR/prepare.env
fi
source $MAIN_DIR/scripts/dependencies.sh

# Define standard compilers if they are not defined as environment variables
if [ -z "$CC" ]; then export CC=$(which mpicc); else export CC; fi
if [ -z "$CXX" ]; then export CXX=$(which mpicxx); else export CXX; fi
if [ -z "$FC" ]; then export FC=$(which mpifort); else export FC; fi

# Device specific compilers
if [ "$DEVICE_TYPE" == "CUDA" ]; then
    if [ -z "$NVCC" ]; then export NVCC=$(which nvcc); else export NVCC; fi
fi
if [ "$DEVICE_TYPE" == "HIP" ]; then
    if [ -z "$HIPCC" ]; then export HIPCC=$(which hipcc); else export HIPCC; fi
fi

# Everything past this point should be general across all setups.
# ============================================================================ #
# Install dependencies (See scripts/dependencies.sh for details)

printf "=%.0s" {1..80} && printf "\n"
printf "Setting up external dependencies\n"

check_system_dependencies          # Check for system dependencies.
find_json_fortran                  # Re-defines the JSON_FORTRAN_DIR variable.
find_nek5000                       # Re-defines the NEK5000_DIR variable.
find_neko                          # Re-defines the NEKO_DIR variable.
[ "$TEST" == true ] && find_pfunit # Re-defines the PFUNIT_DIR variable.

# Done settng up external dependencies
# ============================================================================ #
# Compile the Neko-TOP and example codes.

printf "=%.0s" {1..80} && printf "\n"
printf "Compiling the example codes and Neko-TOP\n"

# Set CMAKE_VARIABLES to pass to the cmake command
if [ -z "$CMAKE_VARIABLES" ]; then CMAKE_VARIABLES=(); fi

# If CMAKE_VARIABLES is a string, convert it to an array
if [ -n "$CMAKE_VARIABLES" ] && [ ! -z "$CMAKE_VARIABLES" ]; then
    CMAKE_VARIABLES=($CMAKE_VARIABLES)
fi

# Set the variables for the compilation
[ "$TEST" == true ] && CMAKE_VARIABLES+=("-DBUILD_TESTING=ON")
[ "$TEST" == true ] && CMAKE_VARIABLES+=("-DPFUNIT_DIR=$PFUNIT_DIR/cmake")
[ "$DEVICE_TYPE" != "OFF" ] && CMAKE_VARIABLES+=("-DDEVICE_TYPE=$DEVICE_TYPE")

# Set the documentation flag
if [ "$DOCS" == true ]; then
    CMAKE_VARIABLES+=("-DBUILD_DOCS=ON")
else
    CMAKE_VARIABLES+=("-DBUILD_DOCS=OFF")
fi

cmake -B $MAIN_DIR/build -S $MAIN_DIR "${CMAKE_VARIABLES[@]}"

# Clean the build directory if the clean flag is set
[ "$CLEAN" == true ] && cmake --build $MAIN_DIR/build --target clean
cmake --build $MAIN_DIR/build --parallel
cmake --build $MAIN_DIR/build --target Examples --parallel

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
    ctest -C Debug --output-on-failure --test-dir $MAIN_DIR/build --parallel
fi
