#!/bin/bash
set -e # Exit with non-zero exit code if anything fails
# ============================================================================ #
# Print the help message
function help() {
    echo -e "Usage: $0 [options]"
    echo -e "Options:"
    echo -e "\t-h, --help        Show this help message and exit"
    echo -e "\t-t, --tests       Run the tests after the installation"
    echo -e "\t-c, --clean       Clean the build directory before compiling"
    echo -e "\t-q, --quiet       Suppress output"
    echo -e "\t-d, --device      Device type to compile for (off, CUDA, HIP)"
    echo -e "\t-e, --examples    Build the examples"
    echo -e "\t    --docs        Build the documentation"
    echo -e ""
    echo -e "Compilation and setup of Neko-TOP, this script will install all"
    echo -e "the dependencies and compile the Neko-TOP code."
    echo -e ""
    echo -e "Environment Variables:"
    echo -e "\tNEKO_DIR          The directory where Neko is installed"
    echo -e "\tJSON_FORTRAN_DIR  The directory where JSON-Fortran is installed"
    echo -e "\tADIOS2_DIR        The directory where ADIOS2 is installed"
    echo -e "\tNEK5000_DIR       The directory where Nek5000 is installed"
    echo -e "\tPFUNIT_DIR        The directory where PFUnit is installed"
    echo -e "\tGSLIB_DIR         The directory where GSLIB is installed"
    echo -e "\tCUDA_DIR          The directory where CUDA is installed"
    echo -e "\tCUDA_ARCH         CUDA architecture (required for --device CUDA, e.g. 80)"
    echo -e "\tHIP_DIR           The directory where HIP is installed"
    echo -e "\tBLAS_DIR          The directory where BLAS is installed"
    echo -e "\tCMAKE_VARIABLES   Additional variables to pass to CMake"
    echo -e "\tNEKO_CONFIG_FLAGS Additional features to pass to neko configure"
}

# ============================================================================ #
# Set main directories

export CURRENT_DIR=$(pwd)
export MAIN_DIR=$(dirname $(realpath $0))
export EXTERNAL_DIR="$MAIN_DIR/external"

# ============================================================================ #
# Parse the options

# Assign default values to the options
DEVICE_TYPE="CPU"
CLEAN=false
CLEAN_NEKO=false
QUIET=false
TEST=OFF
DOCS=OFF
EXAMPLES=OFF

# List possible options
OPTIONS=help,tests,clean,clean-neko,quiet,device:,docs,examples
OPT=h,t,c,q,d:,e

# Parse the inputs for options
PARSED=$(getopt --options=$OPT --longoptions=$OPTIONS --name "$0" -- "$@")
eval set -- "$PARSED"

# Loop through the options and set the variables
while true; do
    case "$1" in
    "-h" | "--help") help && exit ;;                  # Print help
    "-t" | "--tests") TEST="ON" && shift ;;            # Build the tests
    "-c" | "--clean") CLEAN=true && shift ;;          # Clean compilation
    "-q" | "--quiet") QUIET=true && shift ;;          # Suppress output
    "-d" | "--device") DEVICE_TYPE="$2" && shift 2 ;; # Device type
    "-e" | "--examples") EXAMPLES="ON" && shift ;;    # Build the examples

    # Purely long settings
    "--docs") DOCS="ON" && shift ;;             # Build the documentation
    "--clean-neko") CLEAN_NEKO=true && shift ;; # Clean Neko

    # End of options
    "--") shift && break ;;
    esac
done

[ "$CLEAN_NEKO" == true ] && CLEAN=true

export TEST CLEAN CLEAN_NEKO QUIET DEVICE_TYPE

# ============================================================================ #
# Execute the preparation script if it exists and prepare the environment

printf "=%.0s" {1..80} && printf "\n"
printf "Preparing environment.\n"

# Execute the preparation script if it exists
if [ -f "$MAIN_DIR/prepare.env" ]; then
    source $MAIN_DIR/prepare.env
fi
source $MAIN_DIR/scripts/dependencies.sh

# Define standard compilers if they are not defined as environment variables
if [ -z "$CC" ]; then export CC=$(which mpicc); else export CC; fi
if [ -z "$CXX" ]; then export CXX=$(which mpicxx); else export CXX; fi
if [ -z "$FC" ]; then export FC=$(which mpifort); else export FC; fi
if [ -z "$MPIFC" ]; then export MPIFC=$(which mpif90); else export MPIFC; fi
if [ -z "$MPICC" ]; then export MPICC=$(which mpicc); else export MPICC; fi
if [ -z "$MPICXX" ]; then export MPICXX=$(which mpicxx); else export MPICXX; fi

# Device specific compilers
if [ "$DEVICE_TYPE" == "CUDA" ]; then
    if [ -z "$NVCC" ]; then export NVCC=$(which nvcc); else export NVCC; fi
    if [ -z "$CUDA_ARCH" ]; then echo >&2 "CUDA_ARCH is not set."; exit 1; fi
elif [ "$DEVICE_TYPE" == "HIP" ]; then
    if [ -z "$HIPCC" ]; then export HIPCC=$(which hipcc); else export HIPCC; fi
fi

# Everything past this point should be general across all setups.
# ============================================================================ #
# Install dependencies (See scripts/dependencies.sh for details)

printf "=%.0s" {1..80} && printf "\n"
printf "Setting up external dependencies\n"

check_system_dependencies                      # Check for system dependencies.
find_json_fortran $JSON_FORTRAN_DIR            # Re-defines the JSON_FORTRAN_DIR variable.
find_adios2 $ADIOS2_DIR                        # Re-defines the ADIOS2_DIR variable.
find_neko $NEKO_DIR                            # Re-defines the NEKO_DIR variable.
[ "$TEST" == "ON" ] && find_pfunit $PFUNIT_DIR # Re-defines the PFUNIT_DIR variable.

# Done setting up external dependencies
# ============================================================================ #
# Compile the Neko-TOP and example codes.

printf "=%.0s" {1..80} && printf "\n"
printf "Compiling the example codes and Neko-TOP\n"

# Clean the build directory if the clean flag is set
[ "$CLEAN" == true ] && rm -fr $MAIN_DIR/build
mkdir -p $MAIN_DIR/build

# Validate and persist the POD Python runtime when example launchers are built.
if [ "$EXAMPLES" == "ON" ]; then
    find_pod_python_runtime $MAIN_DIR
fi

# If CMAKE_VARIABLES is a string, convert it to an array
if [ -n "$CMAKE_VARIABLES" ]; then
    CMAKE_VARIABLES=($CMAKE_VARIABLES)
else
    CMAKE_VARIABLES=()
fi

# Enable desired features
CMAKE_VARIABLES+=("-DBUILD_DOCS=$DOCS")
CMAKE_VARIABLES+=("-DBUILD_TESTING=$TEST")
CMAKE_VARIABLES+=("-DBUILD_EXAMPLES=$EXAMPLES")

# Run the cmake command to configure the build
cmake -B $MAIN_DIR/build -S $MAIN_DIR "${CMAKE_VARIABLES[@]}"

# Clean the build directory if the clean flag is set
cmake --build $MAIN_DIR/build --parallel

# ============================================================================ #
# Print the status of the build

printf "=%.0s" {1..80} && printf "\n"
printf "Neko-TOP Installation Complete\n"
printf "=%.0s" {1..80} && printf "\n"
printf "Neko installed to:\n"
printf "\t$NEKO_DIR\n"
printf "Supported features:\n"
printf "\tMPI:           YES\n"
printf "\tDevice:        $DEVICE_TYPE\n"
printf "\tTests:         " && [[ "$TEST" == "ON" ]] && printf "YES\n" || printf "NO\n"
printf "\tExamples:      " && [[ "$EXAMPLES" == "ON" ]] && printf "YES\n" || printf "NO\n"
printf "\tDocumentation: " && [[ "$DOCS" == "ON" ]] && printf "YES\n" || printf "NO\n"
printf "\tHDF5:          " && [[ -d "$HDF5_DIR" ]] && printf "YES\n" || printf "NO\n"
printf "=%.0s" {1..80} && printf "\n"
