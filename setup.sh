#!/bin/bash
set -e # Exit with nonzero exit code if anything fails
# ============================================================================ #
# Print the help message
function help() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -h, --help        Show this help message and exit"
    echo ""
    echo "Compilation and setup of Neko-TOP, this script will install all the"
    echo "dependencies and compile the Neko-TOP code."
    echo ""
    echo "Environment Variables:"
    echo "  NEKO_DIR          The directory where Neko is installed"
    echo "  JSON_FORTRAN_DIR  The directory where JSON-Fortran is installed"
    echo "  NEK5000_DIR       The directory where Nek5000 is installed"
    echo "  PFUNIT_DIR        The directory where PFUnit is installed"
    echo "  GSLIB_DIR         The directory where GSLIB is installed"
    echo "  CUDA_DIR          The directory where CUDA is installed"
    echo "  BLAS_DIR          The directory where BLAS is installed"
    echo "  CC                The C compiler to use"
    echo "  CXX               The C++ compiler to use"
    echo "  FC                The Fortran compiler to use"
    echo "  NVCC              The CUDA compiler to use"
    exit 0
}

while [ "$1" != "" ]; do
    case $1 in
    -h | --help) help ;;
    esac
    shift
done

# ============================================================================ #
# Set main directories
CURRENT_DIR=$(pwd)
MAIN_DIR=$(dirname $(realpath $0))
EXTERNAL_DIR="$MAIN_DIR/external"

# ============================================================================ #
# Execute the preparation script if it exists and prepare the environment

# Execute the preparation script if it exists
if [ -f "$MAIN_DIR/prepare.sh" ]; then
    source $MAIN_DIR/prepare.sh
fi
source $MAIN_DIR/scripts/dependencies.sh

# Ensure local dependencies are used if they are not defined as environment
[ -z "$NEKO_DIR" ] && NEKO_DIR="$EXTERNAL_DIR/neko"
[ -z "$JSON_FORTRAN_DIR" ] && JSON_FORTRAN_DIR="$EXTERNAL_DIR/json-fortran"
[ -z "$NEK5000_DIR" ] && NEK5000_DIR="$EXTERNAL_DIR/Nek5000"
[ -z "$PFUNIT_DIR" ] && PFUNIT_DIR="$EXTERNAL_DIR/pfunit"
[ -z "$GSLIB_DIR" ] && GSLIB_DIR="$NEK5000_DIR/3rd_party/gslib"

# Define standard compilers if they are not defined as environment variables
if [ -z "$CC" ]; then export CC=$(which gcc); else export CC; fi
if [ -z "$CXX" ]; then export CXX=$(which g++); else export CXX; fi
if [ -z "$FC" ]; then export FC=$(which gfortran); else export FC; fi
if [ -z "$NVCC" ]; then export NVCC=$(which nvcc); else export NVCC; fi

# Everything past this point should be general across all setups.
# ============================================================================ #
# Install dependencies (See scripts/dependencies.sh for details)

find_json_fortran $JSON_FORTRAN_DIR # Defines the JSON_FORTRAN variable.
find_gslib $GSLIB_DIR               # Defines the GSLIB variable.
find_pfunit $PFUNIT_DIR             # Defines the PFUNIT variable.

# Define Neko features
FEATURES="--with-gslib=$GSLIB --with-pfunit=$PFUNIT"

# Define optional features
[ ! -z "$CUDA_DIR" ] && FEATURES+=" --with-cuda=$CUDA_DIR"
[ ! -z "$BLAS_DIR" ] && FEATURES+=" --with-blas=$BLAS_DIR"

# Done settng up external dependencies
# ============================================================================ #
# Install Neko

# Setup Neko
cd $NEKO_DIR
if ! $(make --quiet install -j); then
    ./regen.sh
    ./configure --prefix=$NEKO_DIR $FEATURES

    if ! make --quiet -j install; then
        printf "Neko installation failed\n" >&2
        exit 1
    fi
fi
cd $CURRENT_DIR

# ============================================================================ #
# Compile the example codes.

printf "Compiling the example codes and Neko-TOP\n"
cmake -B $MAIN_DIR/build/ -S $MAIN_DIR
cmake --build $MAIN_DIR/build/ --parallel

# ============================================================================ #
# Print the status of the build

printf "Neko-TOP Installation Complete\n"
printf "===============================\n"
printf "Neko installed to: $NEKO_DIR\n"
printf "Supported features:\n"
printf "\tCUDA: " && [[ $FEATURES == *"cuda"* ]] && printf "YES\n" || printf "NO\n"
printf "\tMPI: YES\n"
printf "\tOpenCL: NO\n"
