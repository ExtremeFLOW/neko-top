#!/bin/bash
set -e # Exit with nonzero exit code if anything fails
# ============================================================================ #
# Print the help message
function help() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -h, --help        Show this help message and exit"
    echo "  -t, --test        Run the tests after the installation"
    echo "  -c, --clean       Clean the build directory before compiling"
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
    -t | --test) TEST=1 ;;
    -c | --clean) CLEAN=1 ;;
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

# Check system dependencies
check_system_dependencies

# Ensure local dependencies are used if they are not defined as environment
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

find_json_fortran $JSON_FORTRAN_DIR # Defines the JSON_FORTRAN variable.
find_gslib $GSLIB_DIR               # Defines the GSLIB variable.

# Define optional features
[ ! -z "$GSLIB" ] && FEATURES+="--with-gslib=$GSLIB"
[ ! -z "$CUDA_DIR" ] && FEATURES+=" --with-cuda=$CUDA_DIR"
[ ! -z "$BLAS_DIR" ] && FEATURES+=" --with-blas=$BLAS_DIR"

# Done settng up external dependencies
# ============================================================================ #
# Install Neko

if [ ! -z "$CLEAN" ]; then
    printf "Cleaning the build directory\n"
    rm -rf $MAIN_DIR/build
    cd $NEKO_DIR
    git clean -dfx
    cd $CURRENT_DIR
fi

# Ensure Neko is installed, if not install it.
if [[ -z "$(find $NEKO_DIR -name libneko.a)" ]]; then
    printf "Neko is not installed, installing it now\n"
    if [ ! -f "$NEKO_DIR/regen.sh" ]; then
        git submodule update --init external/neko
    fi
    cd $NEKO_DIR
    ./regen.sh
    ./configure --prefix=$NEKO_DIR $FEATURES
    make -j install
    cd $CURRENT_DIR
fi
export PKG_CONFIG_PATH=$NEKO_DIR/lib/pkgconfig:$PKG_CONFIG_PATH

# ============================================================================ #
# Install PFUnit

if [[ "$TEST" && -z "$(find $PFUNIT_DIR -name libpfunit.a)" ]]; then

    if [ ! -f "$PFUNIT_DIR/CMakeLists.txt" ]; then
        printf "Installing PFUnit\n"
        git submodule update --init external/pFUnit
    fi

    cmake -B $PFUNIT_DIR/build -S $PFUNIT_DIR -DCMAKE_INSTALL_PREFIX=$PFUNIT_DIR
    cmake --build $PFUNIT_DIR/build --parallel
    cmake --install $PFUNIT_DIR/build
fi
export PFUNIT_DIR=$(find $PFUNIT_DIR -type d -exec test -f '{}'/lib/libpfunit.a \; -print)

# ============================================================================ #
# Compile the example codes.

VARIABLES="-DJSON_FORTRAN_DIR=$JSON_FORTRAN"
[ "$TEST" ] && VARIABLES+=" -DBUILD_TESTING=ON -DPFUNIT_DIR=$PFUNIT_DIR/cmake"

printf "Compiling the example codes and Neko-TOP\n"
cmake -B $MAIN_DIR/build -S $MAIN_DIR $VARIABLES
cmake --build $MAIN_DIR/build --parallel

# ============================================================================ #
# Print the status of the build

printf "Neko-TOP Installation Complete\n"
printf "===============================\n"
printf "Neko installed to: $NEKO_DIR\n"
printf "Supported features:\n"
printf "\tCUDA: " && [[ $FEATURES == *"cuda"* ]] && printf "YES\n" || printf "NO\n"
printf "\tMPI: YES\n"
printf "\tOpenCL: NO\n"
printf "\tTests: " && [[ $TEST ]] && printf "YES\n" || printf "NO\n"
printf "===============================\n"
printf "To run the tests, execute the following command:\n"
printf "\tctest --test-dir $MAIN_DIR/build\n"
