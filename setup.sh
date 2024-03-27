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
    echo -e "\tCC                The C compiler to use"
    echo -e "\tCXX               The C++ compiler to use"
    echo -e "\tFC                The Fortran compiler to use"
    echo -e "\tNVCC              The CUDA compiler to use"
    exit 0
}

for in in $@; do

    # Opions and flags
    if [[ ${in:0:2} == "--" ]]; then
        case "${in:2}" in
        "help") help ;;        # Print help
        "test") TEST=true ;;   # Build the tests
        "clean") CLEAN=true ;; # Clean compilation

        *)
            printf '  %-10s %-67s\n' "Invalid option:" "$in"
            exit 1
            ;;
        esac

    elif [[ ${in:0:1} == "-" ]]; then
        for ((i = 1; i < ${#in}; i++)); do
            case "${in:$i:1}" in
            "h") help ;;       # Print help
            "t") TEST=true ;;  # Build the tests
            "c") CLEAN=true ;; # Clean compilation

            *)
                printf '  %-10s %-67s\n' "Invalid option:" "${in:$i:1}"
                exit 1
                ;;
            esac
        done
    fi
done
export TEST=$TEST

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
find_pfunit $PFUNIT_DIR             # Defines the PFUNIT variable.

# Define optional features
[ ! -z "$GSLIB" ] && FEATURES+="--with-gslib=$GSLIB"
[ ! -z "$CUDA_DIR" ] && FEATURES+=" --with-cuda=$CUDA_DIR"
[ ! -z "$BLAS_DIR" ] && FEATURES+=" --with-blas=$BLAS_DIR"
[ $TEST ] && FEATURES+=" --with-pfunit=$PFUNIT_DIR"

# Done settng up external dependencies
# ============================================================================ #
# Install Neko

# Ensure Neko is installed, if not install it.
if [[ -z "$(find $NEKO_DIR -name libneko.a)" || $CLEAN ]]; then
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

# Run Tests if the flag is set
if [ $TEST ]; then
    printf "Running Neko tests\n"
    cd $NEKO_DIR
    make check
    cd $CURRENT_DIR
fi
export PKG_CONFIG_PATH=$NEKO_DIR/lib/pkgconfig:$PKG_CONFIG_PATH

# ============================================================================ #
# Compile the example codes.

# Clean the build directory if the clean flag is set
[ $CLEAN ] && rm -rf $MAIN_DIR/build

VARIABLES="-DJSON_FORTRAN_DIR=$JSON_FORTRAN"
[ $TEST ] && VARIABLES+=" -DBUILD_TESTING=ON -DPFUNIT_DIR=$PFUNIT_DIR/cmake"

printf "Compiling the example codes and Neko-TOP\n"
cmake -B $MAIN_DIR/build -S $MAIN_DIR $VARIABLES
cmake --build $MAIN_DIR/build --parallel

# ============================================================================ #
# Print the status of the build

printf "Neko-TOP Installation Complete\n"
printf "=%.0s" {1..80} && printf "\n"
printf "Neko installed to:\n"
printf "\t$NEKO_DIR\n"
printf "Supported features:\n"
printf "\tCUDA: " && [[ $FEATURES == *"cuda"* ]] && printf "YES\n" || printf "NO\n"
printf "\tMPI: YES\n"
printf "\tOpenCL: NO\n"
printf "\tTests: " && [[ $TEST ]] && printf "YES\n" || printf "NO\n"
printf "=%.0s" {1..80} && printf "\n"
if [ $TEST ]; then
    printf "To run the tests, execute the following command:\n"
    printf "\tctest --test-dir $MAIN_DIR/build\n"
fi
