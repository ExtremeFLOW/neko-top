#!/bin/bash
set -e # Exit with nonzero exit code if anything fails
# ============================================================================ #
# This section contain all the settings which are unique for this specific
# setup. This should allow us to have a relatively platform independent way of
# running the neko installation.

# Set main directories
CURRENT_DIR=$(pwd)
MAIN_DIR=$(dirname $(realpath $0))
EXTERNAL_DIR="$MAIN_DIR/external"
FEATURES=""

# Load modules when module command is available
if [ "$(which module)" ]; then
    module --silent load mpi/4.1.4-gcc-12.2.0-binutils-2.39
    module --silent load openblas/0.3.23 cuda/12.2

    BLAS_DIR="/appl/OpenBLAS/0.3.23/$CPUTYPEV/gcc-12.2.0/lib/"
    CUDA_DIR="/appl/cuda/12.2.0"
elif [ "$(which spack)" ]; then
    spack env activate neko-top
fi

# Look for CUDA if it is not defined as environment variable
if [[ -z "$CUDA_DIR" && -d "/usr/local/cuda" ]]; then
    CUDA_DIR="/usr/local/cuda"
    NVCC="$CUDA_DIR/bin/nvcc"
fi

# Ensure local dependencies are used if they are not defined as environment
if [ -z "$NEKO_DIR" ]; then
    NEKO_DIR="$EXTERNAL_DIR/neko"
fi
if [ -z "$JSON_FORTRAN_DIR" ]; then
    JSON_FORTRAN_DIR="$EXTERNAL_DIR/json-fortran"
fi
if [ -z "$NEK5000_DIR" ]; then
    NEK5000_DIR="$EXTERNAL_DIR/Nek5000"
fi

if [ -z "$PFUNIT_DIR" ]; then
    PFUNIT_DIR="$EXTERNAL_DIR/pFUnit"
fi

# Define standard compilers if they are not defined as environment variables
if [ -z "$CC" ]; then export CC=$(which gcc); else export CC; fi
if [ -z "$CXX" ]; then export CXX=$(which g++); else export CXX; fi
if [ -z "$FC" ]; then export FC=$(which gfortran); else export FC; fi
if [ -z "$NVCC" ]; then export NVCC=$(which nvcc); else export NVCC; fi

# Everything past this point should be general across all setups.
# ============================================================================ #
# Install dependencies

cmake -S $JSON_FORTRAN_DIR \
    -B $JSON_FORTRAN_DIR/build \
    --install-prefix $JSON_FORTRAN_DIR \
    -Wno-dev \
    -DUSE_GNU_INSTALL_CONVENTION=ON \
    -DSKIP_DOC_GEN=ON

cmake --build $JSON_FORTRAN_DIR/build --parallel
cmake --install $JSON_FORTRAN_DIR/build

JSON_FORTRAN_LIB=""
if [ -d "$JSON_FORTRAN_DIR/lib" ]; then
    JSON_FORTRAN_LIB="$JSON_FORTRAN_DIR/lib"
elif [ -d "$JSON_FORTRAN_DIR/lib64" ]; then
    JSON_FORTRAN_LIB="$JSON_FORTRAN_DIR/lib64"
fi

# Done settng up external dependencies
# ============================================================================ #
# Define features available to neko

# Setup GSLIB
if [[ -z "$GSLIB_DIR" && -d "$NEK5000_DIR/3rd_party/gslib" ]]; then
    GSLIB_DIR="$NEK5000_DIR/3rd_party/gslib"
fi

if [ ! -f "$GSLIB_DIR/lib*/libgs.a" ]; then
    if [ -f "$GSLIB_DIR/install" ]; then
        cd $GSLIB_DIR
        CC=mpicc ./install
        GSLIB_DIR="$GSLIB_DIR/gslib/build/"
        cd $CURRENT_DIR
    else
        printf "GSLIB not found at GSLIB_DIR: \n\t$GSLIB_DIR" >&2
        exit 1
    fi
    FEATURES+="--with-gslib=$GSLIB_DIR "
fi

if [ -d "$CUDA_DIR" ]; then
    FEATURES+="--with-cuda=$CUDA_DIR "
else
    CUDA_DIR=""
fi

if [ ! -z "$BLAS_DIR" ]; then
    FEATURES+="--with-blas=$BLAS_DIR "
fi

# ============================================================================ #
# Install Neko

# Setup environment variables
export PKG_CONFIG_PATH="$JSON_FORTRAN_LIB/pkgconfig:$PKG_CONFIG_PATH"
export LD_LIBRARY_PATH="$JSON_FORTRAN_LIB:$LD_LIBRARY_PATH"

# Setup Neko
cd $NEKO_DIR
if ! make --quiet install -j; then
    ./regen.sh
    ./configure --prefix=$NEKO_DIR $FEATURES
    if ! make --quiet install; then
        printf "Neko installation failed\n" >&2
        exit 1
    fi
fi
cd $CURRENT_DIR

# ============================================================================ #
# Compile the example codes.
cmake -B $MAIN_DIR/build/ -S $MAIN_DIR
cmake --build $MAIN_DIR/build/ --parallel --clean-first

# ============================================================================ #
# Print the status of the build

printf "\n\n"
printf "Neko-TOP Installation Complete\n"
printf "===============================\n"
printf "Neko installed to: $NEKO_DIR\n"
printf "Supported features:\n"
printf "\tCUDA:"
if [ -z "$CUDA_DIR" ]; then printf " NO\n"; else printf " YES\n"; fi
printf "\tMPI: YES\n"
printf "\tOpenCL: NO\n"

printf "\n\n"

# EOF
