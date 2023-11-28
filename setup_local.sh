#!/bin/bash
set -e # Exit with nonzero exit code if anything fails
# ============================================================================ #
# This section contain all the settings which are unique for this specific
# setup. This should allow us to have a relatively platform independent way of
# running the neko installation.

# Ensure the modules are available for the compilation
# module --silent load mpi/4.1.4-gcc-12.2.0-binutils-2.39 openblas/0.3.23 cuda/12.2

# BLAS_DIR="/usr/lib/x86_64-linux-gnu/openblas"
CUDA_DIR="/usr/local/cuda"

# Everything past this point should be general across all setups.
# ============================================================================ #
# Install dependencies

MAIN_DIR=$(dirname $(realpath $0))
EXTERNAL_DIR="$MAIN_DIR/external"

FEATURES=""

cmake -S $EXTERNAL_DIR/json-fortran \
    -B $EXTERNAL_DIR/json-fortran/build \
    --install-prefix $EXTERNAL_DIR/json-fortran \
    -Wno-dev \
    -DUSE_GNU_INSTALL_CONVENTION=ON \
    -DSKIP_DOC_GEN=ON

cmake --build $EXTERNAL_DIR/json-fortran/build --parallel
cmake --install $EXTERNAL_DIR/json-fortran/build

JSON_FORTRAN_LIB=""
if [ -d "$EXTERNAL_DIR/json-fortran/lib" ]; then
    JSON_FORTRAN_LIB="$EXTERNAL_DIR/json-fortran/lib"
elif [ -d "$EXTERNAL_DIR/json-fortran/lib64" ]; then
    JSON_FORTRAN_LIB="$EXTERNAL_DIR/json-fortran/lib64"
fi

# Setup GSLIB
if [ -z "$GSLIB_DIR" ]; then
    GSLIB_DIR="$EXTERNAL_DIR/Nek5000/3rd_party/gslib"
fi

if [ ! -f "$GSLIB_DIR/lib*/libgs.a" ]; then
    if [ -f "$GSLIB_DIR/install" ]; then
        cd $GSLIB_DIR
        ./install
        GSLIB_DIR="$GSLIB_DIR/gslib/build/"
        cd $MAIN_DIR
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

# ============================================================================ #
# Install Neko

# Setup environment variables
export PKG_CONFIG_PATH="$JSON_FORTRAN_LIB/pkgconfig:$PKG_CONFIG_PATH"
export LD_LIBRARY_PATH="$JSON_FORTRAN_LIB:$LD_LIBRARY_PATH"

# Setup Neko
cd $EXTERNAL_DIR/neko
if [ ! -f "Makefile" ]; then
    ./regen.sh
    ./configure --prefix=$EXTERNAL_DIR/neko $FEATURES
fi
make install -j
cd ../../

# ============================================================================ #
# Compile the example codes.
rm -fr $MAIN_DIR/build/
cmake -G Ninja -B $MAIN_DIR/build/ -S $MAIN_DIR
cmake --build $MAIN_DIR/build/ --parallel --clean-first

# ============================================================================ #
# Print the status of the build

printf "\n\n"
printf "Neko-TOP Installation Complete\n"
printf "===============================\n"
printf "Neko installed to: $EXTERNAL_DIR/neko\n"
printf "Supported features:\n"
printf "\tCUDA:"
if [ -z "$CUDA_DIR" ]; then
    printf " NO\n"
else
    printf " YES\n"
fi
printf "\tMPI: YES\n"
printf "\tOpenCL: NO\n"

printf "\n\n"

# EOF
