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
# Cleanup the external dependencies to ensure a clean build

git submodule foreach git clean -fdx
git submodule update --force --remote --init

# ============================================================================ #
# Install dependencies

MAIN_DIR="$PWD"
EXTERNAL_DIR="$PWD/External"

CC=gcc CXX=g++ FC=gfortran \
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

# ============================================================================ #
# Install Neko

# Setup environment variables
export PKG_CONFIG_PATH="$JSON_FORTRAN_LIB/pkgconfig:$PKG_CONFIG_PATH"
export LD_LIBRARY_PATH="$JSON_FORTRAN_LIB:$LD_LIBRARY_PATH"

# Setup Neko
cd $EXTERNAL_DIR/neko
./regen.sh
./configure --prefix=$EXTERNAL_DIR/neko --with-cuda=$CUDA_DIR
make install -j
cd ../../
