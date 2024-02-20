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

# Execute the preparation script if it exists
if [ -f "$MAIN_DIR/prepare.sh" ]; then
    source $MAIN_DIR/prepare.sh
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
    PFUNIT_DIR="$EXTERNAL_DIR/pfunit"
fi
if [[ -z "$GSLIB_DIR" && -d "$NEK5000_DIR/3rd_party/gslib" ]]; then
    GSLIB_DIR="$NEK5000_DIR/3rd_party/gslib"
fi

# Define standard compilers if they are not defined as environment variables
if [ -z "$CC" ]; then export CC=$(which gcc); else export CC; fi
if [ -z "$CXX" ]; then export CXX=$(which g++); else export CXX; fi
if [ -z "$FC" ]; then export FC=$(which gfortran); else export FC; fi
if [ -z "$NVCC" ]; then export NVCC=$(which nvcc); else export NVCC; fi

# Everything past this point should be general across all setups.
# ============================================================================ #
# Install dependencies

# Install Json-Fortran
if [[ ! -f "$JSON_FORTRAN_DIR/lib*/jsonfortran.so" &&
    -f "$JSON_FORTRAN_DIR/CMakeLists.txt" ]]; then

    cmake -S $JSON_FORTRAN_DIR \
        -B $JSON_FORTRAN_DIR/build \
        --install-prefix $JSON_FORTRAN_DIR \
        -Wno-dev \
        -DUSE_GNU_INSTALL_CONVENTION=ON \
        -DSKIP_DOC_GEN=ON

    cmake --build $JSON_FORTRAN_DIR/build --parallel
    cmake --install $JSON_FORTRAN_DIR/build
fi

if [ -d "$JSON_FORTRAN_DIR/lib" ]; then
    JSON_FORTRAN_LIB="$JSON_FORTRAN_DIR/lib"
elif [ -d "$JSON_FORTRAN_DIR/lib64" ]; then
    JSON_FORTRAN_LIB="$JSON_FORTRAN_DIR/lib64"
else
    printf "Json-Fortran not found at JSON_FORTRAN_DIR: \n\t$JSON_FORTRAN_DIR" >&2
    exit 1
fi

# Setup GSLIB
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

# Install pFunit
if [ ! -f "$PFUNIT_DIR/PFUNIT-*/lib*/libpfunit.a" ]; then
    if [ ! -f "$PFUNIT_DIR/CMakeLists.txt" ]; then
        git submodule init $PFUNIT_DIR
    fi

    cmake -S $PFUNIT_DIR -B $PFUNIT_DIR/build -G "Unix Makefiles" \
        --install-prefix $PFUNIT_DIR \
        -DSKIP_MPI=False

    cmake --build $PFUNIT_DIR/build --parallel
    cmake --install $PFUNIT_DIR/build
fi

# Done settng up external dependencies
# ============================================================================ #
# Define features available to neko
FEATURES=""

if [ -f "$GSLIB_DIR/lib*/libgs.a" ]; then
    FEATURES+="--with-gslib=$GSLIB_DIR "
fi

if [ -d "$CUDA_DIR" ]; then
    FEATURES+="--with-cuda=$CUDA_DIR "
else
    CUDA_DIR=""
fi

if [ -d "$BLAS_DIR" ]; then
    FEATURES+="--with-blas=$BLAS_DIR "
fi

if [ -d "$(realpath $PFUNIT_DIR/PFUNIT-*)" ]; then
    PFUNIT_DIR="$(realpath $PFUNIT_DIR/PFUNIT-*)"
    FEATURES+="--with-pfunit=$PFUNIT_DIR "
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
cmake --build $MAIN_DIR/build/ --clean-first

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
