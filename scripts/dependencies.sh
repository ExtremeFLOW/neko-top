#!/bin/bash

# ============================================================================ #
# Check system dependencies
function check_system_dependencies() {
    MISSING=()

    [ -z "$(which cmake)" ] && MISSING+=("CMake")
    [ -z "$(which make)" ] && MISSING+=("Make")
    [ -z "$(which git)" ] && MISSING+=("Git")
    [ -z "$(which mpicc)" ] && MISSING+=("MPICC")
    [ -z "$(which mpif90)" ] && MISSING+=("MPIF90")
    [ -z "$(which aclocal)" ] && MISSING+=("Aclocal")
    [ -z "$(which autoconf)" ] && MISSING+=("Autoconf")
    [ -z "$(which automake)" ] && MISSING+=("Automake")
    [ -z "$(which pkg-config)" ] && MISSING+=("Pkg-config")

    if [ ! -z "$MISSING" ]; then
        printf "The following dependencies are not installed:\n"
        for dep in "${MISSING[@]}"; do
            printf "  - $dep\n"
        done
        exit 1
    fi
}

# ============================================================================ #
# Ensure JSON-Fortran is installed, if not install it.
function find_json_fortran() {

    # Ensure JSON-Fortran is installed, if not install it.
    if [[ -z "$(find $1 -name libjsonfortran.so)" ]]; then
        if [ ! -f "$1/CMakeLists.txt" ]; then
            git clone https://github.com/jacobwilliams/json-fortran $1
        fi

        cmake -S $1 -B $1/build \
            --install-prefix $1 \
            -Wno-dev \
            -DUSE_GNU_INSTALL_CONVENTION=ON \
            -DSKIP_DOC_GEN=ON

        cmake --build $1/build --parallel
        cmake --install $1/build
        rm -fr $1/build
    fi

    JSON_FORTRAN=$(find $1 -type d \
        -exec test -f '{}'/libjsonfortran.so \; -print)
    if [ -z "$JSON_FORTRAN" ]; then
        error "JSON-Fortran not found at:"
        error "\t$1"
        error "Please set JSON_FORTRAN_DIR to the directory containing"
        error "the JSON-Fortran source code."
        error "You can download the source code from:"
        error "\thttps://github.com/jacobwilliams/json-fortran"
        error "Or invoke the git submodule command:"
        error "\tgit submodule update --init --recursive"
        exit 1
    fi

    export JSON_FORTRAN=$(realpath $JSON_FORTRAN)

    # Setup environment variables
    export PKG_CONFIG_PATH="$JSON_FORTRAN/pkgconfig:$PKG_CONFIG_PATH"
    export LD_LIBRARY_PATH="$JSON_FORTRAN:$LD_LIBRARY_PATH"
}

# ============================================================================ #
# Ensure GSLIB is installed, if not install it.
function find_gslib() {
    if [ ! -d $1 ]; then return; fi

    if [ -z "$(find $1 -name libgs.a)" ]; then
        current=$(pwd)
        cd $1
        CC=mpicc ./install
        cd $current
    fi

    GSLIB_DIR="$(find $1 -type d -exec test -f '{}/lib/libgs.a' \; -print)"
    if [ -z "$GSLIB_DIR" ]; then
        error "GSLIB not found at:"
        error "\t$1"
        error "Please set GSLIB_DIR to the directory containing"
        error "the GSLIB source code."
        error "You can download the source code from:"
        error "\thttps://github.com/Nek5000/gslib"
        error "Or invoke the git submodule command:"
        error "\tgit submodule update --init --recursive"
        exit 1
    fi

    export GSLIB_DIR=$(realpath $GSLIB_DIR)
}

# ============================================================================ #
# Ensure PFUnit is installed, if not install it.

function find_pfunit() {

    if [ ! $TEST ]; then return; fi

    if [[ -z "$(find $1 -name libpfunit.a)" ]]; then

        if [ ! -f "$1/CMakeLists.txt" ]; then
            git clone https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git $1
        fi

        cmake -B $1/build -S $1 -G "Unix Makefiles" \
            -DCMAKE_INSTALL_PREFIX=$1
        cmake --build $1/build --parallel
        cmake --install $1/build
    fi

    PFUNIT_DIR=$(find $1 -type d -exec test -f '{}'/lib/libpfunit.a \; -print)
    if [ -z "$PFUNIT_DIR" ]; then
        error "pFUnit not found at:"
        error "\t$1"
        error "Please set PFUNIT_DIR to the directory containing"
        error "the pFUnit source code."
        error "You can download the source code from:"
        error "\thttps://github.com/Goddard-Fortran-Ecosystem/pFUnit.git"
        error "Or invoke the git submodule command:"
        error "\tgit submodule update --init --recursive"
        exit 1
    fi

    export PFUNIT_DIR=$(realpath $PFUNIT_DIR)
}
# ============================================================================ #
# Helper function to print errors
function error() {
    echo -e "$1" >&2
}
