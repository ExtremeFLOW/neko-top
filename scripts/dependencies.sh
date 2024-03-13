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
    if [[ -z "$(find $1 -name libjsonfortran.so)" && -f $1/CMakeLists.txt ]]; then
        cmake -S $1 -B $1/build \
            --install-prefix $1 \
            -Wno-dev \
            -DUSE_GNU_INSTALL_CONVENTION=ON \
            -DSKIP_DOC_GEN=ON

        cmake --build $1/build --parallel
        cmake --install $1/build
        rm -fr $1/build
    fi

    JSON_FORTRAN=$(find $1 -type d -exec test -f '{}'/libjsonfortran.so \; -print)
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
    export LD_LIBRARY_PATH="$JSON_FORTRAN):$LD_LIBRARY_PATH"
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

    GSLIB="$(find $1 -type d -exec test -f '{}/lib/libgs.a' \; -print)"
    export GSLIB=$(realpath $GSLIB)
}

# ============================================================================ #
# Helper function to print errors
function error() {
    echo -e "$1" >&2
}
