#!/bin/bash

# ============================================================================ #
# Ensure JSON-Fortran is installed, if not install it.
function find_json_fortran() {

    # Ensure JSON-Fortran is installed, if not install it.
    if [ -z "$(find $1 -name libjsonfortran.so)" ]; then
        cmake -S $1 -B $1/build \
            --install-prefix $1 \
            -Wno-dev \
            -DUSE_GNU_INSTALL_CONVENTION=ON \
            -DSKIP_DOC_GEN=ON

        cmake --build $1/build --parallel
        cmake --install $1/build

    fi
    rm -fr $1/build

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
    if [ -z "$(find $1 -name libgs.a)" ]; then
        current=$(pwd)
        cd $1
        CC=mpicc ./install
        cd $current
    fi

    GSLIB=$(find $1 -type d -exec test -f '{}/lib/libgs.a' \; -print | head -n 1)
    if [ -z "$GSLIB" ]; then
        error "GSLIB not found at:"
        error "\t$1"
        error "Please set 1 to the directory containing the GSLIB"
        error "source code."
        error "You can download the source code from:"
        error "\thttps://github.com/bsmithyman/gslib"
        error "Or invoke the git submodule command:"
        error "\tgit submodule update --init --recursive"
        exit 1
    fi

    export GSLIB=$(realpath $GSLIB)
}

# ============================================================================ #
# Ensure pFunit is installed, if not install it.
function find_pfunit() {

    # Install pFunit
    PFUNIT=$(find $1 -type d -exec test -f '{}'/libpfunit.a \; -print)

    if [ -z "$(find $1 -name libpfunit.a)" ]; then

        cmake -S $1 -B $1/build -G "Unix Makefiles" \
            --install-prefix $1 \
            -DSKIP_MPI=False

        cmake --build $1/build --parallel
        cmake --install $1/build
        rm -fr $1/build
    fi

    PFUNIT=$(find $1 -type d -exec test -f '{}'/libpfunit.a \; -print)
    if [ -z "$PFUNIT" ]; then
        error "pFunit not found at:"
        error "\t$1"
        error "Please set PFUNIT_DIR to the directory containing the"
        error "pFunit source code."
        error "You can download the source code from:"
        error "\thttps://github.com/Goddard-Fortran-Ecosystem/pFUnit"
        error "Or invoke the git submodule command:"
        error "\tgit submodule update --init --recursive"
        exit 1
    fi

    export PFUNIT=$(realpath $PFUNIT)
}

# ============================================================================ #
# Helper function to print errors
function error() {
    echo -e "$1" >&2
}
