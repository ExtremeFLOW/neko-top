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

    if [[ ! -d $1 || $(ls -A $1 | wc -l) -eq 0 ]]; then
        [ -z "$JSON_FORTRAN_VERSION" ] && JSON_FORTRAN_VERSION="master"

        git clone --depth=1 --branch $JSON_FORTRAN_VERSION \
            https://github.com/jacobwilliams/json-fortran $1
    fi

    # Ensure JSON-Fortran is installed, if not install it.
    if [[ -z "$(find $1 -name libjsonfortran.so)" ]]; then
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
        exit 1
    fi

    export JSON_FORTRAN=$(realpath $JSON_FORTRAN)

    # Setup environment variables
    export PKG_CONFIG_PATH="$JSON_FORTRAN/pkgconfig:$PKG_CONFIG_PATH"
    export LD_LIBRARY_PATH="$JSON_FORTRAN:$LD_LIBRARY_PATH"
}

# ============================================================================ #
# Ensure Nek5000 is installed, if not install it.
function find_nek5000() {

    if [[ ! -d $1 || $(ls -A $1 | wc -l) -eq 0 ]]; then
        [ -z "$NEK5000_VERSION" ] && NEK5000_VERSION="master"

        git clone --depth 1 --branch $NEK5000_VERSION \
            https://github.com/Nek5000/Nek5000.git $1
    fi
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
        exit 1
    fi

    export GSLIB_DIR=$(realpath $GSLIB_DIR)
}

# ============================================================================ #
# Ensure PFUnit is installed, if not install it.
function find_pfunit() {

    if [ ! $TEST ]; then return; fi

    if [[ ! -d $1 || $(ls -A $1 | wc -l) -eq 0 ]]; then
        [ -z "$PFUNIT_VERSION" ] && PFUNIT_VERSION="main"

        git clone --depth=1 --branch $PFUNIT_VERSION \
            https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git $1
    fi

    if [[ -z "$(find $1 -name libpfunit.a)" ]]; then
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
        exit 1
    fi

    export PFUNIT_DIR=$(realpath $PFUNIT_DIR)
}

# ============================================================================ #
# Ensure Neko is installed, if not install it.
find_neko() {

    # Clone Neko from the repository if it does not exist.
    if [[ ! -d $1 || $(ls -A $1 | wc -l) -eq 0 ]]; then
        [ -z "$NEKO_VERSION" ] && NEKO_VERSION="master"

        git clone --depth 1 --branch $NEKO_VERSION \
            https://github.com/ExtremeFLOW/neko.git $1
    fi

    # Determine available features
    [ ! -z "$GSLIB_DIR" ] && FEATURES+="--with-gslib=$GSLIB_DIR"
    [ ! -z "$BLAS_DIR" ] && FEATURES+=" --with-blas=$BLAS_DIR"
    [ "$TEST" == true ] && FEATURES+=" --with-pfunit=$PFUNIT_DIR"

    # Handle device specific features
    if [ "$DEVICE_TYPE" == "CUDA" ]; then
        [ ! -z "$CUDA_DIR" ] && FEATURES+=" --with-cuda=$CUDA_DIR"
    elif [ "$DEVICE_TYPE" != "OFF" ]; then
        printf "Invalid device type: $DEVICE_TYPE\n"
        exit 1
    fi

    if [[ -z "$(find $1 -name libneko.a)" || "$CLEAN" == true ]]; then
        cd $1
        if [[ ! -f "configure" || "$CLEAN" == true ]]; then
            ./regen.sh
        fi
        if [[ ! -f Makefile || "$CLEAN" == true ]]; then
            ./configure --prefix="$(realpath ./)" $FEATURES
        fi

        # Update compile dependencies if makedepf90 is installed
        if [ ! -z "$(which makedepf90)" ]; then
            size_pre=$(stat -c %s src/.depends)
            cd src/ && make depend && cd ../
            if [ "$size_pre" != "$(stat -c %s src/.depends)" ]; then
                automake -a
                rm -fr autom4te.cache
            fi
        fi

        [ "$CLEAN" == true ] && make clean
        [ "$QUIET" == true ] && make -s -j install || make -j install

        # Run Tests if the flag is set
        if [ "$TEST" == true ]; then
            printf "Running Neko tests\n"
            make check
        fi
        cd $CURRENT_DIR
    fi

    NEKO_DIR=$(find $1 -type d -exec test -f '{}'/lib/libneko.a \; -print)
    if [ -z "$NEKO_DIR" ]; then
        error "Neko not found at:"
        error "\t$1"
        error "Please set NEKO_DIR to the directory containing"
        error "the Neko source code."
        error "You can download the source code from:"
        error "\thttps://github.com/ExtremeFLOW/neko.git"
        exit 1
    fi

    export NEKO_DIR=$(realpath $NEKO_DIR)
    export PKG_CONFIG_PATH=$NEKO_DIR/lib/pkgconfig:$PKG_CONFIG_PATH
}

# ============================================================================ #
# Helper function to print errors
function error() {
    echo -e "$1" >&2
}
