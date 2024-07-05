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

    JSON_FORTRAN_LIB=$(find $1 -type d \
        -exec test -f '{}'/libjsonfortran.so \; -print)
    if [ -z "$JSON_FORTRAN_LIB" ]; then
        error "JSON-Fortran not found at:"
        error "\t$1"
        error "Please set JSON_FORTRAN_DIR to the directory containing"
        error "the JSON-Fortran source code."
        error "You can download the source code from:"
        error "\thttps://github.com/jacobwilliams/json-fortran"
        exit 1
    fi

    JSON_FORTRAN_LIB=$(realpath $JSON_FORTRAN_LIB)

    # Setup environment variables
    export JSON_FORTRAN_DIR=$(realpath $JSON_FORTRAN_LIB/../)
    export PKG_CONFIG_PATH="$JSON_FORTRAN_LIB/pkgconfig:$PKG_CONFIG_PATH"
    export LD_LIBRARY_PATH="$JSON_FORTRAN_LIB:$LD_LIBRARY_PATH"
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

    GSLIB_LIB=$(find $1 -type d -name 'lib*' \
        -exec test -f '{}/libgs.a' \; -print)
    if [ -z "$GSLIB_LIB" ]; then
        error "GSLIB not found at:"
        error "\t$1"
        error "Please set GSLIB_DIR to the directory containing"
        error "the GSLIB source code."
        error "You can download the source code from:"
        error "\thttps://github.com/Nek5000/gslib"
        exit 1
    fi

    export GSLIB_DIR=$(realpath $GSLIB_LIB/../)
}

# ============================================================================ #
# Ensure PFUnit is installed, if not install it.
function find_pfunit() {

    if [ ! $TEST ]; then return; fi

    if [[ ! -d $1 || $(ls -A $1 | wc -l) -eq 0 ]]; then
        [ -z "$PFUNIT_VERSION" ] && PFUNIT_VERSION="v4.4.2"

        git clone --depth=1 --branch $PFUNIT_VERSION \
            https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git $1
    fi

    if [[ -z "$(find $1 -name libpfunit.a)" ]]; then
        cmake -B $1/build -S $1 -G "Unix Makefiles" \
            --install-prefix=$(realpath $1)
        cmake --build $1/build --parallel
        cmake --install $1/build
    fi

    PFUNIT_LIB=$(find $1 -type d -name 'lib*' \
        -exec test -f '{}'/libpfunit.a \; -print)
    if [ -z "$PFUNIT_LIB" ]; then
        error "pFUnit not found at:"
        error "\t$1"
        error "Please set PFUNIT_DIR to the directory containing"
        error "the pFUnit source code."
        error "You can download the source code from:"
        error "\thttps://github.com/Goddard-Fortran-Ecosystem/pFUnit.git"
        exit 1
    fi

    export PFUNIT_DIR=$(realpath $PFUNIT_LIB/../)
}

# ============================================================================ #
# Ensure GKlib is installed, if not install it. (Dependency of PARMetis)

function find_gklib() {
    if [[ ! -d $1 || $(ls -A $1 | wc -l) -eq 0 ]]; then
        [ -z "$GKLIB_VERSION" ] && GKLIB_VERSION="master"

        git clone --depth=1 --branch $GKLIB_VERSION \
            https://github.com/KarypisLab/GKlib.git $1
    fi

    if [[ -z "$(find $1 -name libGKlib.a)" ]]; then
        make -C $1 config prefix=$(realpath $1) CFLAGS="-D_POSIX_C_SOURCE=199309L"
        make -C $1 install -j
        # rm -fr $1/build
    fi

    GKLIB_LIB=$(find $1 -type d -name 'lib*' \
        -exec test -f '{}'/libGKlib.a \; -print)
    if [ -z "$GKLIB_LIB" ]; then
        error "GKlib not found at:"
        error "\t$1"
        error "Please set GKLIB_DIR to the directory containing"
        error "the GKlib source code."
        error "You can download the source code from:"
        error "\thttps://github.com/KarypisLab/GKlib"
        exit 1
    fi

    #     export GKLIB_DIR=$(realpath $GKLIB_LIB/../)
}

# ============================================================================ #
# Ensure METIS is installed, if not install it. (Dependency of PARMetis)

function find_metis() {
    if [[ ! -d $1 || $(ls -A $1 | wc -l) -eq 0 ]]; then
        [ -z "$METIS_VERSION" ] && METIS_VERSION="v5.2.1"

        git clone --depth=1 --branch $METIS_VERSION \
            https://github.com/KarypisLab/METIS.git $1
    fi

    if [[ -z "$(find $1 -name libmetis.a)" ]]; then
        cd $1
        make config prefix=$(realpath $1) gklib_path=$GKLIB_DIR
        make install
        # rm -fr $1/build
        cd $CURRENT_DIR
    fi

    METIS_LIB=$(find $1 -type d -name 'lib*' \
        -exec test -f '{}'/libmetis.a \; -print)
    if [ -z "$METIS_LIB" ]; then
        error "METIS not found at:"
        error "\t$1"
        error "Please set METIS_DIR to the directory containing"
        error "the METIS source code."
        error "You can download the source code from:"
        error "\thttps://github.com/KarypisLab/METIS.git"
        exit 1
    fi

    # export METIS_DIR=$(realpath $METIS_LIB/../)
}

# ============================================================================ #
# Ensure ParMETIS is installed, if not install it.

function find_parmetis() {
    # if [[ ! -d $1 || $(ls -A $1 | wc -l) -eq 0 ]]; then
    #     [ -z "$PARMETIS_VERSION" ] && PARMETIS_VERSION="main"

    #     git clone --depth=1 --branch $PARMETIS_VERSION \
    #         https://github.com/KarypisLab/ParMETIS.git $1
    # fi

    # [ -z "$GKLIB_DIR" ] && GKLIB_DIR=$(realpath $1/GKlib)
    # [ -z "$METIS_DIR" ] && METIS_DIR=$(realpath $1/METIS)

    # export IDXTYPEWIDTH=32
    # export REALTYPEWIDTH=32

    # find_gklib $GKLIB_DIR
    # find_metis $METIS_DIR

    # if [[ -z "$(find $1 -name libparmetis.a)" ]]; then
    #     cd $1
    #     make config cc=mpicc cxx=mpicxx prefix=$(realpath $1) \
    #         gklib_path=$GKLIB_DIR metis_path=$METIS_DIR
    #     make install 1>out.log 2>error.log
    #     cd $CURRENT_DIR
    # fi

    mkdir -p $1
    cd $1
    wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz
    tar xzf parmetis-4.0.3.tar.gz && cd parmetis-4.0.3 && make config prefix=$1
    make -j$(nproc) && make install
    cd $CURRENT_DIR

    PARMETIS_LIB=$(find $1 -type d -name 'lib*' \
        -exec test -f '{}'/libparmetis.a \; -print)
    if [ -z "$PARMETIS_LIB" ]; then
        error "ParMETIS not found at:"
        error "\t$1"
        error "Please set PARMETIS_DIR to the directory containing"
        error "the ParMETIS source code."
        error "You can download the source code from:"
        error "\thttps://github.com/KarypisLab/ParMETIS.git"
        exit 1
    fi

    export PARMETIS_DIR=$(realpath $PARMETIS_LIB/../)
}

# ============================================================================ #
# Ensure Neko is installed, if not install it.
function find_neko() {

    # Clone Neko from the repository if it does not exist.
    if [[ ! -d $1 || $(ls -A $1 | wc -l) -eq 0 ]]; then
        [ -z "$NEKO_VERSION" ] && NEKO_VERSION="master"

        git clone --depth 1 --branch $NEKO_VERSION \
            https://github.com/ExtremeFLOW/neko.git $1
    fi

    # Determine available features
    FEATURES="--enable-contrib"
    [ ! -z "$GSLIB_DIR" ] && FEATURES+=" --with-gslib=$GSLIB_DIR"
    [ ! -z "$BLAS_DIR" ] && FEATURES+=" --with-blas=$BLAS_DIR"
    [ ! -z "$PARMETIS_DIR" ] && FEATURES+=" --with-parmetis=$PARMETIS_DIR"
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
