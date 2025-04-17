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
    check_external_dir

    # Determine the JSON-Fortran installation directory
    if [ ! -z "$1" ]; then
        JSON_FORTRAN_DIR="$(realpath $1)"
    elif [ -z "$JSON_FORTRAN_DIR" ]; then
        JSON_FORTRAN_DIR="$(realpath $EXTERNAL_DIR/json-fortran)"
    fi

    # Clone JSON-Fortran from the repository if it does not exist.
    if [[ ! -d $JSON_FORTRAN_DIR || $(ls -A $JSON_FORTRAN_DIR | wc -l) -eq 0 ]]; then
        [ -z "$JSON_FORTRAN_VERSION" ] && JSON_FORTRAN_VERSION="master"

        git clone --depth=1 --branch $JSON_FORTRAN_VERSION \
            https://github.com/jacobwilliams/json-fortran $JSON_FORTRAN_DIR
    fi

    # Ensure JSON-Fortran is installed, if not install it.
    if [[ -z "$(find $JSON_FORTRAN_DIR -name libjsonfortran.so)" ]]; then
        cmake -S $JSON_FORTRAN_DIR -B $JSON_FORTRAN_DIR/build \
            --install-prefix $JSON_FORTRAN_DIR \
            -Wno-dev \
            -DUSE_GNU_INSTALL_CONVENTION=ON \
            -DSKIP_DOC_GEN=ON

        cmake --build $JSON_FORTRAN_DIR/build --parallel
        cmake --install $JSON_FORTRAN_DIR/build
        rm -fr $JSON_FORTRAN_DIR/build
    fi

    # Add JSON-Fortran to the environment variables
    JSON_FORTRAN_LIB=$(find $JSON_FORTRAN_DIR -type d \
        -exec test -f '{}'/libjsonfortran.so \; -print)
    if [ -z "$JSON_FORTRAN_LIB" ]; then
        error "JSON-Fortran not found at:"
        error "\t$JSON_FORTRAN_DIR"
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
    check_external_dir

    # Determine the Nek5000 installation directory
    if [ ! -z "$1" ]; then
        NEK5000_DIR="$(realpath $1)"
    elif [ -z "$NEK5000_DIR" ]; then
        NEK5000_DIR="$(realpath $EXTERNAL_DIR/nek5000)"
    fi

    if [[ ! -d $NEK5000_DIR || $(ls -A $NEK5000_DIR | wc -l) -eq 0 ]]; then
        [ -z "$NEK5000_VERSION" ] && NEK5000_VERSION="master"

        git clone --depth 1 --branch $NEK5000_VERSION \
            https://github.com/Nek5000/Nek5000.git $NEK5000_DIR
    fi
}

# ============================================================================ #
# Ensure GSLIB is installed, if not install it.
function find_gslib() {
    check_external_dir

    # Determine the GSLib installation directory
    if [ ! -z "$1" ]; then
        GSLIB_DIR="$(realpath $1)"
    elif [ -z "$GSLIB_DIR" ]; then
        GSLIB_DIR="$(realpath $EXTERNAL_DIR/gslib)"
    fi

    # Clone GSLIB from the repository if it does not exist.
    if [ ! -d $GSLIB_DIR ]; then
        git clone --depth 1 --branch master \
            https://github.com/nek5000/gslib.git $GSLIB_DIR
    fi

    # Ensure GSLIB is installed, if not install it.
    if [ -z "$(find $GSLIB_DIR -name libgs.a)" ]; then
        echo "Building GSLIB"

        [ -z "$CURRENT_DIR" ] && CURRENT_DIR=$(pwd)
        cd $GSLIB_DIR

        make CC=mpicc
        make install DESTDIR=.
        rm -fr build

        cd $CURRENT_DIR
        echo "GSLIB built"
    fi

    # Add GSLIB to the environment variables
    GSLIB_LIB=$(find $GSLIB_DIR -type d -name 'lib*' \
        -exec test -f '{}/libgs.a' \; -print)
    if [ -z "$GSLIB_LIB" ]; then
        error "GSLIB not found at:"
        error "\t$GSLIB_DIR"
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
    check_external_dir

    # Determine the pFUnit installation directory
    if [ ! -z "$1" ]; then
        PFUNIT_DIR="$(realpath $1)"
    elif [ -z "$PFUNIT_DIR" ]; then
        PFUNIT_DIR="$(realpath $EXTERNAL_DIR/pfunit)"
    fi

    # Clone pFUnit from the repository if it does not exist.
    if [[ ! -d $PFUNIT_DIR || $(ls -A $PFUNIT_DIR | wc -l) -eq 0 ]]; then
        [ -z "$PFUNIT_VERSION" ] && PFUNIT_VERSION="v4.4.2"

        git clone --depth=1 --branch $PFUNIT_VERSION \
            https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git $PFUNIT_DIR

        # Patch pFUnit to work with Neko
        [ -z "$CURRENT_DIR" ] && CURRENT_DIR=$(pwd)
        cd $PFUNIT_DIR
        cat >>pfunit_error_stop.patch <<_ACEOF
diff --git a/src/funit/FUnit.F90 b/src/funit/FUnit.F90
index 7df7b65..4f7dbf5 100644
--- a/src/funit/FUnit.F90
+++ b/src/funit/FUnit.F90
@@ -168,6 +168,6 @@ contains
 #if defined(PGI)
          call exit(-1)
 #else
-         stop '*** Encountered 1 or more failures/errors during testing. ***'
+         error stop '*** Encountered 1 or more failures/errors during testing. ***'
 #endif
       end if
_ACEOF
        git apply pfunit_error_stop.patch
        cd $CURRENT_DIR
    fi

    if [[ -z "$(find $PFUNIT_DIR -name libpfunit.a)" ]]; then
        cmake -B $PFUNIT_DIR/build -S $PFUNIT_DIR -G "Unix Makefiles" \
            -DCMAKE_INSTALL_PREFIX=$PFUNIT_DIR
        cmake --build $PFUNIT_DIR/build --parallel
        cmake --install $PFUNIT_DIR/build
    fi

    PFUNIT_LIB=$(find $PFUNIT_DIR -type d -name 'lib*' \
        -exec test -f '{}'/libpfunit.a \; -print)
    if [ -z "$PFUNIT_LIB" ]; then
        error "pFUnit not found at:"
        error "\t$PFUNIT_DIR"
        error "Please set PFUNIT_DIR to the directory containing"
        error "the pFUnit source code."
        error "You can download the source code from:"
        error "\thttps://github.com/Goddard-Fortran-Ecosystem/pFUnit.git"
        exit 1
    fi

    export PFUNIT_DIR=$(realpath $PFUNIT_LIB/../)
}

# ============================================================================ #
# Ensure HDF5 is installed, if not install it.
function find_hdf5() {
    check_external_dir

    # Determine the HDF5 installation directory
    if [ ! -z "$1" ]; then
        HDF5_DIR="$(realpath $1)"
    elif [ -z "$HDF5_DIR" ]; then
        HDF5_DIR="$(realpath $EXTERNAL_DIR/hdf5)"
    fi

    # Clone HDF5 from the repository if it does not exist.
    if [ ! -d $HDF5_DIR ]; then
        [ -z "$HDF5_VERSION" ] && HDF5_VERSION="hdf5_1.14.6 "
        git clone --depth 1 --branch $HDF5_VERSION \
            https://github.com/HDFGroup/hdf5.git $HDF5_DIR
    fi

    # Ensure HDF5 is installed, if not install it.
    if [[ -z "$(find $HDF5_DIR -name libhdf5.so)" ]]; then
        cmake -B $HDF5_DIR/build -S $HDF5_DIR --install-prefix $HDF5_DIR \
            -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx \
            -DCMAKE_Fortran_COMPILER=mpifort -DHDF5_ENABLE_PARALLEL=ON \
            -DHDF5_BUILD_FORTRAN=ON -DHDF5_ENABLE_SZIP_SUPPORT:BOOL=OFF \
            -DCMAKE_BUILD_TYPE=Release
        cmake --build $HDF5_DIR/build/ --config Release --parallel
        cmake --install $HDF5_DIR/build/
        rm -fr $HDF5_DIR/build
    fi

    # Add HDF5 to the environment variables
    HDF5_LIB=$(find $HDF5_DIR -type d -name 'lib*' \
        -exec test -f '{}'/libhdf5_fortran.so \; -print)
    if [ -z "$HDF5_LIB" ]; then
        error "HDF5 not found at:"
        error "\t$HDF5_DIR"
        error "Please set HDF5_DIR to the directory containing"
        error "the HDF5 source code."
        error "You can download the source code from:"
        error "\thttps://github.com/HDFGroup/hdf5.git"
        exit 1
    fi

    export HDF5_DIR=$(realpath $HDF5_LIB/../)
    export LD_LIBRARY_PATH="$HDF5_LIB:$LD_LIBRARY_PATH"
    export PKG_CONFIG_PATH="$HDF5_LIB/pkgconfig:$PKG_CONFIG_PATH"
}

# ============================================================================ #
# Ensure Neko is installed, if not install it.
function find_neko() {
    check_external_dir

    # Find the required dependencies for Neko
    find_json_fortran $JSON_FORTRAN_DIR
    find_gslib
    find_hdf5
    [ "$TEST" == true ] && find_pfunit

    # Determine the Neko installation directory
    if [ ! -z "$1" ]; then
        NEKO_DIR="$(realpath $1)"
    elif [ -z "$NEKO_DIR" ]; then
        NEKO_DIR="$(realpath $EXTERNAL_DIR/neko)"
    fi

    # Clone Neko from the repository if it does not exist.
    if [[ ! -d $NEKO_DIR || $(ls -A $NEKO_DIR | wc -l) -eq 0 ]]; then
        [ -z "$NEKO_VERSION" ] && NEKO_VERSION="develop"

        git clone --depth 1 --branch $NEKO_VERSION \
            https://github.com/ExtremeFLOW/neko.git $NEKO_DIR
    fi

    # Check if Neko is installed, if not install it.
    if [[ -z "$(find $NEKO_DIR/lib*/ -name libneko.a)" || "$CLEAN_NEKO" == true ]]; then

        # Determine available features
        FEATURES="--enable-contrib "
        [ ! -z "$GSLIB_DIR" ] && FEATURES+="--with-gslib=$GSLIB_DIR"
        [ ! -z "$BLAS_DIR" ] && FEATURES+=" --with-blas=$BLAS_DIR"
        [ ! -z "$HDF5_DIR" ] && FEATURES+=" --with-hdf5=$HDF5_DIR"
        [ "$TEST" == true ] && FEATURES+=" --with-pfunit=$PFUNIT_DIR"

        # Handle device specific features
        if [ "$DEVICE_TYPE" == "CUDA" ]; then
            if [ -d "$CUDA_DIR" ]; then
                FEATURES+=" --with-cuda=$CUDA_DIR"
            else
                error "CUDA_DIR is not set."
                error "Please set CUDA_DIR to the directory containing"
                error "the CUDA installation."
                exit 1
            fi
        elif [ "$DEVICE_TYPE" == "HIP" ]; then
            if [ -d "$HIP_DIR" ]; then
                FEATURES+=" --with-hip=$HIP_DIR"
            else
                error "HIP_DIR is not set."
                error "Please set HIP_DIR to the directory containing"
                error "the HIP installation."
                exit 1
            fi
        elif [ "$DEVICE_TYPE" != "NONE" ]; then
            printf "Device type not recognized: $DEVICE_TYPE\n"
            printf "\tValid options are: CUDA, HIP or NONE\n"
            printf "\tPlease submit an issue if you would like to see additional options.\n"
            exit 1
        fi

        # Add additional features to be applied by the user
        FEATURES+=" ${EXTRA_FEATURES[@]}"

        [ -z "$CURRENT_DIR" ] && CURRENT_DIR=$(pwd)
        cd $NEKO_DIR
        if [[ ! -f "configure" || "$CLEAN_NEKO" == true ]]; then
            ./regen.sh
        fi
        if [[ ! -f Makefile || "$CLEAN_NEKO" == true ]]; then
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
        [ "$CLEAN_NEKO" == true ] && make clean
        [ "$QUIET" == true ] && make -s -j install || make -j install
        [ "$TEST" == true ] && make check

        # Verify installation device type
        if [ "$DEVICE_TYPE" == "CUDA" ]; then
            # Look for the line "  integer, parameter :: NEKO_BCKND_CUDA = 1"
            if [ -z "$(grep "NEKO_BCKND_CUDA = 1" src/config/neko_config.f90)" ]; then
                error "CUDA backend not found in Neko."
                error "Please ensure that the CUDA installation is correct."
                exit 1
            fi
        elif [ "$DEVICE_TYPE" == "HIP" ]; then
            # Look for the line "  integer, parameter :: NEKO_BCKND_CUDA = 1"
            if [ -z "$(grep "NEKO_BCKND_HIP = 1" src/config/neko_config.f90)" ]; then
                error "HIP backend not found in Neko."
                error "Please ensure that the HIP installation is correct."
                exit 1
            fi
        fi

        cd $CURRENT_DIR
    fi

    NEKO_DIR=$(find $NEKO_DIR -type d -exec test -f '{}'/lib/libneko.a \; -print)
    if [ -z "$NEKO_DIR" ]; then
        error "Neko not found at:"
        error "\$tNEKO_DIR"
        error "Please set NEKO_DIR to the directory containing"
        error "the Neko source code."
        error "You can download the source code from:"
        error "\thttps://github.com/ExtremeFLOW/neko.git"
        exit 1
    fi

    export NEKO_DIR=$(realpath $NEKO_DIR)
    export PKG_CONFIG_PATH=$NEKO_DIR/lib/pkgconfig:$PKG_CONFIG_PATH
    export LD_LIBRARY_PATH=$NEKO_DIR/lib:$LD_LIBRARY_PATH
    export PATH=$NEKO_DIR/bin:$PATH
}

# ============================================================================ #
# Enssure Cubit is installed, if not install it.
function find_cubit() {

    # Check if Cubit are available
    if [ ! -z "$(which cubit)" ]; then
        cubit=$(which cubit)
    elif [ ! -z "$(which coreform_cubit)" ]; then
        cubit=$(which coreform_cubit)
    elif [ ! -z "$(which trelis)" ]; then
        cubit=$(which trelis)
    else
        error "Cubit not found."
        error "Please ensure it is installed and available in the PATH."
        exit 1
    fi

    # Setup the alias and export the variable
    export cubit

}

# ============================================================================ #
# Ensure ExodusII to Nek5000 is installed, if not install it.
function find_exo2nek() {
    find_nek5000

    # Check if exo2nek is available
    if [ ! -z "$(which exo2nek)" ]; then
        exo2nek=$(which exo2nek)
    elif [ -f "$NEK5000_DIR/bin/exo2nek" ]; then
        export PATH=$NEK5000_DIR/bin:$PATH
    elif [ -f "$NEK5000_DIR/tools/maketools" ]; then

        [ -z "$CURRENT_DIR" ] && CURRENT_DIR=$(pwd)
        cd $NEK5000_DIR/tools

        ./maketools exo2nek
        cd $CURRENT_DIR
        export PATH=$NEK5000_DIR/bin:$PATH
    else
        error "exo2nek not found."
        error "Please ensure it is installed and available in the PATH."
        exit 1
    fi

}

# ============================================================================ #
# Ensure Rea2Nbin is installed.
function find_rea2nbin() {
    find_json_fortran $JSON_FORTRAN_DIR

    # Check if rea2nbin is available
    if [ ! -z "$(which rea2nbin)" ]; then
        rea2nbin=$(which rea2nbin)
    elif [ -f "$NEKO_DIR/bin/rea2nbin" ]; then
        rea2nbin="$NEKO_DIR/bin/rea2nbin"
    else
        error "rea2nbin not found."
        error "Please ensure it is installed and available in the PATH."
        exit 1
    fi
    export rea2nbin
}

# ============================================================================ #
# Ensure Gmsh is installed.
function find_gmsh() {
    # Check if gmsh is available
    if [ ! -z "$(which gmsh)" ]; then
        gmsh=$(which gmsh)
    else
        error "Gmsh not found."
        error "Please ensure it is installed and available in the PATH."
        exit 1
    fi
    export gmsh
}

# ============================================================================ #
# Ensure Gmsh2Nek is installed.
function find_gmsh2nek() {

    # Check if gmsh2nek is available
    if [ ! -z "$(which gmsh2nek)" ]; then
        gmsh2nek=$(which gmsh2nek)
    elif [ -f "$NEKO_DIR/bin/gmsh2nek" ]; then
        gmsh2nek="$NEKO_DIR/bin/gmsh2nek"
    elif [ -f "$NEKO_DIR/contrib/gmsh2nek/compile.sh" ]; then
        [ -z "$CURRENT_DIR" ] && CURRENT_DIR=$(pwd)
        cd $NEKO_DIR/contrib/gmsh2nek
        ./compile.sh
        cd $CURRENT_DIR
        gmsh2nek="$NEKO_DIR/contrib/gmsh2nek/gmsh2nek"
    else
        error "gmsh2nek not found."
        error "Please ensure it is installed and available in the PATH."
        exit 1
    fi

    export gmsh2nek
}

# ============================================================================ #
# Helper function

function error() {
    echo -e "$1" >&2
}

function check_external_dir() {
    if [ -z "$EXTERNAL_DIR" ]; then
        echo "Environment EXTERNAL_DIR is not set."
        echo "Default path will be used: ~/tmp/external"
        export EXTERNAL_DIR=$(realpath ~/tmp/external)
    fi

    mkdir -p $EXTERNAL_DIR

}
