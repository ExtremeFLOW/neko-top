#!/bin/bash

# ============================================================================ #
# Check system dependencies
function check_system_dependencies() {
    MISSING=()

    if ! command -v cmake 2>&1 1>/dev/null; then MISSING+=("CMake"); fi
    if ! command -v make 2>&1 1>/dev/null; then MISSING+=("Make"); fi
    if ! command -v git 2>&1 1>/dev/null; then MISSING+=("Git"); fi
    if ! command -v aclocal 2>&1 1>/dev/null; then MISSING+=("Aclocal"); fi
    if ! command -v autoconf 2>&1 1>/dev/null; then MISSING+=("Autoconf"); fi
    if ! command -v automake 2>&1 1>/dev/null; then MISSING+=("Automake"); fi
    if ! command -v pkg-config 2>&1 1>/dev/null; then MISSING+=("Pkg-config"); fi
    [ -z "$MPICC" ] && MISSING+=("MPICC Environment Variable")
    [ -z "$MPICXX" ] && MISSING+=("MPICXX Environment Variable")
    [ -z "$MPIFC" ] && MISSING+=("MPIFC Environment Variable")

    if [ ! -z "$MISSING" ]; then
        printf "The following dependencies are not installed:\n"
        for dep in "${MISSING[@]}"; do
            printf "  - $dep\n"
        done
        exit 1
    fi

    # Check for correct version of cmake executable (>= 3.21)
    CMAKE_VERSION=$(cmake --version | grep -oP '(?<=version )\d+\.\d+')
    if [ $(echo "$CMAKE_VERSION < 3.21" | bc) -eq 1 ]; then
        printf "CMake version >= 3.21 is required.\n"
        printf "Please update your CMake installation.\n"
        exit 1
    fi

}

# ============================================================================ #
# Ensure JSON-Fortran is installed, if not install it.
function find_json_fortran() {
    check_external_dir

    # Determine the JSON-Fortran installation directory
    if [[ $# -ge 1 ]]; then
        JSON_FORTRAN_DIR="$1"
    elif [ -z "$JSON_FORTRAN_DIR" ]; then
        JSON_FORTRAN_DIR="json-fortran"
    fi

    if [ "${JSON_FORTRAN_DIR:0:1}" != "/" ]; then
        JSON_FORTRAN_DIR="$EXTERNAL_DIR/$JSON_FORTRAN_DIR"
    fi

    # Ensure JSON-Fortran is installed, if not install it.
    JSON_FORTRAN_LIB=$(find $JSON_FORTRAN_DIR -type d -name 'lib*' \
        -exec test -f '{}'/libjsonfortran.so \; -print 2>/dev/null) || true
    if [[ ! -d "$JSON_FORTRAN_LIB" ]]; then

        # Clone JSON-Fortran from the repository if it does not exist.
        if [[ ! -d "$JSON_FORTRAN_DIR" || $(ls -A $JSON_FORTRAN_DIR | wc -l) -eq 0 ]]; then
            [ -z "$JSON_FORTRAN_VERSION" ] && JSON_FORTRAN_VERSION="master"

            git clone --depth=1 --branch $JSON_FORTRAN_VERSION \
                https://github.com/jacobwilliams/json-fortran $JSON_FORTRAN_DIR
        fi

        # Install JSON-Fortran
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
    JSON_FORTRAN_LIB=$(find $JSON_FORTRAN_DIR -type d -name 'lib*' \
        -exec test -f '{}'/libjsonfortran.so \; -print 2>/dev/null) || true
    if [ ! -d "$JSON_FORTRAN_LIB" ]; then
        error "JSON-Fortran not found at:"
        error "\t$JSON_FORTRAN_DIR"
        error "Please set JSON_FORTRAN_DIR to the directory containing"
        error "the JSON-Fortran source code."
        error "You can download the source code from:"
        error "\thttps://github.com/jacobwilliams/json-fortran"
        exit 1
    fi

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
    if [[ $# -ge 1 ]]; then
        NEK5000_DIR="$1"
    elif [ -z "$NEK5000_DIR" ]; then
        return
    fi

    if [ "${NEK5000_DIR:0:1}" != "/" ]; then
        NEK5000_DIR="$EXTERNAL_DIR/$NEK5000_DIR"
    fi

    if [[ ! -d "$NEK5000_DIR" || $(ls -A $NEK5000_DIR | wc -l) -eq 0 ]]; then
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
    if [[ $# -ge 1 ]]; then
        GSLIB_DIR="$1"
    elif [ -z "$GSLIB_DIR" ]; then
        return
    fi

    if [ "${GSLIB_DIR:0:1}" != "/" ]; then
        GSLIB_DIR="$EXTERNAL_DIR/$GSLIB_DIR"
    fi

    # Ensure GSLIB is installed, if not install it.
    GSLIB_LIB=$(find $GSLIB_DIR -type d -name 'lib*' \
        -exec test -f '{}/libgs.a' \; -print 2>/dev/null) || true
    if [ ! -d "$GSLIB_LIB" ]; then

        # Clone GSLIB from the repository if it does not exist.
        if [ ! -d "$GSLIB_DIR" ]; then
            git clone --depth 1 --branch master \
                https://github.com/nek5000/gslib.git $GSLIB_DIR
        fi

        # Install GSLIB
        echo "Building GSLIB"

        [ -z "$CURRENT_DIR" ] && CURRENT_DIR=$(pwd)
        cd $GSLIB_DIR

        make CC=$MPICC
        make install DESTDIR=.
        rm -fr build

        cd $CURRENT_DIR
        echo "GSLIB built"
    fi

    # Add GSLIB to the environment variables
    GSLIB_LIB=$(find $GSLIB_DIR -type d -name 'lib*' \
        -exec test -f '{}/libgs.a' \; -print 2>/dev/null) || true
    if [ ! -d "$GSLIB_LIB" ]; then
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
    if [[ $# -ge 1 ]]; then
        PFUNIT_DIR="$1"
    elif [ -z "$PFUNIT_DIR" ]; then
        return
    fi

    if [ "${PFUNIT_DIR:0:1}" != "/" ]; then
        PFUNIT_DIR="$EXTERNAL_DIR/$PFUNIT_DIR"
    fi

    # Clone pFUnit from the repository if it does not exist.
    if [[ ! -d "$PFUNIT_DIR" || $(ls -A $PFUNIT_DIR | wc -l) -eq 0 ]]; then
        [ -z "$PFUNIT_VERSION" ] && PFUNIT_VERSION="v4.12.0"

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
            -DCMAKE_INSTALL_PREFIX=$PFUNIT_DIR \
            -DCMAKE_C_COMPILER=$CC \
            -DCMAKE_Fortran_COMPILER=$FC
        cmake --build $PFUNIT_DIR/build
        cmake --install $PFUNIT_DIR/build
    fi

    PFUNIT_LIB=$(find $PFUNIT_DIR -type d -name 'lib*' \
        -exec test -f '{}'/libpfunit.a \; -print 2>/dev/null) || true
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

    # Determine the HDF5 installation directory. HDF5_ROOT is the name CMake
    # and the module systems use; HDF5_DIR is accepted as a legacy spelling.
    if [[ $# -ge 1 ]]; then
        HDF5_ROOT="$1"
    elif [ -n "$HDF5_ROOT" ]; then
        : # already set in the environment
    elif [ -n "$HDF5_DIR" ]; then
        HDF5_ROOT="$HDF5_DIR"
    else
        return
    fi

    if [ "${HDF5_ROOT:0:1}" != "/" ]; then
        HDF5_ROOT="$EXTERNAL_DIR/$HDF5_ROOT"
    fi

    # Ensure HDF5 is installed, if not install it.
    HDF5_LIB=$(find "$HDF5_ROOT" -type d -name 'lib*' \
        -exec test -f '{}'/libhdf5_fortran.so \; -print 2>/dev/null) || true
    if [[ ! -d "$HDF5_LIB" ]]; then

        # Never try to build into a read-only prefix, such as one provided by
        # a module system.
        if [[ -d "$HDF5_ROOT" && ! -w "$HDF5_ROOT" ]]; then
            error "HDF5 not found under the read-only prefix:"
            error "\t$HDF5_ROOT"
            error "It looks module-provided. Load a module that supplies the"
            error "Fortran bindings, or set HDF5_ROOT to a writable path to"
            error "have one built there."
            exit 1
        fi

        # Clone HDF5 from the repository if it does not exist.
        if [[ ! -d "$HDF5_ROOT" || $(ls -A $HDF5_ROOT | wc -l) -eq 0 ]]; then
            [ -z "$HDF5_VERSION" ] && HDF5_VERSION="hdf5_2.0.0"
            git clone --depth 1 --branch $HDF5_VERSION \
                https://github.com/HDFGroup/hdf5.git $HDF5_ROOT
        fi

        # Build and install HDF5
        cmake -B $HDF5_ROOT/build -S $HDF5_ROOT \
            --install-prefix $HDF5_ROOT -DCMAKE_BUILD_TYPE=Release \
            -DCMAKE_C_COMPILER=$MPICC -DCMAKE_CXX_COMPILER=$MPICXX \
            -DCMAKE_Fortran_COMPILER=$MPIFC -DHDF5_ENABLE_PARALLEL=ON \
            -DHDF5_BUILD_FORTRAN=ON -DHDF5_ENABLE_SZIP_SUPPORT:BOOL=OFF \
            -DHDF5_BUILD_TOOLS:BOOL=ON
        cmake --build $HDF5_ROOT/build/ --config Release --parallel
        cmake --install $HDF5_ROOT/build/ --config Release
        rm -fr $HDF5_ROOT/build
    fi

    # Add HDF5 to the environment variables
    HDF5_LIB=$(find "$HDF5_ROOT" -type d -name 'lib*' \
        -exec test -f '{}'/libhdf5_fortran.so \; -print 2>/dev/null) || true
    if [[ ! -d "$HDF5_LIB" ]]; then
        error "HDF5 not found at:"
        error "\t$HDF5_ROOT"
        error "Please set HDF5_ROOT to the directory containing"
        error "the HDF5 installation."
        error "You can download the source code from:"
        error "\thttps://github.com/HDFGroup/hdf5.git"
        exit 1
    fi

    export HDF5_ROOT=$(realpath $HDF5_LIB/../)
    export HDF5_DIR=$HDF5_ROOT   # legacy spelling, consumed by find_neko
    export LD_LIBRARY_PATH="$HDF5_LIB:$LD_LIBRARY_PATH"
    export PKG_CONFIG_PATH="$HDF5_LIB/pkgconfig:$PKG_CONFIG_PATH"
}

# ============================================================================ #
# Ensure ParMETIS is installed, if not install it.

function find_parmetis() {

    # Determine the Parmetis installation directory
    check_external_dir
    if [[ $# -ge 1 ]]; then
        PARMETIS_DIR="$1"
    elif [ -z "$PARMETIS_DIR" ]; then
        PARMETIS_DIR="parmetis"
    fi

    if [[ "${PARMETIS_DIR:0:1}" != "/" && "${PARMETIS_DIR:0:1}" != "~" ]]; then
        PARMETIS_DIR="$(realpath $EXTERNAL_DIR/$PARMETIS_DIR)"
    fi

    if [[ -z "$(find $PARMETIS_DIR -name libparmetis.a)" ]]; then
        [ -z "$CURRENT_DIR" ] && CURRENT_DIR=$(pwd)
        CMAKE_GENERATOR_OLD=$CMAKE_GENERATOR
        CMAKE_GENERATOR="Unix Makefiles"

        # Download and install ParMETIS
        mkdir -p $PARMETIS_DIR && cd $PARMETIS_DIR
        wget https://github.com/mfem/tpls/raw/refs/heads/gh-pages/parmetis-4.0.3.tar.gz
        tar xzf parmetis-4.0.3.tar.gz
        cd parmetis-4.0.3

        # Modify the minimum requirement of cmake
        cmake_lists=$(find . -name CMakeLists.txt)
        for file in $cmake_lists; do
            sed -i 's/cmake_minimum_required(VERSION 2.8)/cmake_minimum_required(VERSION 3.11)/g' $file
        done

        # Compile the bundled metis library
        cd metis
        make config prefix=${PARMETIS_DIR}
        make -j && make install
        cd ../

        # Compile parmetis
        make config prefix=${PARMETIS_DIR}
        make -j && make install
        cd ../
        rm -rf parmetis-4.0.3 parmetis-4.0.3.tar.gz
        CMAKE_GENERATOR=$CMAKE_GENERATOR_OLD
        cd $CURRENT_DIR
    fi

    PARMETIS_LIB=$(find $PARMETIS_DIR -type d -name 'lib*' \
        -exec test -f '{}'/libparmetis.a \; -print)
    if [ -z "$PARMETIS_LIB" ]; then
        error "ParMETIS not found at:"
        error "\t$PARMETIS_DIR"
        error "Please set PARMETIS_DIR to the directory containing"
        error "the ParMETIS source code."
        error "You can download the source code from:"
        error "\thttps://github.com/KarypisLab/ParMETIS.git"
        exit 1
    fi

    export PARMETIS_DIR=$(realpath $PARMETIS_DIR)
}

# ============================================================================ #
# Ensure Neko is installed, if not install it.
function find_neko() {
    check_external_dir

    # Find the required dependencies for Neko
    find_json_fortran $JSON_FORTRAN_DIR
    find_gslib $GSLIB_DIR
    find_hdf5 $HDF5_ROOT
    find_parmetis $PARMETIS_DIR
    [ -n "$PFUNIT_DIR" ] && find_pfunit $PFUNIT_DIR

    # Determine the Neko installation directory
    if [[ $# -ge 1 ]]; then
        NEKO_DIR="$1"
    elif [ -z "$NEKO_DIR" ]; then
        NEKO_DIR="neko"
    fi

    if [ "${NEKO_DIR:0:1}" != "/" ]; then
        NEKO_DIR="$EXTERNAL_DIR/$NEKO_DIR"
    fi

    # Check if Neko is installed, if not install it.
    NEKO_LIB=$(find $NEKO_DIR -type d -name 'lib*' -maxdepth 1 \
        -exec test -f '{}'/libneko.a \; -print 2>/dev/null) || true
    if [[ ! -d "$NEKO_LIB" || "$CLEAN_NEKO" == true ]]; then

        # Clone Neko from the repository if it does not exist.
        if [[ ! -d "$NEKO_DIR" || $(ls -A $NEKO_DIR | wc -l) -eq 0 ]]; then
            [ -z "$NEKO_VERSION" ] && NEKO_VERSION="neko-top"

            git clone --depth 1 --branch $NEKO_VERSION \
                https://github.com/ExtremeFLOW/neko.git $NEKO_DIR

        fi

        # Apply Cray-specific patches before building on Cray systems
        if [[ -n "${CRAYPE_VERSION:-}" || "${PE_ENV:-}" == "CRAY" || -d "/opt/cray" ]]; then
            cray_patches=(
                "patches/cce_stack.patch"
                "patches/cce_time_state.patch"
                "patches/cce_openmp.patch"
            )
            for patch in "${cray_patches[@]}"; do
                if git -C "$NEKO_DIR" apply --check "$patch" 2>/dev/null; then
                    git -C "$NEKO_DIR" apply "$patch"
                fi
            done
        fi

        # Determine available features
        FEATURES="--enable-contrib"
        [ -n "$GSLIB_DIR" ] && FEATURES+=" --with-gslib=$GSLIB_DIR"
        [ -n "$BLAS_DIR" ] && FEATURES+=" --with-blas=$BLAS_DIR"
        [ -n "$HDF5_DIR" ] && FEATURES+=" --with-hdf5=$HDF5_DIR"
        [ -n "$PARMETIS_DIR" ] && FEATURES+=" --with-parmetis=$PARMETIS_DIR"
        [ -n "$PFUNIT_DIR" ] && FEATURES+=" --with-pfunit=$PFUNIT_DIR"

        # Handle device specific features
        if [ "$DEVICE_TYPE" == "CUDA" ]; then
            if [ -n "$CUDA_DIR" ]; then
                FEATURES+=" --with-cuda=$CUDA_DIR"
            else
                error "CUDA_DIR is not set."
                error "Please set CUDA_DIR to the directory containing"
                error "the CUDA installation."
                exit 1
            fi

            if [ -n "$NEKO_CUDA_ARCH" ]; then
                FEATURES+=" CUDA_ARCH=$NEKO_CUDA_ARCH"
            elif [ -n "$CUDA_ARCH" ]; then
                FEATURES+=" CUDA_ARCH=-arch=sm_$CUDA_ARCH"
            else
                error "CUDA architecture not set."
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
        elif [ "$DEVICE_TYPE" != "CPU" ]; then
            printf "Device type not recognized: $DEVICE_TYPE\n"
            printf "\tValid options are: CUDA, HIP or NONE\n"
            printf "\tPlease submit an issue if you would like to see"
            printf " additional options.\n"
            exit 1
        fi

        # Add additional features to be applied by the user
        FEATURES+=" ${NEKO_CONFIG_FLAGS[@]}"

        [ -z "$CURRENT_DIR" ] && CURRENT_DIR=$(pwd)
        cd $NEKO_DIR
        if [[ ! -f "configure" || "$CLEAN_NEKO" == true ]]; then
            ./regen.sh
        fi
        if [[ ! -f Makefile || "$CLEAN_NEKO" == true ]]; then
            ./configure --prefix="$(realpath ./)" $FEATURES \
                FC=$FC MPIFC=$MPIFC FCFLAGS="$NEKO_FCFLAGS" \
                CC=$CC MPICC=$MPICC MPICXX=$MPICXX CFLAGS="$NEKO_CFLAGS" \
                HIPCC=$HIPCC HIP_HIPCC_FLAGS="$NEKO_HIPCC_FLAGS" \
                CUDA_CFLAGS="$NEKO_CUDA_CFLAGS"
        fi

        # Update compile dependencies if makedepf90 is installed
        if command -v makedepf90 2>&1 1>/dev/null; then
            size_pre=$(stat -c %s src/.depends)
            cd src/ && make depend && cd ../
            if [ "$size_pre" != "$(stat -c %s src/.depends)" ]; then
                automake -a
                rm -fr autom4te.cache
            fi
        fi

        [ "$CLEAN_NEKO" == true ] && make clean
        [ "$QUIET" == true ] && make -s -j || make -j
        [ "$NEKO_TEST" == true ] && make check -j
        make install

        cd $CURRENT_DIR

        # Revert the patches to keep the repository clean
        if [[ -n "${CRAYPE_VERSION:-}" || "${PE_ENV:-}" == "CRAY" || -d "/opt/cray" ]]; then
            for patch in "${cray_patches[@]}"; do
                if git -C "$NEKO_DIR" apply --reverse --check "$patch" 2>/dev/null; then
                    git -C "$NEKO_DIR" apply --reverse "$patch"
                fi
            done
        fi
    fi

    NEKO_LIB=$(find $NEKO_DIR -type d -name 'lib*' -maxdepth 1 \
        -exec test -f '{}'/libneko.a \; -print 2>/dev/null) || true
    if [ ! -d "$NEKO_LIB" ]; then
        error "Neko not found at:"
        error "\t$NEKO_DIR"
        error "Please set NEKO_DIR to the directory containing"
        error "the Neko source code."
        error "You can download the source code from:"
        error "\thttps://github.com/ExtremeFLOW/neko.git"
        exit 1
    fi

    # Check the device type supported by neko
    if [[ -n "$DEVICE_TYPE" && -f "$NEKO_LIB/pkgconfig/neko.pc" ]]; then
        PATTERN="(?<=backend=).*"
        NEKO_DEVICE=$(grep -oP "$PATTERN" $NEKO_LIB/pkgconfig/neko.pc) || true

        if [[ -n "$NEKO_DEVICE" && "$DEVICE_TYPE" != "$NEKO_DEVICE" ]]; then
            error "Neko device type does not match the requested device type."
            error "Please ensure that the Neko installation is correct."
            error "Requested device type: $DEVICE_TYPE"
            error "Neko device type: $NEKO_DEVICE"
            exit 1
        fi
    fi

    export NEKO_DIR=$(realpath $NEKO_LIB/../)
    export PKG_CONFIG_PATH=$NEKO_LIB/pkgconfig:$PKG_CONFIG_PATH
    export LD_LIBRARY_PATH=$NEKO_LIB:$LD_LIBRARY_PATH
    export PATH=$NEKO_DIR/bin:$PATH
}

# ============================================================================ #
# Enssure Cubit is installed, if not install it.
function find_cubit() {

    # Check if Cubit are available
    if command -v cubit 2>&1 1>/dev/null; then
        export cubit=$(command -v cubit)
    elif command -v coreform_cubit 2>&1 1>/dev/null; then
        export cubit=$(command -v coreform_cubit)
    elif command -v trelis 2>&1 1>/dev/null; then
        export cubit=$(command -v trelis)
    else
        error "Cubit not found."
        error "Please ensure it is installed and available in the PATH."
        exit 1
    fi
}

# ============================================================================ #
# Ensure ExodusII to Nek5000 is installed, if not install it.
function find_exo2nek() {
    find_nek5000 $NEK5000_DIR

    # Check if exo2nek is available
    if command -v exo2nek 2>&1 1>/dev/null; then
        export exo2nek=$(command -v exo2nek)
    elif [ -x "$NEK5000_DIR/bin/exo2nek" ]; then
        export exo2nek="$NEK5000_DIR/bin/exo2nek"
    elif [ -x "$NEK5000_DIR/tools/maketools" ]; then

        [ -z "$CURRENT_DIR" ] && CURRENT_DIR=$(pwd)
        cd $NEK5000_DIR/tools

        ./maketools exo2nek
        cd $CURRENT_DIR
        export exo2nek="$NEK5000_DIR/bin/exo2nek"
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
    find_hdf5 $HDF5_ROOT

    # Check if rea2nbin is available
    if command -v rea2nbin 2>&1 1>/dev/null; then
        export rea2nbin=$(command -v rea2nbin)
    elif [ -x "$NEKO_DIR/bin/rea2nbin" ]; then
        export rea2nbin="$NEKO_DIR/bin/rea2nbin"
    else
        error "rea2nbin not found."
        error "Please ensure it is installed and available in the PATH."
        exit 1
    fi
}

# ============================================================================ #
# Ensure Gmsh is installed.
function find_gmsh() {
    # Check if gmsh is available
    if ! command -v gmsh 2>&1 1>/dev/null; then
        error "Gmsh not found."
        error "Please ensure it is installed and available in the PATH."
        exit 1
    fi
}

# ============================================================================ #
# Ensure Gmsh2Nek is installed.
function find_gmsh2nek() {

    # Check if gmsh2nek is available
    if command -v gmsh2nek 2>&1 1>/dev/null; then
        export gmsh2nek=$(command -v gmsh2nek)
    elif [ -x "$NEKO_DIR/bin/gmsh2nek" ]; then
        export gmsh2nek="$NEKO_DIR/bin/gmsh2nek"
    elif [ -f "$NEKO_DIR/contrib/gmsh2nek/compile.sh" ]; then
        tmp_dir=$(pwd)
        cd $NEKO_DIR/contrib/gmsh2nek
        ./compile.sh
        cd $tmp_dir
        export gmsh2nek="$NEKO_DIR/contrib/gmsh2nek/gmsh2nek"
    else
        error "gmsh2nek not found."
        error "Please ensure it is installed and available in the PATH."
        exit 1
    fi
}

# ============================================================================ #
# Helper function

function error() {
    echo -e "$1" >&2
}

function check_external_dir() {
    if [ -z "$EXTERNAL_DIR" ]; then
        echo "Environment EXTERNAL_DIR is not set."
        echo "Default path will be used: $HOME/tmp/external"
        export EXTERNAL_DIR=$HOME/tmp/external
    fi

    mkdir -p $EXTERNAL_DIR

}
