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

# Print the resolved executable path for command-substitution callers.
function find_python_executable() {
    local candidate
    local python_path

    if [ -n "${PYTHON_BIN:-}" ]; then
        if [ -x "${PYTHON_BIN}" ]; then
            if [[ "${PYTHON_BIN}" = /* ]]; then
                printf '%s\n' "${PYTHON_BIN}"
            elif [[ "${PYTHON_BIN}" == */* ]]; then
                python_path=$(cd "$(dirname "${PYTHON_BIN}")" && pwd)
                printf '%s/%s\n' "${python_path}" "$(basename "${PYTHON_BIN}")"
            else
                command -v "${PYTHON_BIN}"
            fi
            return 0
        elif command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
            command -v "${PYTHON_BIN}"
            return 0
        fi

        error "PYTHON_BIN is set but not executable:"
        error "\t${PYTHON_BIN}"
        return 1
    fi

    for candidate in python3 python; do
        if command -v "${candidate}" >/dev/null 2>&1; then
            command -v "${candidate}"
            return 0
        fi
    done

    return 1
}

# ============================================================================ #
# Ensure JSON-Fortran is installed, if not install it.
function find_json_fortran() {
    check_external_dir

    # Determine the JSON-Fortran installation directory
    if [[ $# -ge 1 ]]; then
        JSON_FORTRAN_DIR="$(realpath $1)"
    elif [ -z "$JSON_FORTRAN_DIR" ]; then
        JSON_FORTRAN_DIR="$(realpath $EXTERNAL_DIR/json-fortran)"
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
        if [[ "${1:0:1}" != "/" && "${1:0:1}" != "~" ]]; then
            NEK5000_DIR="$(realpath $EXTERNAL_DIR/$1)"
        else
            NEK5000_DIR="$(realpath $1)"
        fi
    else
        export NEK5000_DIR=""
        return
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
        if [[ "${1:0:1}" != "/" && "${1:0:1}" != "~" ]]; then
            GSLIB_DIR="$(realpath $EXTERNAL_DIR/$1)"
        else
            GSLIB_DIR="$(realpath $1)"
        fi
    else
        export GSLIB_DIR=""
        return
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

    if [[ "${PFUNIT_DIR:0:1}" != "/" && "${PFUNIT_DIR:0:1}" != "~" ]]; then
        PFUNIT_DIR="$(realpath $EXTERNAL_DIR/$PFUNIT_DIR)"
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

    # Determine the HDF5 installation directory
    check_external_dir
    if [[ $# -ge 1 ]]; then
        HDF5_DIR="$1"
    elif [ -z "$HDF5_DIR" ]; then
        return
    fi

    if [[ "${HDF5_DIR:0:1}" != "/" && "${HDF5_DIR:0:1}" != "~" ]]; then
        HDF5_DIR="$(realpath $EXTERNAL_DIR/$HDF5_DIR)"
    fi

    # Ensure HDF5 is installed, if not install it.
    HDF5_LIB=$(find $HDF5_DIR -type d -name 'lib*' \
        -exec test -f '{}'/libhdf5_fortran.so \; -print 2>/dev/null) || true
    if [[ ! -d "$HDF5_LIB" ]]; then

        # Clone HDF5 from the repository if it does not exist.
        if [ ! -d "$HDF5_DIR" ]; then
            [ -z "$HDF5_VERSION" ] && HDF5_VERSION="hdf5_2.0.0"
            git clone --depth 1 --branch $HDF5_VERSION \
                https://github.com/HDFGroup/hdf5.git $HDF5_DIR
        fi

        # Build and install HDF5
        cmake -B $HDF5_DIR/build -S $HDF5_DIR \
            --install-prefix $HDF5_DIR -DCMAKE_BUILD_TYPE=Release \
            -DCMAKE_C_COMPILER=$MPICC -DCMAKE_CXX_COMPILER=$MPICXX \
            -DCMAKE_Fortran_COMPILER=$MPIFC -DHDF5_ENABLE_PARALLEL=ON \
            -DHDF5_BUILD_FORTRAN=ON -DHDF5_ENABLE_SZIP_SUPPORT:BOOL=OFF \
            -DHDF5_BUILD_TOOLS:BOOL=ON
        cmake --build $HDF5_DIR/build/ --config Release --parallel
        cmake --install $HDF5_DIR/build/ --config Release
        rm -fr $HDF5_DIR/build
    fi

    # Add HDF5 to the environment variables
    HDF5_LIB=$(find $HDF5_DIR -type d -name 'lib*' \
        -exec test -f '{}'/libhdf5_fortran.so \; -print 2>/dev/null) || true
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
# Ensure ADIOS2 is installed, if not install it.
function find_adios2() {
    check_external_dir
    find_hdf5 $HDF5_DIR
    local pyexe
    local pyver
    local cmake_args=()

    if [[ $# -ge 1 && -n "$1" ]]; then
        ADIOS2_DIR="$1"
    elif [ -z "${ADIOS2_DIR:-}" ]; then
        return
    fi

    if ! pyexe=$(find_python_executable); then
        echo "Error: could not find python3 or python in PATH." >&2
        return 1
    fi
    pyver=$("${pyexe}" -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')

    if [[ "${ADIOS2_DIR:0:1}" != "/" && "${ADIOS2_DIR:0:1}" != "~" ]]; then
        ADIOS2_DIR="$EXTERNAL_DIR/$ADIOS2_DIR"
    fi

    ADIOS2_CONFIG="$ADIOS2_DIR/bin/adios2-config"

    if [[ ! -x "${ADIOS2_CONFIG}" ]]; then
        [ -z "${ADIOS2_VERSION:-}" ] && ADIOS2_VERSION="2.10.1"
        [ -z "${ADIOS2_ENABLE_FORTRAN:-}" ] && ADIOS2_ENABLE_FORTRAN="ON"
        [ -z "${ADIOS2_ENABLE_PYTHON:-}" ] && ADIOS2_ENABLE_PYTHON="ON"
        [ -z "${ADIOS2_ENABLE_SST:-}" ] && ADIOS2_ENABLE_SST="ON"

        if [ ! -d "$ADIOS2_DIR/.git" ]; then
            git clone --depth 1 --branch "v${ADIOS2_VERSION}" \
                https://github.com/ornladios/ADIOS2.git "$ADIOS2_DIR"
        fi

        cmake_args=(
            -DCMAKE_BUILD_TYPE=RelWithDebInfo
            -DCMAKE_INSTALL_PREFIX="$ADIOS2_DIR"
            -DCMAKE_INSTALL_PYTHONDIR="lib/python${pyver}/site-packages"
            -DADIOS2_BUILD_EXAMPLES=OFF
            -DADIOS2_USE_MPI=ON
            -DADIOS2_USE_SST="$ADIOS2_ENABLE_SST"
            -DADIOS2_USE_Python="$ADIOS2_ENABLE_PYTHON"
            -DADIOS2_USE_Fortran="$ADIOS2_ENABLE_FORTRAN"
            -DADIOS2_USE_BZip2=OFF
            -DBUILD_TESTING=OFF
            -DPython3_EXECUTABLE="$pyexe"
            -DPython_EXECUTABLE="$pyexe"
            -DPYTHON_EXECUTABLE="$pyexe"
            -DPython3_FIND_STRATEGY=LOCATION
            -DPython_FIND_STRATEGY=LOCATION
            -DCMAKE_C_COMPILER="${MPICC:-${CC:-cc}}"
            -DCMAKE_CXX_COMPILER="${MPICXX:-${CXX:-CC}}"
        )

        if [ -n "${HDF5_DIR:-}" ]; then
            cmake_args+=(
                -DADIOS2_USE_HDF5=ON
                -DHDF5_ROOT="$HDF5_DIR"
            )
        else
            cmake_args+=(
                -DADIOS2_USE_HDF5=OFF
            )
        fi

        cmake -S "$ADIOS2_DIR" -B "$ADIOS2_DIR/build" "${cmake_args[@]}"
        cmake --build "$ADIOS2_DIR/build" --parallel
        cmake --install "$ADIOS2_DIR/build"
        rm -rf "$ADIOS2_DIR/build"

        ADIOS2_CONFIG="$ADIOS2_DIR/bin/adios2-config"
    fi

    if [ ! -x "${ADIOS2_CONFIG}" ]; then
        error "ADIOS2 not found at:"
        error "\t$ADIOS2_DIR"
        error "Please set ADIOS2_DIR to the directory containing"
        error "the ADIOS2 installation."
        exit 1
    fi

    export ADIOS2_DIR="$(realpath "$ADIOS2_DIR")"
    export ADIOS2_PATH="$ADIOS2_DIR"
    export ADIOS2_FORTRAN_DIR="$ADIOS2_DIR"
    export PATH="$ADIOS2_DIR/bin:$PATH"

    [ -d "$ADIOS2_DIR/lib/pkgconfig" ] && \
        export PKG_CONFIG_PATH="$ADIOS2_DIR/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
    [ -d "$ADIOS2_DIR/lib64/pkgconfig" ] && \
        export PKG_CONFIG_PATH="$ADIOS2_DIR/lib64/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"

    [ -d "$ADIOS2_DIR/lib" ] && \
        export LD_LIBRARY_PATH="$ADIOS2_DIR/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
    [ -d "$ADIOS2_DIR/lib64" ] && \
        export LD_LIBRARY_PATH="$ADIOS2_DIR/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

    if [ -n "${pyver:-}" ] && [ -d "$ADIOS2_DIR/lib/python${pyver}/site-packages" ]; then
        export PYTHONPATH="$ADIOS2_DIR/lib/python${pyver}/site-packages${PYTHONPATH:+:$PYTHONPATH}"
    fi
    if [ -n "${pyver:-}" ] && [ -d "$ADIOS2_DIR/lib64/python${pyver}/site-packages" ]; then
        export PYTHONPATH="$ADIOS2_DIR/lib64/python${pyver}/site-packages${PYTHONPATH:+:$PYTHONPATH}"
    fi

    echo "Using ADIOS2_DIR=$ADIOS2_DIR"
    echo "Using Python=${pyexe:-<not found>}"
    echo "done"
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

        # Modify the bundled CMake files to satisfy newer toolchains.
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
    find_hdf5 $HDF5_DIR
    find_adios2 $ADIOS2_DIR
    find_parmetis $PARMETIS_DIR
    [ -n "$PFUNIT_DIR" ] && find_pfunit $PFUNIT_DIR

    # Determine the Neko installation directory
    if [[ $# -ge 1 ]]; then
        NEKO_DIR="$(realpath $1)"
    elif [ -z "$NEKO_DIR" ]; then
        NEKO_DIR="$(realpath $EXTERNAL_DIR/neko)"
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
        if [ -n "$ADIOS2_DIR" ]; then
            FEATURES+=" --with-adios2=$ADIOS2_DIR"
            FEATURES+=" --with-adios2-fortran=$ADIOS2_DIR"
        fi
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
                CC=$CC MPICC=$MPICC CFLAGS="$NEKO_CFLAGS" \
                CXX=$CXX MPICXX=$MPICXX CXXFLAGS="$NEKO_CXXFLAGS" \
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
    find_hdf5 $HDF5_DIR

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
        echo "Default path will be used: ~/tmp/external"
        EXTERNAL_DIR=~/tmp/external
    fi

    mkdir -p "$EXTERNAL_DIR"
    export EXTERNAL_DIR=$(realpath "$EXTERNAL_DIR")

}
