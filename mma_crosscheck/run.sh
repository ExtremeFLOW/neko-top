#!/bin/sh
# Cross-check Neko-TOP's CPU MMA ("dip" subsolver) against Niels Aage's MMA.cc
# (topopt_in_petsc) on the 1D cantilever beam from tests/regression/mma.
#
# usage: ./run.sh /path/to/neko-top /path/to/topopt_in_petsc [iterations] [move_limit]
# needs: mpif90, mpicxx, LAPACK/BLAS, PETSc visible to pkg-config (module name PETSc)
#        e.g. export PKG_CONFIG_PATH=$PETSC_DIR/$PETSC_ARCH/lib/pkgconfig
set -e
NEKOTOP=$(cd "${1:?path to neko-top}" && pwd)
PETSCREF=$(cd "${2:?path to topopt_in_petsc}" && pwd)
NIT=${3:-20}
MOVE=${4:--1.0}          # <= 0: no move limit (as in 1d_beam_dip.case)
HERE=$(cd "$(dirname "$0")" && pwd)
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
FFLAGS="-O2 -cpp -ffree-line-length-none"
mkdir -p "$HERE/work" && cd "$HERE/work"

build_neko() {   # $1 = variant name, $2 = preprocessor flags ("" = upstream)
    rm -rf "build_$1" && mkdir "build_$1" && cd "build_$1"
    cp "$HERE/neko_stubs.f90" "$HERE/mma_device_stub.f90" "$HERE/beam_driver.f90" .
    cp "$NEKOTOP/sources/mma/mma.f90" "$NEKOTOP/sources/mma/bcknd/cpu/mma_cpu.f90" .
    if [ -n "$2" ]; then patch -s mma_cpu.f90 < "$HERE/aage_variants.patch"; fi
    for f in neko_stubs mma mma_cpu mma_device_stub beam_driver; do
        mpif90 $FFLAGS $2 -c $f.f90
    done
    mpif90 -o "../beam_neko_$1" beam_driver.o mma.o mma_cpu.o mma_device_stub.o \
        neko_stubs.o -llapack -lblas
    cd ..
}

build_neko upstream ""
build_neko aage "-DAAGE_REG -DAAGE_TERM"
mpicxx -O2 -I"$PETSCREF" $(pkg-config --cflags PETSc) "$HERE/beam_ref.cc" \
    "$PETSCREF/MMA.cc" -o beam_ref $(pkg-config --libs PETSc)

# Parameters of tests/regression/mma/cases/1d_beam_dip.case
./beam_ref ref "$NIT" "$MOVE"
./beam_neko_upstream dip neko_upstream "$NIT" 0.5 1.2 0.7 1000.0 "$MOVE"
./beam_neko_aage dip neko_aage "$NIT" 0.5 1.2 0.7 1000.0 "$MOVE"

python3 "$HERE/compare.py" neko_upstream ref "upstream Neko-TOP vs MMA.cc"
python3 "$HERE/compare.py" neko_aage ref "Neko-TOP with Aage p/q coefficients and stopping rules vs MMA.cc"
