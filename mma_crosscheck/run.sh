#!/bin/sh
# Cross-check Neko-TOP's CPU MMA ("dip" subsolver) against Niels Aage's MMA.cc
# (topopt_in_petsc) on the 1D cantilever beam from tests/regression/mma.
#
# usage: ./run.sh /path/to/neko-top /path/to/topopt_in_petsc [iterations]
# needs: mpif90, mpicxx, LAPACK/BLAS, python3 with numpy. PETSc is used when
#        pkg-config finds it (module name PETSc), otherwise the reference is
#        built against petsc_shim/petsc.h (Vec API only, real MPI).
set -e
NEKOTOP=$(cd "${1:?path to neko-top}" && pwd)
PETSCREF=$(cd "${2:?path to topopt_in_petsc}" && pwd)
NIT=${3:-20}
HERE=$(cd "$(dirname "$0")" && pwd)
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
FFLAGS="-O2 -cpp -ffree-line-length-none"
mkdir -p "$HERE/work" && cd "$HERE/work"

# Neko-TOP: mma.f90 + mma_cpu.f90 from the given tree against stub Neko modules
rm -rf build_neko && mkdir build_neko && cd build_neko
cp "$HERE/neko_stubs.f90" "$HERE/mma_device_stub.f90" "$HERE/beam_driver.f90" .
cp "$NEKOTOP/sources/mma/mma.f90" "$NEKOTOP/sources/mma/bcknd/cpu/mma_cpu.f90" .
for f in neko_stubs mma mma_cpu mma_device_stub beam_driver; do
    mpif90 $FFLAGS -c $f.f90
done
mpif90 -o ../beam_neko beam_driver.o mma.o mma_cpu.o mma_device_stub.o \
    neko_stubs.o -llapack -lblas
cd ..

# Reference: the unmodified MMA.cc
if pkg-config --exists PETSc 2>/dev/null; then
    PETSC_FLAGS="$(pkg-config --cflags PETSc)"
    PETSC_LIBS="$(pkg-config --libs PETSc)"
else
    echo "PETSc not found by pkg-config, using petsc_shim/petsc.h"
    PETSC_FLAGS="-I$HERE/petsc_shim"
    PETSC_LIBS=""
fi
mpicxx -O2 -I"$PETSCREF" $PETSC_FLAGS "$HERE/beam_ref.cc" "$PETSCREF/MMA.cc" \
    -o beam_ref $PETSC_LIBS

# Parameters of tests/regression/mma/cases/1d_beam_dip.case (= MMA.cc defaults)
P="0.5 1.2 0.7 1000.0"
#            name       iterations  move limit  a
./beam_ref   ref        "$NIT"      -1.0        0.0
./beam_ref   ref_ml     "$NIT"      0.2         0.0
./beam_ref   ref_a1     "$NIT"      -1.0        1.0
#                       subsolver name     iterations  asymptotes, c  move limit  a
./beam_neko             dip  neko          "$NIT"      $P             -1.0        0.0
./beam_neko             dip  neko_ml       "$NIT"      $P             0.2         0.0
./beam_neko             dip  neko_a1       "$NIT"      $P             -1.0        1.0
mpirun -n 3 ./beam_neko dip  neko_np3      "$NIT"      $P             -1.0        0.0

C="python3 $HERE/compare.py -q"
$C neko     ref    "dip vs MMA.cc"
$C neko_ml  ref_ml "dip vs MMA.cc, move limit 0.2"
$C neko_a1  ref_a1 "dip vs MMA.cc, a = 1"
$C neko_np3 ref    "dip on 3 ranks vs MMA.cc"
