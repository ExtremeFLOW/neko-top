# Neko-TOP (Topology Optimization in Neko)

[![Documentation](https://github.com/ExtremeFLOW/neko-top/actions/workflows/documentation.yml/badge.svg?branch=develop&event=push)](https://github.com/ExtremeFLOW/neko-top/actions/workflows/documentation.yml)

The Neko-TOP library is an extension of the Neko library, which is a high-order
spectral element solver. The Neko-TOP library is designed to solve topology
optimization problems using an immersed boundary method. The library is written
in Fortran and is designed to be used in combination with the Neko library.

Details can be found in the documentation on github pages,
https://extremeflow.github.io/neko-top/.

## Dependencies

The Neko-TOP library is dependent on the following libraries:

- Fortran 2008  
    We assume gfortran, use `FC` environment variable to override.
- MPI: We tested with OpenMPI 3.1.
- [Neko](https://github.com/ExtremeFlow/Neko)
    - [![Neko: Neko-TOP branch](https://github.com/ExtremeFLOW/neko-top/actions/workflows/neko_neko-top.yml/badge.svg?event=schedule)](https://github.com/ExtremeFLOW/neko-top/actions/workflows/neko_neko-top.yml)
    - [![Neko: Develop branch](https://github.com/ExtremeFLOW/neko-top/actions/workflows/neko_develop.yml/badge.svg?event=schedule)](https://github.com/ExtremeFLOW/neko-top/actions/workflows/neko_develop.yml)
    - [![Neko: Master branch](https://github.com/ExtremeFLOW/neko-top/actions/workflows/neko_master.yml/badge.svg?event=schedule)](https://github.com/ExtremeFLOW/neko-top/actions/workflows/neko_master.yml)

- [JSON-Fortran](https://github.com/jacobwilliams/json-fortran).
- CUDA; Optional for GPU acceleration in Neko.

The Neko-TOP library is also dependent on the following libraries for testing:

- pFUnit         (Built through CMake if unavailable)

## Quick-start compilation

To compile the library and all external dependencies, the user can run the
`setup.sh` script. This script will download and compile all dependencies and
the Neko-TOP library. The script will also compile all the advanced examples and
run the unit tests if desired.

```sh
git clone https://github.com/ExtremeFlow/Neko-TOP.git neko-top
cd neko-top
./setup.sh
```

For CUDA builds, `CUDA_ARCH` must be set explicitly before running setup, for
example:

```sh
export CUDA_ARCH=80
./setup.sh -d CUDA
```

## Python and MPI build order

For the ADIOS2-enabled Python workflow, the build order matters. The intended
order is:

1. Create and activate a fresh Python environment.
2. Install `mpi4py` into that environment with the MPI compiler wrapper that
   will be used for the rest of the build.
3. Build ADIOS2 against that same active Python environment.
4. Build Neko against that ADIOS2 installation.
5. Build Neko-TOP on top of that Neko build.

In practice, `./setup.sh` performs steps 3-5, so the critical requirement is
that steps 1-2 are done first in the shell where setup is invoked. This avoids
mixing one Python environment for `mpi4py` with another Python environment for
ADIOS2 and the runtime scripts.

The recommended workflow is therefore:

```sh
python -m venv PATH_TO_ENV
source PATH_TO_ENV/bin/activate
MPICC=mpicc python -m pip install --no-binary=mpi4py mpi4py
./setup.sh
```

If the active Python environment or MPI toolchain changes after ADIOS2 has been
built, rebuild from ADIOS2 onward so that `mpi4py`, ADIOS2, Neko, and
Neko-TOP all agree on the same Python and MPI stack.

## Example execution

The run.sh script is the main driver for managing example execution. The run
script will construct a temporary folder system in `logs` for execution of any
example defined in the examples folder. The can then run a given example or list
of examples by invoking the run script with the name of the example as an
argument.

``` sh
./run.sh EXAMPLE_NAME
```

After successful execution of Neko, the results will be moved to the `results`
folder and the log folder will be cleaned.

The `status.sh` script provide a simple way of probing the current status of an
example.
