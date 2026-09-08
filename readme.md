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
