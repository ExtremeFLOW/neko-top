# Installation {#installation}
\tableofcontents

Neko-TOP is an extension of the Neko library, and as such it requires the Neko
library to be installed. The Neko library also has a number
of dependencies, such as the JSON-Fortran library, which is used for reading and
writing JSON files, and the PFUnit library, which is used for unit testing.

Here we provide a brief overview of the dependencies and how to compile the
library. For more detailed information, please see the documentation for the
individual libraries. In particular, the Neko library has a detailed
documentation on how to compile and install the library, with support for
various accelerators. [Neko homepage](https://neko.cfd).


## Dependencies

We group the dependencies into two categories: system dependencies and external
libraries. The system dependencies are the dependencies that are required to
compile the Neko library and the Neko-TOP library. The external libraries are
the libraries that are required to run the examples and the tests.

### System dependencies

The system dependencies are the dependencies that are required to compile the
Neko library and the Neko-TOP library. The system dependencies are:

1. A Fortran 2008 compliant compiler
2. CMake
5. A C compiler
6. A C++ compiler
7. A MPI library
8. A BLAS/LAPACK library
9. Autotools
10. PKG-Config

### External libraries {#installation-external}

The external libraries are required to be present for the compilation of
Neko-TOP. However, these can be compiled and installed by the `setup.sh` script
if they are not already present on the system. The external libraries are:

1. Neko
2. JSON-Fortran (Required by Neko)
3. GSLib (Required by Neko)
4. Nek5000 (optional)
5. PFUnit (optional)
6. CUDA (optional)

## Quick-start compilation {#installation-quick}

To compile the library and all external dependencies, the user can run the
`setup.sh` script. This script will download and compile all dependencies and
the Neko-TOP library. The script will also compile all the advanced examples and
run the unit tests if desired.

```sh
git clone https://github.com/ExtremeFlow/Neko-TOP.git neko-top
cd neko-top
./setup.sh
```

### Setup script

The `setup.sh` script is an automated setup of Neko-TOP along with Neko and other
dependencies. The script relies in a number of environment variables, which can
be used to modify the behaviour of the system and allow the user to specify
custom install locations for the given dependencies.

| Variable           | Description                                                          | Default               |
| ------------------ | -------------------------------------------------------------------- | --------------------- |
| `NEKO_DIR`         | Location of the Neko library.                                        | external/neko         |
| `JSON_FORTRAN_DIR` | JSON-Fortran library, required dependency of Neko.                   | external/json-fortran |
| `ADIOS2_DIR`       | ADIOS2 installation; setting it enables ADIOS2 support.              | -                     |
| `NEKO_LIBS`        | Extra linker flags for Neko, passed to Neko's `configure` script.    | -                     |
| `NEK5000_DIR`      | Nek5000, primarily used for meshing and for GSLib.                   | external/Nek5000      |
| `PFUNIT_DIR`       | Unit testing library used in Neko.                                   | -                     |
| `CUDA_DIR`         | Location of the CUDA library folders, needed for Nvidia GPU support. | -                     |
| `CUDA_ARCH`        | CUDA architecture for GPU builds (for example `80` for sm_80).       | Required for CUDA     |

These can be defined either on the command line by the user or in a
`prepare.env` file which is loaded by the setup script if it exists in the root
of Neko-TOP. This preparation script will also be loaded by the `run.sh` script,
so it is possible to define environment variables for the execution of the
examples as well. The prepare script provide a convenient way to use module
systems such as `spack` or similar to activate environments and such before
compilation.

An example of a `prepare.env` file is shown below:

```bash
#!/bin/bash
module load cuda/10.1

export CUDA_DIR=$CUDA_HOME
export CUDA_ARCH=80
export NEKO_DIR=$HOME/neko
```

### ADIOS2-enabled Python workflow

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
export ADIOS2_DIR=adios2
./setup.sh
```

The relative value `adios2` installs ADIOS2 in `external/adios2`; an absolute
path can be used to install it elsewhere.

If the active Python environment or MPI toolchain changes after ADIOS2 has been
built, rebuild from ADIOS2 onward so that `mpi4py`, ADIOS2, Neko, and
Neko-TOP all agree on the same Python and MPI stack.

\warning Some Fortran MPI wrappers do not automatically link the C++ runtime
required by Neko's ADIOS2 C++ interface. If the Neko link step reports undefined
C++ symbols such as `std::`, `__gxx_personality_v0`, or `adios2::`, set
`NEKO_LIBS` in `prepare.env` before rebuilding Neko. It must include both the
C++ runtime and ADIOS2 libraries after `libneko.a`. For GNU compilers with an
ADIOS2 installation selected by `ADIOS2_DIR`, use:

```bash
adios2_config="$EXTERNAL_DIR/$ADIOS2_DIR/bin/adios2-config"
adios2_link_libs=""
for adios2_flag in $("$adios2_config" --cxx-libs); do
    case "$adios2_flag" in
    */lib*.so*|*/lib*.a)
        adios2_lib_dir=$(dirname "$adios2_flag")
        adios2_lib_name=$(basename "$adios2_flag")
        adios2_lib_name=${adios2_lib_name#lib}
        adios2_lib_name=${adios2_lib_name%%.so*}
        adios2_lib_name=${adios2_lib_name%%.a}
        adios2_link_libs="$adios2_link_libs -L$adios2_lib_dir -l$adios2_lib_name"
        ;;
    *) adios2_link_libs="$adios2_link_libs $adios2_flag" ;;
    esac
done
NEKO_LIBS="$adios2_link_libs -lstdc++"
```

`NEKO_LIBS` is passed unchanged to Neko's `configure` script. Use the
equivalent compiler-runtime flag for a non-GNU toolchain.

For CUDA builds, `CUDA_ARCH` must be explicitly specified before running
`setup.sh` (for example `export CUDA_ARCH=80`).

Additional examples of the preparation script for specific systems can be seen
in the section on clusters: \subpage clusters.


### Notes on linking against CUDA on WSL.

Look through the following documentations:

1. https://learn.microsoft.com/en-us/windows/ai/directml/gpu-cuda-in-wsl
2. https://docs.nvidia.com/cuda/wsl-user-guide/index.html#getting-started-with-cuda-on-wsl
3. https://developer.nvidia.com/cuda-downloads?target_os=Linux&target_arch=x86_64&Distribution=WSL-Ubuntu&target_version=2.0&target_type=deb_local

Link 1 is the microsoft description of getting started with WSL 2. Link 2 is the
NVidia guideline to how to correctly use WSL and CUDA together. Link 3 is the
link to download instructions for CUDA toolkit and drivers to WSL. Remember to
update NVidia graphics drivers on the windows side as well.
