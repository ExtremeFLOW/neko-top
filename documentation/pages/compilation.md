# Compilation of Neko-TOP {#compilation}
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

### External libraries

The external libraries are required to be present for the compilation of
Neko-TOP. However, these can be compiled and installed by the `setup.sh` script
if they are not already present on the system. The external libraries are:

1. Neko
2. JSON-Fortran (Required by Neko)
3. GSLib (Required by Neko)
4. Nek5000 (optional)
5. PFUnit (optional)
6. CUDA (optional)

## Quick-start compilation

To compile the library and all external dependencies, the user can run the
`setup.sh` script. This script will download and compile all dependencies and
the Neko-TOP library. The script will also compile all the advanced examples and
run the unit tests if desired.

```sh
git clone --recursive https://github.com/ExtremeFlow/Neko-TOP.git neko-top
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
| `NEK5000_DIR`      | Nek5000, primarily used for meshing and for GSLib.                   | external/Nek5000      |
| `PFUNIT_DIR`       | Unit testing library used in Neko.                                   | -                     |
| `CUDA_DIR`         | Location of the CUDA library folders, needed for Nvidia GPU support. | -                     |

These can be defined either on the command line by the user or in a `prepare.sh`
file which is loaded by the setup script if it exists in the root of Neko-TOP.
This preparation script will also be loaded by the `run.sh` script, so it is
possible to define environment variables for the execution of the examples as
well.
The prepare script provide a convenient way to use module systems such as
`spack` or similar to activate environments and such before compilation.

An example of a `prepare.sh` file is shown below:

```bash
#!/bin/bash
module load cuda/10.1

export CUDA_DIR=$CUDA_HOME
export NEKO_DIR=$HOME/neko
```

### Notes on linking against CUDA on WSL.

Look through the following documentations:

1. https://learn.microsoft.com/en-us/windows/ai/directml/gpu-cuda-in-wsl
2. https://docs.nvidia.com/cuda/wsl-user-guide/index.html#getting-started-with-cuda-on-wsl
3. https://developer.nvidia.com/cuda-downloads?target_os=Linux&target_arch=x86_64&Distribution=WSL-Ubuntu&target_version=2.0&target_type=deb_local

Link 1 is the microsoft description of getting started with WSL 2. Link 2 is the
NVidia guideline to how to correctly use WSL and CUDA together. Link 3 is the
link to download instructions for CUDA toolkit and drivers to WSL. Remember to
update NVidia graphics drivers on the windows side as well.
