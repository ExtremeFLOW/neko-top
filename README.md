# Neko-TOP (Topology Optimization in Neko)

Short description of what we are trying to accomplish

Additional information can be found in the general
[documentation](./Documentation/neko-top.md).

## Quick-start

To compile and run the Taylor-Green-Vortex example from Neko, please execute the
following commands.

```sh
./setup.sh
./run.sh neko_examples/tgv
```

### Long-start

The setup script is an automated setup of Neko, Neko-TOP and other required
dependencies. The script relies in a number of environment variables, which can
be used to modify the behaviour of the system and allow the user to specify
custom install locations for the given dependencies.

| Variable           | Description                                                         | Default               |
| ------------------ | ------------------------------------------------------------------- | --------------------- |
| `NEKO_DIR`         | Location of the Neko library.                                       | external/neko         |
| `JSON_FORTRAN_DIR` | JSON-Fortran library, required dependency of Neko.                  | external/json-fortran |
| `NEK5000_DIR`      | Nek5000, primarily used for meshing and for GSLib.                  | external/Nek5000      |
| `PFUNIT_DIR`       | Unit testing library used in Neko.                                  | external/pFUnit       |
| `CUDA_DIR`         | Location of the CUDA library folders, needed for Nvidia GPU support | -                     |

These can be defined either on the command line by the user or in a `prepare.sh`
file which is loaded by the setup script if it exists in the root of Neko-TOP.
The prepare script provide a convenient way to use module systems such as
`spack` or similar to activate environments and such before compilation.

### Compilation of examples

When one have run the appropriate setup script, then one should use cmake to
compile all example files.

### Notes on linking against CUDA on WSL.

Look through the following documentations:

1. https://learn.microsoft.com/en-us/windows/ai/directml/gpu-cuda-in-wsl
2. https://docs.nvidia.com/cuda/wsl-user-guide/index.html#getting-started-with-cuda-on-wsl
3. https://developer.nvidia.com/cuda-downloads?target_os=Linux&target_arch=x86_64&Distribution=WSL-Ubuntu&target_version=2.0&target_type=deb_local

Link 1 is the microsoft description of getting started with WSL 2. Link 2 is the
NVidia guideline to how to correctly use WSL and CUDA together. Link 3 is the
link to download instructions for CUDA toolkit and drivers to WSL. Remember to
update NVidia graphics drivers on the windows side as well.

## Job execution

The run.sh script is the main driver for any example.
The run script will construct a local system for execution of any example
defined in the examples folder. Each example should contain a .case file and a
compiled version of Neko ready to be executed. This is in general provided by
our CMake setup.

The contents of Example/EXAMPLE_NAME will be copied to a temporary
log/EXAMPLE_NAME along with a job_script. This job_script can be user defined by
adding a bash file in the Scripts folder. If no specific file is found, the 
default is used.

After successful execution of Neko, the results will be moved to the Results
folder and the log folder will be cleaned.

The status script provide a simple way of probing the current status of an
example.
