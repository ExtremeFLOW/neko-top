# Neko-TOP (Topology Optimization in Neko)

Short description of what we are trying to accomplish

Additional information can be found in the general
[documentation](./Documentation/neko-top.md).

## Quick-start

```sh
./setup_local.sh
```

### Compilation of Neko and JSON-Fortran

We have setup a few scripts which should help to setup the neko and
JSON-Fortran.

Both of these libraries are included as submodules to the system. 

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

### Notes on the InnoTop machine

On the InnoTop machine a spack environment have been constructed which can be
loaded with:

```bash
spack env activate neko-top
```

This will load all required packages from spack, after which the
`setup_local.sh` script should work as intended.

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
