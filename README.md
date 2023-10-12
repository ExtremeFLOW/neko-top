# Neko - TOP (Topology Optimization in Neko)

Short description of what we are trying to accomplish

## Compilation


### Compilation of Neko and JSON-Fortran

We have setup a few scripts which should help to setup the neko and
JSON-Fortran.

Both of these libraries are included as submodules to the system. 

### Compilation of examples

When one have run the appropriate setup script, then one should use cmake to
compile all example files.

## Job execution

The run.sh script is the main driver for any example.
The run script will construct a local system for execution of any example
defined in the examples folder. Each example should contain a .case file and a
compiled version of neko ready to be executed. This is in general provided by
our CMake setup.

The contents of Example/EXAMPLE_NAME will be copied to a temporary
log/EXAMPLE_NAME along with a job_script. This job_script can be user defined by
adding a bash file in the Scripts folder. If no specific file is found, the
default is used.

After successful execution of neko, the results will be moved to the Results
folder and the log folder will be cleaned.

The status script provide a simple way of probing the current status of an
example.
