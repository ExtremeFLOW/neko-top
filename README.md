# Neko-TOP (Topology Optimization in Neko)

Short description of what we are trying to accomplish

Additional information can be found in the doxygen based documentation in the
`documentation` folder. To view the documentation, please open the
[documentation/html/index.html](./documentation/html/index.html) file in a web
browser. 

## Quick-start

To compile and run the Taylor-Green-Vortex example from Neko, please execute the
following commands. Please note that any custom paths can be read from the
environment or loaded from a `prepare.sh` script.

```sh
./setup.sh
./run.sh neko_examples/tgv
```

### Compilation of examples

When one have run the appropriate setup script, then one should use cmake to
compile all example files.

## Job execution

The run.sh script is the main driver for any example.
The run script will construct a local system for execution of any example
defined in the examples folder. Each example should contain a .case file and a
compiled version of Neko ready to be executed. This is in general provided by
our CMake setup.

The contents of `example/EXAMPLE_NAME` will be copied to a temporary
`log/EXAMPLE_NAME` along with a job_script. This job_script can be user defined by
adding a bash file in the Scripts folder. If no specific file is found, the 
default is used. If the user wish to execute one of the Neko examples, then the
`run.sh` script can be invoked with the name of the example as an argument along
the switch `-n` or `--neko`.

```sh
./run.sh --neko tgv
```

After successful execution of Neko, the results will be moved to the `results`
folder and the log folder will be cleaned.

The `status.sh` script provide a simple way of probing the current status of an
example.
