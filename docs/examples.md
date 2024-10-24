# Examples {#examples}
\tableofcontents

The execution of examples have been automated to allow for easy testing and
development of new features. The examples are located in the `examples` folder
and are self-contained with case files and the nmsh files required (or details
on how to get them). The examples are compiled and executed using the `run.sh`
script.

For inline help on the `run.sh` script, use the `-h` or `--help` switch.

```sh
./run.sh [-h] [--help]
```

## List of examples

1. \subpage permeability_block
2. \subpage easy-E

## Execution of examples

The [run.sh](../../run.sh) script is the main driver for all examples.
The run script will construct a local system for execution of any example
defined in the examples folder. Each example should at least contain a .case
file with the specifications to Neko. The run script will copy the contents of
the example folder to a temporary folder and execute the example from there.

The contents of `example/EXAMPLE_NAME` will be copied to a temporary
`log/EXAMPLE_NAME`. If the user wish to execute one of the Neko examples, then the
`run.sh` script can be invoked with the name of the example as an argument along
the switch `-n` or `--neko`.

```sh
./run.sh --neko tgv
```

### Submission of jobs to a cluster.

Currently the `run.sh` script supports submission of jobs to a cluster using the
`-s` or `--submit` switch. The script will submit the job to the cluster using a
batch script associated with the desired cluster. The batch script is located in
the folder `scripts/jobscripts` under a folder associated with the cluster name.

The script will copy the contents of the example folder to a temporary folder as
is done when running the example locally. The job will be submitted to the
cluster and the output will be stored in the `results` folder after completed
execution.

```sh
./run.sh --submit MN5 --neko tgv
```

Currently the following clusters are supported:

- MN5: The MareNostrum 5 cluster at BSC. (Utilizes the SLURM scheduler)
- DTU: The local cluster at the Technical University of Denmark. (Utilizes the
  LSF10 scheduler)

Unsupported clusters will result in an error message, but, we recommend to
execute the `run.sh` script with the `--dry-run` switch to organize the job
files. Then one can manually submit the job to the cluster.

## Adding examples

To construct new examples, place a folder in the `examples` folder. Each example
should be self-contained with case files, the nmsh required. Any additional
source files should ideally be placed in that folder as well.

Additionally the following constraints on the structure of the example folder
are required:

- The example is required to only contain a single `run.sh` script. This script
  indicate the root of the example.
- The example may contain any number of case files, but if there is no `run.sh`
  script, all of them should be placed in the root of the example folder.
  Otherwise, each folder containing a case file will be considered an example.

## Case files and meshes.

When running an example a link is made to the folders `data` and `data_local`
both sitting at the root of Neko-TOP. The `data` folder contain some official
meshes used in the predefined examples. The `data_local` allow the user to save
meshes to a folder accessible to all examples to reduce need for redundancy.
This folder can include anything needed by the examples to run. In general it is
used for local copies of meshes which can be accessed by all examples to avoid
redundancy and copying massive folders around.

Additionally, we avoid copying meshes contained in the examples by creating a
link to any `.nmsh` files contained in the example folder. This is done to avoid
copying large meshes around and to allow for easy access to the meshes.

## Advanced example setups

Elaborate examples can be constructed by using the CMakeLists.txt file to
specify the source files and the driver type. The driver type is used to
determine how the example should be compiled and executed. The CMakeLists.txt in
the example folder must call the `build_example` function to setup the example. 

Current driver types:

- "default": Pure neko with no user defined source files.
- "user":    Equivalent to using makeneko to generate the executable.
- "topopt":  Topology optimization driver defined in the Neko-TOP library.
- "custom":  Use the driver.f90 provided in the current example folder or read
             the cmake variable DRIVER.

Additional source files can be added to the example by setting the EXTRA_SOURCES
variable in the CMakeLists.txt file. Preferably these files should be specified
with a full path or using the CMake variables.

Please note that the final entry in the CMakeLists.txt file should be the
call to the `build_example` function. This function will take care of the
rest of the setup.

Example of a CMakeLists.txt file:

```cmake
set(DRIVER_TYPE "default")

# If example require a custom driver, define the source file here.
# set(DRIVER PATH/TO/CUSTOM/DRIVER.f90)

# Extra sources that need to be compiled into the driver.
# set(EXTRA_SOURCES "")

build_example()

```
