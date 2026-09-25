# Running examples {#examples}
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

## Command-line options {#examples-options}

The run script supports the following switches:

- `-a`, `--all`: Queue all discovered examples.
- `-c`, `--clean`: Remove existing logs for queued examples before setup.
- `-d`, `--delete`: Remove all stored results and logs (with confirmation).
- `-h`, `--help`: Print inline help.
- `-n`, `--neko`: Run examples from the Neko repository (`$NEKO_DIR/examples`) instead of Neko-TOP examples.
- `-s`, `--submit CLUSTER`: Submit prepared examples to a batch system.
- `--dry-run`: Prepare logs/job files and stop before execution/submission.
- `-r`, `--re-run`: Re-run examples even if results already exist.
- `-p`, `--procs N`: Set number of MPI processes for local execution.
- `--sequential`: Chain submitted examples through scheduler dependencies.
- `--njobs N`: Submit `N` jobs per queued example in sequence.

Examples:

```sh
# Run one local Neko-TOP example
./run.sh passive_scalar

# Run Neko examples from NEKO_DIR
./run.sh --neko tgv

# Submit to a cluster with three sequential jobs per example
./run.sh --submit MN5 --njobs 3 --sequential static_mixer

# Prepare everything without running/submitting
./run.sh --submit LUMI-G --dry-run thermal_mix_channel
```

## Execution of examples {#examples-running}

The [run.sh](../../run.sh) script is the main driver for all examples. It
constructs a local run directory in `logs/` for each queued example and then
either executes locally or submits with a scheduler.

Input selection supports names, simple patterns, and direct paths to example
folders/case files. Case files can be either `.case` or `.json`. When multiple
files match, all valid matches are queued.

Before execution/submission, run setup does the following for each queued
example:

- Creates `logs/EXAMPLE_NAME`.
- Copies case files and auxiliary files from the example folder.
- Keeps `.nmsh` files as symbolic links (to avoid large mesh copies).
- Links shared `data/` and `data_local/` directories.
- Copies runtime helper scripts and (for cluster mode) a selected job script.

If a previous run failed, existing log files are archived under `run_XX`
subfolders and a fresh run state is prepared.

For local Neko-TOP examples, the script also builds the `Examples` target from
the project build tree before running.

```sh
./run.sh --neko tgv
```

### Submission of jobs to a cluster.

Cluster submission is enabled with `-s` / `--submit`. For each queued example,
`run.sh` looks up a batch script under:

`scripts/jobscripts/CLUSTER/`

and first tries a case-specific script (`path/to/case.sh`), then recursively
falls back to `default.sh` in parent folders.

Submission behavior:

- If `--njobs N` is set, `N` jobs are submitted per example in sequential order.
- If `--sequential` is set, dependencies are also chained across examples, so
  the next example starts after previously submitted jobs complete.
- Existing queued jobs with the same name are skipped.

As for local runs, staged input files live in `logs/`, and completed outputs are
placed under `results/`.

```sh
./run.sh --submit MN5 --neko tgv
```

Currently the following clusters are supported:

- MN5: The MareNostrum 5 cluster at BSC. (Utilizes the SLURM scheduler)
- LUMI-C: The LUMI CPU partition. (Utilizes the SLURM scheduler)
- LUMI-G: The LUMI GPU partition. (Utilizes the SLURM scheduler)
- DTU: The local cluster at the Technical University of Denmark. (Utilizes the
  LSF10 scheduler; use matching job scripts in `scripts/jobscripts/DTU`)

Unsupported clusters will result in an error message, but, we recommend to
execute the `run.sh` script with the `--dry-run` switch to organize the job
files. Then one can manually submit the job to the cluster.

## Adding examples. {#examples-new}

To construct new examples, place a folder in the `examples` folder. Each example
should be self-contained with case files, the nmsh required. Any additional
source files should ideally be placed in that folder as well.

Additionally the following constraints on the structure of the example folder
are required:

- The example is required to contain at most a single `run.sh` script. This
  script indicate the root of the example.
- The example may contain any number of case files, but if there is no `run.sh`
  script, all of them should be placed in the same subfolder of the example.
  Otherwise, each folder containing a case file will be considered an example.

## Case files and meshes. {#examples-data}

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

## Advanced example setups. {#examples-advanced}

Elaborate examples can be constructed by using the CMakeLists.txt file to
specify the source files and the driver type. The driver type is used to
determine how the example should be compiled and executed. The CMakeLists.txt in
the example folder must call the `build_example` function to setup the example. 

Current driver types:

- "neko":        Pure neko with no user defined source files.
- "neko-user":   Equivalent to using makeneko to generate the executable.
- "topopt":      Topology optimization driver defined in the Neko-TOP library.
- "topopt-user": Topology optimization driver with loading of a Neko user
                 module.
- "custom":      Use the driver.f90 provided in the current example folder or
                 read the cmake variable DRIVER for a fully custom driver.

Additional source files can be added to the example by setting the EXTRA_SOURCES
variable in the CMakeLists.txt file. Preferably these files should be specified
with a full path or using the CMake variables.

Please note that the final entry in the CMakeLists.txt file should be the
call to the `build_example` function. This function will take care of the
rest of the setup.

Example of a CMakeLists.txt file:

```cmake
set(DRIVER_TYPE "custom")

# If example require a custom driver, define the source file here.
# set(DRIVER PATH/TO/CUSTOM/DRIVER.f90)

# Extra sources that need to be compiled into the driver.
# set(EXTRA_SOURCES "")

build_example()

```
