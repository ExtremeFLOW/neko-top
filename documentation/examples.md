## Examples

To construct new examples, place a folder in the Examples root. Each example
should be self-contained with a case file, the nmsh required and a simple
CMakeLists.txt. Any additional source files should ideally be placed in that
folder as well, however, the EXTRA_SOURCES variable should just point to the
correct files.

The CMakeLists.txt should follow the structure of the provided template.
Therefore, one should set the extra sources required along with the type of
driver needed. An error is thrown if a non-supported driver type is requested.

Current driver types:

- "user":    Equivalent to using makeneko to generate the executable.
- "custom":  Use the driver.f90 provided in the current example folder or read
             the cmake variable DRIVER.
- "default": Pure neko with no user defined source files.
- "topopt":  Topology optimization driver defined in the Neko-TOP library.

## Execution

The run.sh script is the main driver for any example. The run script will
construct a local system for execution of any example defined in the examples
folder. Each example should contain a .case file and a compiled version of Neko
ready to be executed. This is in general provided by our CMake setup.

The contents of `example/EXAMPLE_NAME` will be copied to a temporary
`log/EXAMPLE_NAME` along with a job_script. This job_script can be user defined
by
adding a bash file in the Scripts folder. If no specific file is found, the
default is used. If the user wish to execute one of the Neko examples, then the
`run.sh` script can be invoked with the name of the example as an argument along
the switch `-n` or `--neko`.

```sh
./run.sh --neko tgv
```

## Case files and meshes.

When running an example a link is made to the local data folder if it exists.
This folder can include anything needed by the examples to run. In general it is
used for local copies of meshes which can be accessed by all examples to avoid
redundancy and copying massive folders around.

#### Meshing

A system have been setup to allow the user to create new meshes using Cubit.
A journal fle should be placed in the example folder to which it applies. The
`mesh.sh` script can then be used to generate the mesh and convert it to a Neko
supported format. It is assumed the mesh is exported as an exodus file with same
basename as the journal file.

To avoid cluttering the example folder with unnecessary and large files, the
generated mesh is placed in the [data](data/) folder under the same
folder structure as the example that was run. This is checked if no .nmsh file
is found in the example folder by run.sh.
