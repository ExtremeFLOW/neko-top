## Examples

## Example guide:

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

### Case files and meshes.

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
