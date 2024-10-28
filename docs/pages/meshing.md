# Meshing {#meshing}

A system have been setup to allow the user to create new meshes using Cubit.
A journal fle should be placed in the example folder to which it applies. The
`mesh.sh` script can then be used to generate the mesh and convert it to a Neko
supported format. It is assumed the mesh is exported as an exodus file with same
basename as the journal file.

To avoid cluttering the example folder with unnecessary and large files, the
generated mesh is placed in the [data](data/) folder under the same
folder structure as the example that was run. This is checked if no .nmsh file
is found in the example folder by run.sh.

> [!note]
> The meshing tool does currently not support periodic boundary conditions.
> These must be manually added to the mesh file. We recommend manually adding
> periodic boundary conditions to the mesh file using the `Nek5000` tools,
> `exo2nek` and `gmsh2nek`. After which the meshing tool can manage the output
> `re2` file.

## Supported formats

The following mesh formats are supported by our meshing tool:

- Exodus II (.e .exo)
- Gmsh      (.geo .msh)
- Nek5000   (.rea .re2)

