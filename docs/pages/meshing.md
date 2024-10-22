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
