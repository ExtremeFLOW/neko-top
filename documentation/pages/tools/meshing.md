# Meshing {#meshing}

The Neko-TOP library provides a meshing tool that allows users to create meshes
for their cases. The meshing tool is designed to work with the Neko-TOP library
and the Neko library.

The meshing tool is a wrapper around a variety of meshing tools, such as
Cubit, Gmsh, and Nek5000. It allows users to create meshes in a variety of
formats, including Exodus II, Gmsh, and Nek5000, before ultimately converting it
to the Neko supported format.

## Supported meshing tools
The meshing tool supports the following meshing tools:
- [Cubit](https://www.3ds.com/products-services/catia/products/cubit/)
- [Gmsh](http://gmsh.info/)
- [genmeshbox](https://github.com/ExtremeFLOW/neko/tree/develop/contrib/genmeshbox)

## Usage
The meshing tool command line interface depends on the selected meshing tool.

### genmeshbox
The `genmeshbox` tool is a simple meshing tool that generates a mesh for a
rectangular box. It is a part of the Neko library and can be used to
generate meshes for simple geometries. The command line interface is as follows:

```bash
mesh.sh --box [options] x0 x1 y0 y1 z0 z1 nx ny nz [px py pz]
```

Where:
- `x0`, `x1`, `y0`, `y1`, `z0`, `z1` are the coordinates of the box
- `nx`, `ny`, `nz` are the number of elements in the x, y, and z directions
- `px`, `py`, `pz` determine if a periodic boundary condition should be applied
  in the x, y, and z directions (optional, default is no periodicity)

### Other meshing tools
The other meshing tools, such as Cubit and Gmsh, rely on input files that
define the geometry and mesh parameters. The command line interface for these
tools is typically more complex and requires the user to create input files
that specify the geometry, mesh size, and other parameters. The Neko-TOP library
provides a set of example input files which can be used as a starting point for
creating your own input files.

```bash
mesh.sh [options] input_file
```

By default the input file is expected to be in the data directory of the
Neko-TOP library and an output file will be created in the data_local
directory.