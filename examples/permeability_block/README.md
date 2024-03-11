# Permeability block example {#permeability_block}

The permeability block example is a simple example of a flow through a porous
media. The example is designed to show how to setup a simple case file and how
to run a simulation using the Neko library. The example is located in the
`examples/permeability_block` folder.

The example consists of a meshed square duct where a wall is placed inside. The
wall is defined as a porous media with a given permeability. The flow is driven
by a pressure difference between the inlet and the outlet.

The permeability is added as a volume force with a given value proportional to
the local velocity. The permeability is defined as a scalar value and is added
to the case file as a parameter.

$$
    f(x) = - perm * u(x)
$$

1. `permeability_3`: The permeability is set to 1000.
2. `permeability_5`: The permeability is set to 100000.
3. `permeability_wall`: The mesh is modified to include
   a wall boundary condition instead of a porous media.

Please note that due to size constraints the mesh is not included in the
repository. However, the Cubit journal file is included in the
`data` folder. The mesh can be generated using
the journal file and the Cubit software, along with the `Exo2Nek` plugin from
Nek5000. Please see the meshing documentation for details on how to generate the
mesh: \ref meshing.
