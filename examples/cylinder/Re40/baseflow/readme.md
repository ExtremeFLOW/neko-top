# Re40 fixed point Newton solve

Runnable translation of the Re40 steady cylinder baseflow example.
Run `./prepare.sh` first to generate `newton.nmsh` with gmsh/gmsh2nek/rea2nbin.
The initial Newton guess is defined by `case.fluid.initial_condition`.

The driver stops after the Newton fixed point and writes
`viz/baseflow_newton_Re40*.fld`. Run the separate Re40 `LNS` or `adjoint`
example to compute eigenvalues.
