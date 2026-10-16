# Re60 artificially stabilized baseflow

Artificially stabilized steady cylinder baseflow setup at Re60.
Run `./prepare.sh` first to generate `newton.nmsh` with gmsh/gmsh2nek/rea2nbin.
The initial Newton guess is defined by `case.fluid.initial_condition`.

By default this driver stops after the Newton fixed point and writes
`viz/baseflow_newton_Re60*.fld`. Set
`case.lightkrylov.eigensolver.enabled` to `true` in `1cyl_newton.case` to run
the same eigensolver stage used by `examples/newton`.
