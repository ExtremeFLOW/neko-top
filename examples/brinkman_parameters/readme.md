# Parameter study for Brinkman flow around a cylinder {#brinkman-parameters}

To run this example please execute the following commands from the root of
Neko-TOP:

```bash
./data/brinkman_parameters/meshing_periodic.sh
./examples/brinkman_parameters/generate_cases.sh
./run.sh brinkman_parameters
```

## Neko

The example builds and runs against current Neko `develop`.

The `brinkman` default case template uses a `distance_transform` of type
`smooth_step` rather than `step`, because Neko's `step` transform currently
produces an all-zero indicator.

The `idw` template needs the source term type `idw`, which is only available
on Neko's unmerged `feature/idw_wip` branch
(`src/source_terms/idw_source_term.f90`); it is not yet on `develop`.