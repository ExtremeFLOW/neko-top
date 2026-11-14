## `unsteady_mixer`

Six topology-optimization cases for the same Brinkman-mapped mixer geometry as
[Static mixer](@ref static-mixers), at Re = 200: the `steady_200*`/
`unsteady_200*` pairs select the steady or unsteady adjoint formulation, and
the `_short`/`_short_small` suffixes shorten the run time and/or coarsen the
mesh.

### GPU-only velocity solver

All six `*.case` files set the velocity solver to `fused_cg`:

```json
"velocity_solver": {
    "type": "fused_cg",
    "preconditioner": { "type": "jacobi" },
    ...
}
```

`fused_cg` is implemented only for the CUDA and HIP backends and has no CPU
fallback, so a CPU build stops at once during fluid setup:

```
*** ERROR: FusedCG only supported for CUDA/HIP ***
```

For a CPU run, change the velocity solver type to `cg` in a working copy of
the case file instead:

```json
"velocity_solver": {
    "type": "cg",
    "preconditioner": { "type": "jacobi" },
    ...
}
```

The tracked case files deliberately keep `fused_cg` - they exist to exercise
the GPU-tuned solver - so make this substitution in a copy rather than editing
them in place.
