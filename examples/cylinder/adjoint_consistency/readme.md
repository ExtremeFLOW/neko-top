# Adjoint consistency check

This mini example checks whether the transient LNS propagator and its adjoint
satisfy

```text
<A x, y> = <x, A* y>
```

for a small number of fitted Floquet timesteps. The timestep is fitted exactly
as in the full Floquet examples, but the interval is truncated to
`case.lightkrylov.adjoint_test.number_of_timesteps`.

Run from this directory with:

```bash
./neko 1cyl_adjoint_consistency.case
```

The scalar result is written to `adjoint_consistency.dat`.
