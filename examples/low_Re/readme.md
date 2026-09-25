# Low-Re POD {#low-re-pod}

This example is the POD state-recovery low-Re mixer case used on LUMI. It is
the cluster-oriented POD example in this split stack.

For manual testing outside Slurm, launch it from the repository root with
`./run.sh low_Re`, overriding `NEKO_RANKS` and `PY_RANKS` if needed.

For a LUMI submission, run `./run.sh --submit LUMI-G low_Re`. The
example-specific job script under `scripts/jobscripts/LUMI-G/low_Re` sets up
the full-node layout used by the coupled run:

- 8 Neko ranks per node, one GPU-backed rank per MI250x GCD
- 48 Python ranks per node, using the remaining CPU-only tasks

When building Neko on LUMI with CCE, disable Coarray Fortran in the compiler
flags. The coupled MPMD launch has a split Neko/Python communicator, and the
CAF runtime can otherwise trip Cray PMI/SMA during startup. Add `-hnocaf` to the
Fortran flags used by `prepare.env`, for example:

```bash
export NEKO_FCFLAGS="-O3 -m4 -hnocaf"
```

That job path generates `select_gpu` and `mpmd.conf` automatically before
calling `srun --multi-prog`, so the cluster-specific placement stays with the
example rather than in the generic MPMD helpers.
