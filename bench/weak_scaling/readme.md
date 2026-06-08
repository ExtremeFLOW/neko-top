# Weak scaling benchmark for static mixer

This directory contains the weak-scaling benchmark workflow for the static mixer
optimization example.

## Example executed

The benchmark cases are generated from the static mixer setup through
bench/weak_scaling/create_examples.sh.

- Mesh pattern: mixer.
- Base data path: data_local/static_mixer.
- Benchmark case path: examples/benchmark/<tag>/.
- Experiment metadata path: results/benchmark/experiments/<tag>/.

For weak scaling, the scripted cases keep work per GPU approximately constant
while increasing global mesh size and GPU count. On LUMI-G, the case generator
uses 8 GPUs per node, so the mesh-to-node mapping depends on the selected
checkpoint-memory setting (n_memory).

The weak-scaling suites are:

| Suite | n_memory | Elements per GPU |
| ----- | -------- | ---------------- |
| A     | 500      | 8192             |
| B     | 250      | 16384            |

Mesh progression by node count on LUMI-G:

| Nodes | GPUs | Mesh (Suite A: n_memory=500) | Mesh (Suite B: n_memory=250) |
| ----- | ---- | ---------------------------- | ---------------------------- |
| 1     | 8    | 64x32x32                     | 128x32x32                    |
| 2     | 16   | 128x32x32                    | 256x32x32                    |
| 4     | 32   | 256x32x32                    | 128x64x64                    |
| 8     | 64   | 128x64x64                    | 256x64x64                    |
| 16    | 128  | 256x64x64                    | 512x64x64                    |

Each generated case is tracked by a tag, which is typically aligned with a git
history step used in performance comparisons.

## Metrics analyzed in the notebook

The analysis notebook in bench/weak_scaling/analyse.ipynb processes experiment
logs and runtime statistics and focuses on:

- Time to solution.
- Time-Step.
- Time-Step Adjoint.
- Checkpoint save.
- Checkpoint restore.
- Optimizer iteration.
- MMA update (with alias support for MMA gensub).
- MMA KKT computation (with alias support for MMA subsolve).

From these measurements, the notebook computes per-element timing summaries and
weak-scaling efficiency curves (including ideal-efficiency references), and it
compares efficiency evolution across tags.

## Benchmark cases constructed

The executed case order is defined by bench/weak_scaling/order.txt and must be
preserved. The cases are grouped by theme below while keeping strict ordering.

### Solver and threading progression

1. baseline
	- Baseline benchmark tag used as reference for later comparisons.
	- No direct commit title match in current local history.
2. fused_cg
	- Solver-focused update represented by fused_cg tag.
	- Git history match: bb9b8781 (Update to fused_cg).
3. OpenMP
	- OpenMP-enabled variant for threaded execution changes.
	- Git history match: c53b374a (OpenMP).

### Memory and checkpointing progression

4. LimitRAMCheckpoints
	- Memory-limiting checkpoint variant used for RAM footprint studies.
	- No direct commit title match in current local history.
5. trento_module
	- trento_module tagged implementation variant.
	- Git history match: 7c1354b2 (trento_module).
6. load_balancing
	- Load-balancing benchmark variant.
	- No direct commit title match in current local history.
7. phmg
	- PHMG-tagged variant used in the benchmark sequence.
	- Git history match: a1bcb44f (phmg).
8. hdf5_checkpoint
	- HDF5 checkpointing variant.
	- Git history match: 42a80911 (hdf5_checkpoint).

### IO striping and output format progression

9. stripe_by_tasks
	- Stripe-by-task IO layout variant.
	- Git history match: c8b032a5 (stripe_by_tasks).
10. stripe_minus_one
	- Alternative striping strategy variant.
	- Git history match: 31fe759f (stripe_minus_one).
11. chunk_size_2M
	- Chunk-size tuning variant (2M).
	- Git history match: cbbbf7d4 (chunk_size_2M).
12. chunk_size_8M
	- Chunk-size tuning variant (8M).
	- Git history match: 08eadee9 (chunk_size_8M).
13. vtkhdf
	- VTKHDF output variant.
	- Git history match: 8c15d988 (vtkhdf).

Cases without a direct commit-title match above are still part of the executed
benchmark order and are kept as explicit benchmark tags.