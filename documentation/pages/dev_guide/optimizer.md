# Optimizer {#optimizer}

\tableofcontents

## Overview {#optimizer-overview}

An optimizer drives the outer iteration of a topology-optimization run: given
a \ref problem (objective, constraints and their sensitivities) and a
\ref design "design" (the variables being optimized), it repeatedly updates
the design until one of its stopping criteria is met. The abstract base type
`optimizer_t` (`sources/optimizer/optimizer.f90`) implements the parts common
to every optimizer — the run loop, checkpointing, the convergence log — and
leaves the update rule itself to a concrete extension.

A run selects its optimizer through `optimization.solver.type` in the case
file; `optimizer_factory` (`sources/optimizer/optimizer_factory.f90`)
allocates and initializes the matching concrete type. See
\ref configuration-solver "Optimization Solver settings" for the case-file
side of this.

## The `optimizer_t` interface {#optimizer-interface}

A concrete optimizer extends `optimizer_t` and implements its deferred
procedures:

- `init_from_json(parameters, problem, design, simulation)` — reads
  `optimization.solver` (and the optimizer's own nested object within it)
  from the case file, validates it, and prepares any internal state. Must
  call `this%init_base(...)` (see below).
- `free()` — releases the optimizer's own state, then calls
  `this%free_base()`.
- `initialize(problem, design, simulation)` — evaluates the problem (and its
  sensitivity) at the starting design, before the first iteration.
- `step(iter, problem, design, simulation)` — performs one optimization
  iteration: computes an update from the current sensitivity, applies it to
  `design` (`design%update_design`), and re-evaluates the problem at the new
  design. Returns `.true.` only if the optimizer has itself converged (the
  meaning of "converged" is specific to the method); returning `.false.`
  leaves convergence to the base run loop's own criteria.
  \warning `step` **must** set `this%max_design_change` and
  `this%norm2_design_change` on every call, including any early-return path.
  The base run loop reads both immediately after `step` returns — to decide
  the `stop_design_change` criterion and to log them — so an iteration that
  does not set them stalls that criterion or logs stale values.
- `validate(problem, design)` — checked once, after the loop ends, against
  the final design.
- `write(iter, problem)` — appends the optimizer's convergence log entry
  (`optimization_data.csv`) and any optimizer-specific output, normally via
  `this%write_log(iter, problem, extra_values)`.
- `save_checkpoint_components(filename, overwrite)` /
  `load_checkpoint_components(filename)` — save/restore whatever internal
  state the optimizer needs to resume correctly (e.g. history vectors). The
  base `optimizer_t` state and the design are handled separately, by
  `optimizer_t` and `design_t` themselves.

`init_from_json` should call, in this order: `this%free()` (which resets the
base settings read by `read_base_settings` — notably `stop_design_change` and
the checkpoint fields — to their un-set defaults); `this%read_base_settings
(solver_parameters)`; the optimizer's own validated settings; and finally
`this%init_base(optimizer_type, max_iterations, ...)`, which stores the
mandatory `optimizer_type`/`max_iterations` and any of the optional base
settings the concrete optimizer passes through. Reversing this order silently
loses settings, because `free`/`init_base` reset or overwrite fields that
`read_base_settings` has already populated.

## The base run loop {#optimizer-run-loop}

`optimizer_t%run` (`optimizer_run` in `optimizer.f90`) implements the loop
shared by every optimizer:

1. **Restart.** If an explicit `restart_file` was given, or a runtime
   checkpoint (`optimizer_rt_checkpoint.<ext>`) exists at the default
   checkpoint location, `load_checkpoint` restores the base state, calls
   `load_checkpoint_components`, and restores the design — setting
   `this%current_iteration` to the checkpointed iteration.
2. **Continuation restore.** If restarting (`current_iteration > 0`), the
   design's and simulation's output counters are set to
   `current_iteration - 1`, and the \ref continuation "continuation"
   scheduler is updated to `current_iteration`, so the problem below is
   re-evaluated with the same continuation parameter values as the
   checkpointed iteration. A fresh start skips this (`update(0)` would index
   the scheduler's value array out of range).
3. **Initialize.** `simulation%set_design_iteration` (if a simulation is
   present) and then `this%initialize(problem, design, simulation)` evaluate
   the problem and its sensitivity at the starting design.
4. **Initial log entry.** On a fresh start only, `this%write(0, problem)` logs
   the starting point; `design%write(current_iteration)` always runs.
5. **Per-iteration loop**, while `current_iteration < max_iterations`:
   increment `current_iteration`; update the continuation scheduler and
   `simulation%set_design_iteration` for the new iteration; call
   `this%step(current_iteration, problem, design, simulation)`; log with
   `this%write` and `design%write`; save a checkpoint if
   `checkpoint.interval` divides the iteration count.
6. **Stop checks**, evaluated in this order after each `step`: `step`
   reported convergence (stop reason 0); `current_iteration >=
   max_iterations` (reason 1); `max_design_change < stop_design_change`
   (reason 2); out of `max_runtime` (reason 3, which also saves a runtime
   checkpoint before exiting).
7. **Finish.** `this%validate(problem, design)` runs once the loop exits, then
   the stop reason is logged.

## Concrete optimizers {#optimizer-concrete}

### MMA {#optimizer-concrete-mma}

`mma_optimizer_t` (`sources/optimizer/mma_optimizer.f90`) implements the
Method of Moving Asymptotes, and handles both constrained and unconstrained
problems. See \subpage MMA for the algorithm and
\ref configuration-solver-mma "MMA Options" for its case-file settings.

### Steepest descent with move limits {#optimizer-concrete-steepest-descent}

`steepest_descent_optimizer_t`
(`sources/optimizer/steepest_descent_optimizer.f90`) implements the steepest
descent method with move limits (SDM), Appendix D of da Silva, Aage, Beck &
Sigmund (2021). It handles **unconstrained** problems only — `init_from_json`
raises an error pointing to `mma` if the problem has constraints. See
\ref configuration-solver-steepest-descent "Steepest descent options" for the
algorithm and its case-file settings.

Its elementwise update rules live in `steepest_descent_math.f90`, which
dispatches each one (by `NEKO_BCKND_DEVICE`) to the CPU implementation in
`bcknd/cpu/steepest_descent_cpu.f90` or to the device one in
`bcknd/device_steepest_descent.f90` — this is the same CPU/device mirroring
pattern described below for a new optimizer's own kernels.

## Adding a new optimizer {#optimizer-adding}

1. Define a type `extends(optimizer_t)` in its own module under
   `sources/optimizer/`, implementing every deferred procedure above (see
   `steepest_descent_optimizer.f90` for a complete, currently-unconstrained
   example, or `mma_optimizer.f90` for one that also supports constraints).
2. Register it in the factory submodule,
   `sources/optimizer/optimizer_factory.f90`: add its case-file type string
   to the `KNOWN_TYPES` array (used in the error message when `type` is
   unrecognized) and add a matching `case` branch in `optimizer_factory` that
   allocates the new type before the shared `call
   object%init_from_json(...)`.
3. Add its sources to `sources/optimizer/CMakeLists.txt`. Any CUDA/HIP
   sources and their include directories go behind `if(DEVICE_TYPE STREQUAL
   "CUDA")` / `elseif(DEVICE_TYPE STREQUAL "HIP")` guards, as the steepest
   descent kernels already do there.
4. If the optimizer needs state beyond what `optimizer_t` already
   checkpoints, add an HDF5 IO submodule,
   `sources/IO/hdf5/<name>_hdf5_io.f90`, as `submodule (<name>_optimizer)
   <name>_hdf5_io` implementing `save_checkpoint_components`/
   `load_checkpoint_components`'s `module subroutine` bodies (declared as
   `interface` blocks, plus `#if !HAVE_HDF5` stubs that call `neko_error`, in
   the parent module — see `optimizer.f90`'s own
   `optimizer_save/load_checkpoint_hdf5` for the pattern, and
   `steepest_descent_hdf5_io.f90` for a full worked example, including the
   collective, hyperslab-based read/write of a distributed `vector_t`). Wire
   it into `sources/IO/CMakeLists.txt` behind
   `$<$<BOOL:${HAVE_HDF5}>:...>`.
5. If the update rule needs device-accelerated kernels, mirror
   `sources/optimizer/bcknd/{cpu,device/{cuda,hip}}/`: CPU routines operating
   on plain arrays, CUDA/HIP kernels operating on device pointers, and one
   dispatch module `bcknd/device_<name>.f90` selecting the `HAVE_HIP`/
   `HAVE_CUDA` `bind(c)` interface at compile time — follow
   `device_steepest_descent.f90` (or
   `mapping_functions/bcknd/device_SIMP_mapping.f90` for an older instance of
   the same pattern).
6. Document the new case-file keys in `user_guide/configuration.md`, in the
   same table format as the existing optimizers, and give every new public
   procedure and keyword argument its own `@param` — this page and those
   docstrings are what a user or reviewer is pointed back to.
