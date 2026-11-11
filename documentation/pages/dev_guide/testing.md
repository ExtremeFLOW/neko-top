# Testing {#testing}

\note This section is under construction.


\warning The use of pFUnit requires an up to date version as there are issues in
the pFUnit CTest macros that are fixed in the latest version. In particular the
fix was included in version 4.7.1 of pFUnit.


## Unit Testing

## Build types {#testing-build-types}

`./setup.sh -t` builds tests and, unless a build type is already chosen
(the `CMAKE_BUILD_TYPE` environment variable), defaults `CMAKE_BUILD_TYPE` to
`Testing`. `Testing` is an optimised configuration (`-g -O2` for GNU, debug
symbols kept) with the same warning set as `Debug`, so the test suite runs at
a realistic speed while still surfacing compiler warnings. See the
[setup script variables](@ref installation-setup) for the other build types
and `CMAKE_BUILD_TYPE`.

## Runtime checks {#testing-runtime-checks}

GNU builds can additionally enable bounds/do/pointer/recursion checking,
floating-point exception traps and poisoned values for uninitialised variables
via `-DNEKO_TOP_RUNTIME_CHECKS=ON` (default `AUTO`). Setting `AUTO` enables
runtime known passing checks for `Debug` builds and disables them otherwise.
Setting `ON` enables all runtime checks regardless of the build type, and
setting `OFF` disables all runtime checks.

### Known failing tests with runtime checks enabled

- Several tests currently fail with it enabled, due to collapsed-index loops in
the sources (for example `x(i,1,1,1)` with `i` running over the whole array,
which `-fcheck=bounds` rejects).

## Validation tests {#testing-validation}

Validation tests reproduce published results and run for hours, so they are
not registered by default. Configure with `-DNEKO_TOP_VALIDATION_TESTS=ON` to
register them, and select them with `ctest -L validation`. They do not carry
the `regression` label, so `ctest -L regression` never runs them.

- `borrvall_rugby_ball_validation_0p8` runs the Borrvall & Petersson (2003)
  rugby ball with a fluid volume fraction of 0.8 (about 80 minutes on eight
  MPI ranks; set `NEKO_TOP_RUGBY_BALL_RANKS` to use another count). It checks
  that the final objective is within 2% of the paper's 31.75, that the final
  design is feasible, and that the paper's stopping rule is met before the
  last iteration.

The fast `borrvall_rugby_ball_regression` test runs the same physics on a
small mesh for a few iterations and compares the whole optimization history
with a stored reference. It is part of the default suite.
