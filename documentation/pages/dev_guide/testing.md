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

## Regression testing {#testing-regression}

`tests/regression/mma/run.sh` and `tests/regression/sensitivity/run.sh` are
developer scripts that run their regression cases directly against an
already-built `reg_mma_bin` / `sensitivity_regression_driver`, independently
of `ctest`. Both locate the binary in the build tree (its
`RUNTIME_OUTPUT_DIRECTORY`) via `NEKO_TOP_BUILD_DIR`, which defaults to
`<repo root>/build` -- the directory `./setup.sh -t` creates -- and can be
set to a differently-named build directory.

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

## Keeping test executables in sync {#testing-build-fixtures}

CTest only runs tests; it does not build anything. After editing the sources,
build before testing, or `ctest` runs the executables from the previous build:

```bash
cmake --build build
ctest --test-dir build
```

Some regression tests register their build step as a test of its own, so that
running `ctest` alone keeps their executable in sync with the sources. The CI
job for the valgrind regression tests relies on this, as it runs `ctest`
without a separate build step. In `tests/regression/valgrind/CMakeLists.txt`:

- `regression_valgrind_build` runs `cmake --build` for `reg_valgrind_bin` and
  is the setup test (`FIXTURES_SETUP`) of the `regression_valgrind_bin`
  fixture.
- `regression_valgrind_steady_mixer` and `regression_valgrind_unsteady_mixer`
  require that fixture (`FIXTURES_REQUIRED`).

Concurrent `cmake --build` calls on one build tree race on shared outputs, so
every test that builds in the tree should take the `neko_top_build_tree`
resource lock (`RESOURCE_LOCK`). CTest then never runs two of them at once
under `ctest -j`.
