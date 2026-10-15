# Testing {#testing}

\note This section is under construction.


\warning The use of pFUnit requires an up to date version as there are issues in
the pFUnit CTest macros that are fixed in the latest version. In particular the
fix was included in version 4.7.1 of pFUnit.


## Unit Testing

## Build types {#testing-build-types}

`./setup.sh -t` builds tests and, unless a build type is already chosen
(either `-DCMAKE_BUILD_TYPE=...` in `CMAKE_VARIABLES`, or the
`CMAKE_BUILD_TYPE` environment variable), defaults `CMAKE_BUILD_TYPE` to
`Testing`. `Testing` is an optimised configuration (`-g -O2` for GNU, debug
symbols kept) with the same warning set as `Debug`, so the test suite runs at
a realistic speed while still surfacing compiler warnings. See the
[setup script variables](@ref installation-setup) for the other build types
and `CMAKE_BUILD_TYPE`.

## Runtime checks {#testing-runtime-checks}

GNU builds can additionally enable bounds/do/pointer/recursion checking,
floating-point exception traps and poisoned values for uninitialised
variables via `-DNEKO_TOP_RUNTIME_CHECKS=ON` (default `OFF`;
`AUTO` is reserved for future use and currently behaves as `OFF`). This is
experimental: several tests currently fail with it enabled, due to
collapsed-index loops in the sources (for example `x(i,1,1,1)` with `i`
running over the whole array, which `-fcheck=bounds` rejects). Only use it to
investigate a specific suspected memory or FPE issue, not as part of a normal
test run.

