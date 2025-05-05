# Using `neko-top` and `LightKrylov` to perform stability analysis
This example demonstrates how the modulararity of both `neko-top` and
[`LightKrylov`](https://github.com/nekStab/LightKrylov) can be utilized to
perform linear stability analysis, attempting to replicate some of the
functionality of either the
[KTH_Framework](https://github.com/KTH-Nek5000/KTH_Framework)
or
[nekStab](https://github.com/nekStab/nekStab) projects with the addition of
GPU acceleration.

## guide for `LightKrylov` and `neko-top` people
`neko-top` is a framework intended to use `neko` as a library to perform topology
optimization. Currently, the linear and adjoint solvers only exist in `neko-top`
and are yet to be migrated to `neko`, hence, this example being run through the
`neko-top` framework. The `setup.sh` script should handle all the external
dependencies related to `neko` compiling and installing them in `external`. These
include
- `neko`
- `json-fortran`
- `gslib`
- `hdf5`

and relies on
- `CMake`

`LightKrylov` requires the additional dependency
- `fpm`

Finally, I used `gmsh` to generate the mesh, so that's one more dependency.

hopefully :pray:, if all goes well, one can run
```bash
git clone --recursive https://github.com/ExtremeFlow/Neko-TOP.git neko-top
cd neko-top
git checkout LightKrylov/linking_dependencies
./setup.sh
./run.sh ext_cyl
```



## Issues on the `neko` / `neko-top` side
The following issues should be addressed before this can be used for production.

- `neko` currently doesn't have a true 2D solver, so there were a lot of hacky
additions that were required to force this flow 2D.
- The linear/adjoint and even non-linear solvers currently still rely heavily
on the field registry, so spawning multiple solvers is still cumbersome.
- Currently this only works for first order time integration, implying something
has gone wrong with the way I'm resetting the adjoint solver between runs.
- A sponge was implemented (in a very hardcoded way) but this should be extended
to a `neko` source term that utilizes the `.case` file properly.
- We currently don't have SFD or `bostconv`, so we can't find our own baseflows
using `neko`.
- If we want to this available to be used by the `LightKrylov` team we need to
update the linear/adjoint solver and migrate them from `neko-top` to `neko`. We
then need to find some kind of `fpm` compatibility.

## Issues on the `LightKrylov` side
There are a few "issues" with `LightKrylov` that have forced me to fork it as
opposed to using the `dev` branch along.

This is primarily with the frequent use of `intent(out)` for operations using
`abstract_vector` types. For the `field_t` in `neko` one must first initialize
using a `dofmap_t`, implying the `abstract_vector` we would use would ideally
have an `init` functionality which could allocate the required fields. Indeed,
before Krylov subspace is allocated, it is followed by a `%zero()`, allowing
us to sneak in an initialization during this call. However, the use of
`intent(out)` in certain subroutines will deallocate the `abstract_vector` upon
entry to the subroutine, requiring it to be re-initialized with the subroutine.

On top of this, when allocating in `LightKrylov`, for example in `linear_combination` 
it is common to use 
```fortran
allocate(Y(size(B, 2)), source=X(1))
```

Given fields are pointers in the `abstract_vector`, this `=` assignment passes
the pointers across, implying we edit the field in `X(1)`. There was an attempt
to reassign `=` with

```fortran
   interface assignment(=)
    module procedure state_vector_assignment
   end interface
```

But it still took the abstract class' assignment, not the derived type.

Instead A VERY hacky solution of a deep copy was implemented in `abpby`.

I would recommend adjusting the `abstract_vector` to do a deep copy on `=`.

***Edit*** This was overcome by first free'ing and then initializing every time
`%zero()` was called to avoid a shallow copy. But now I'm suspicious that I'm
never actually free'ing any fields, just nullifying their pointers, so there
is potentially a massive memory leak in here...