# Passive scalar example (moving towards the mixer) {#passive-scalar}
So in order to replicate the passive mixer problem by 
[C. S. Andreasen et al. 2009](https://doi.org/10.1002/fld.1964)
we need adjoint passive scalar capability.

Their objective function was essentially

min \f$\int_{\Gamma_{out}} (\phi - \bar{\phi})^2 d\Gamma\f$ 
(with some normalization)

s.t. \f$\Delta P \leq \beta \Delta P_{ref}\f$

They set \f$P_{out} = 0\f$ implying \f$\Delta P = \int_{\Gamma_{in}}p d \Gamma\f$

So their workflow involves:
- Solve the steady forward problem (velocity and passive scalar)
- Solve an adjoint problem (velocity and passive scalar) for the objective
- Solve an adjoint problem (only velocity) for the constraint

So for the passive scalar we need
- `adjoint_scalar_scheme` :white_check_mark:
- `adjoint_scalar_pnpn` :white_check_mark:
- update `adv_adjoint_dealias.f90` :white_check_mark:
- update `adv_adjoint_no_dealias.f90` :white_check_mark:
- `adjoint_scalar_convection_source_term` :x:

That last one I don't know a good name for... but it's the term arising from
linearizing the adjoint passive scalar convective term which then enters the
adjoint velocity equation. The term looks like \f$\nabla \phi \phi^\dagger\f$.
The reason it's not fully completed with a :white_check_mark: is because I
beleive if one is running with `dealiased=true` then this term should be 
evaluated on the dealiased mesh, which I currently haven't implemented.

Then for this case specifically we need
- `enhanced_mixing_objective_function` :white_check_mark: 
- `pressure_drop_constraint` :x:
- `adjoint_mixing_scalar_source_term` (this enters the adjoint
passive scalar equation.) :white_check_mark:
- `adjoint_enhanced_mixing_scalar_bc` (this enters the adjoint
passive scalar equation.) :x:
- `adjoint_pressure_drop_BC` This enters the adjoint velocity equation to 
account for the pressure drop constraint. :x:
- `multi_objective_problem_type` :x: This will be a way of handling the
additional constraint and subsequent adjoint solve that comes with it.

The adjoint BCs are, at the time of writing, very poorly handled, so 
implementing those last two is going to be hard, and I'll save it for another
PR.

Often times one can replace an objective on a BC with a volume integral, 
ie, instead of of evaluated the "mixedness" on the outlet surface, measure
the mixedness in a volume downstream close to the outlet.
This turns the surface integral (and subsequent adjoint BCs) into a 
volume integral (and subsequent adjoint source terms) which may be easier to
work with.

Hence, for testing, until we fix the handing of the adjoint BCs, I've 
implemented a `adjoint_enhanced_mixing_scalar_source_term` instead.

> I also feel the design needs to enter the passive scalar
> equation. Probably not in the diffusive term, but I think it should show up
> in the convective term. ie, not
>
> \f$(\mathbf{u} \cdot \nabla ) \phi\f$
>
> but \f$C(\rho) (\mathbf{u} \cdot \nabla ) \phi\f$ such that \f$C(\rho) \in [0,1]\f$
>
> I could be wrong... I mean... I suppose if the Brinkman term is doing its
> job correctly we should have \f$\mathbf{u}= \mathbf{0}\f$ in the solid anyways..
> So maybe it's ok.

# Progress

## Solver
We have basically cloned 
- `scalar_scheme.f90`
- `scalar_pnpn.f90`

and adjust the BCs such that `'v' -> 'w'` and changed the order,
because in principal when you do an adjoint the order of everything reverses. 
Which means we'll go:
- step fluid forward
- step scalar forward
- step adjoint scalar backward
- step adjoint fluid backward

So this means the `adjoint_scalar_convection_source_term` is evaluated 
explicitly but on the correct timestep because we've stepped the scalar 
backwards before the adjoint fluid. So this is nice.

A little caveat worth mentioning: Due to the one way coupling, I believe one 
could do this a bit more efficiently in a steady context by converging only
the fluid to a steady state and only then start the passive scalar, etc for
the adjoint (instead of having the fluid and the scalar solved together at
each timestep)
However, for unsteady simulations you would want to solve it the way it's being
done right now.

## Convective term
I've added the adjoint convective term in 
- `sources/adjoint/adv_adjoint_dealias.f90`
- `sources/adjoint/adv_adjoint_no_dealias.f90`
- `sources/adjoint/adjoint_scalar_convection_source_term.f90`

## Optimization
Just for this case I've implemented
- `scalar_mixing_objective_function.f90`
- `adjoint_mixing_scalar_source_term.f90`

to test out an unconstrained version of the mixing problem.

There's also been a change to `objective_function_t` such that the inputs are
primal and adjoint `case_t` as opposed to `fluid_scheme_t`. This means we have
access to both the fluid and the scalar.

## general todo's:
- obviously testing... but we don't really have a good case to test yet
When I did my derivation (I'll write it up neatly at some point...) 
- Double check with Casper what the BC's are, it looks like he keeps drawing
parabolic profiles, is there a periodic or symmetric directions?
- The JSON stuff isn't finished
- The adjoint BCs aren't finished
- We need to look into non-dimensionalization etc more carefully, 
to make sure everything is consistent with the standard passive scalar. 
  (I think it is... but I would like to double check)
- boundary conditions (especially user) are not handled well for the adjoint!

