# Passive scalar example (moving towards the mixer)
So in order to replicate the passive mixer problem by 
[C. S. Andreasen et al. 2009](https://doi.org/10.1002/fld.1964)
we need adjoint passive scalar capability.

Their objective function was essentially

min $\int_{\Gamma_{out}} (\phi - \bar{\phi})^2 d\Gamma$ 
(with some normalization)

s.t. $\Delta P \leq \beta \Delta P_{ref}$

They set $P_{out} = 0$ implying $\Delta P = \int_{\Gamma_{in}}p d \Gamma$

So there workflow involves:
- Solve the steady forward problem (velocity and passive scalar)
- Solve an adjoint problem (velocity and passive scalar) for the objective
- Solve an adjoint problem (only velocity) for the constraint

So for the passive scalar we need
- `adjoint_scalar_scheme`
- `adjoint_scalar_pnpn`
- update `adv_adjoint_dealias.f90`
- update `adv_adjoint_no_dealias.f90`
- `adjoint_scalar_convection_source_term`?

That last one I don't know a good name for, but it's the term arising from
linearizing the adjoint passive scalar convective term which then enters the
adjoint velocity equation. The term looks like $\nabla \phi \phi^\dagger$.

Then for this case specifically we need
- `enhanced_mixing_objective_function`
- `pressure_drop_constraint`
- `adjoint_enhanced_mixing_scalar_bc` (this enters the adjoint
passive scalar equation.)
- `adjoint_pressure_drop_BC` This enters the adjoint velocity equation to 
account for the pressure drop constraint.

The adjoint BCs are, at the time of writing, very poorly handled, so 
implementing those last two is going to be hard.

Often times one can replace an objective on a BC with a volume integral, 
ie, instead of of evaluated the "mixedness" on the outlet surface, measure
the mixedness in a volume downsteam close to the outlet.
This turns the surface integral (and subsequent adjoint BCs) into a 
volume integral (and subseqent adjoint source terms) which may be easier to
work with.

Hence, for testing, until we fix the handing of the adjoint BCs, I've 
implemented a `adjoint_enhanced_mixing_scalar_source_term` instead.

> I also feel the design needs to enter the passive scalar 
> equation. Probably not in the difusive term, but I think it should show up
> in the convective term. ie, not
>
> $(\mathbf{u} \cdot \nabla ) \phi$
>
> but $C(\rho) (\mathbf{u} \cdot \nabla ) \phi$ such that $C(\rho) \in [0,1]$
> 
> I could be wrong... I mean... I suppose if the Brinkman term is doing its
> job correctly we should have $\mathbf{u}= \mathbf{0}$ in the solid anyways..
> So maybe it's ok.


# Progress
OK so far NOTHING has been tested aside from compilation

We have basically cloned 
	`scalar_scheme.f90`
	`scalar_pnpn.f90`

Now we need to uncomment all the parts of `adjoint_simulation.f90` that 
correspond to the a passive scalar

This is something worth talking about actually.
Because in principal when you do an adjoint the order of everything reverses. 
Which means we'll go:
- step fluid forward
- step scalar forward
- step adjoint scalar backward
- step adjoint fluid backward

So this means the `adjoint_scalar_convection_source_term` is evaluated 
explicitly but on the correct timestep because we've stepped the scalar 
backwards before the adjoint fluid. So this is nice and I don't think
`Nek5000` did this (I could be wrong, I should double check)

I've added the adjoint convective term in 
	`sources/adjoint/adv_adjoint_dealias.f90`
	`sources/adjoint/adv_adjoint_no_dealias.f90`
But they still require testing.

TODO:
- obviously testing... but we don't really have a good case to test yet
When I did my derivation (I'll write it up neatly at some point...) 
- Double check with Casper what the BC's are, it looks like he keeps drawing
parabolic profiles, is there a periodic or symetric directions?

- The JSON stuff isn't finished
- The adjoint BCs aren't finished
- We need to look into non-dimensionalization etc more carefully, 
to make sure everything is consistent with the standard passive scalar. 
  (I think it is... but I would like to double check)
- boundary conditions (especially user) are not handled well for the adjoint!

