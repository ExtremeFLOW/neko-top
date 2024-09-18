OK so far we NOTHING has been tested aside from compilation

We have basically cloned 
	scalar_scheme.f90
	scalar_pnpn.f90

I've added most of the adjoint convective term in 
	sources/adjoint/adv_adjoint_dealias.f90
	sources/adjoint/adv_adjoint_no_dealias.f90
But they still require testing.

TODO:
	- obviously testing... but we don't really have a good case to test yet
	- When I did my derivation (I'll write it up neatly at some point...) we also get a source term in velocity equation of the form:
	$\nabla s s_adj$
	So I'll add another source term type to handle this term, and I guess we just treat it explicitly to avoid too much coupling.

	This is something worth talking about actually.
	Because in principal when you do an adjoint the order of everything reverses. Which means we'll go:
	step fluid forward
	step scalar forward
	step adjoint scalar backward
	step adjoint fluid backward

	So this means it's not explicit in the same way as the brinkman term, where it's evaulated on the previous timestep.
	It's explicit and evaluated on the correct timestep because we've stepped the scalar backwards before the adjoint fluid.
	So this is nice!!!

	- The JSON stuff isn't finished
	- The adjoint BCs aren't finished

