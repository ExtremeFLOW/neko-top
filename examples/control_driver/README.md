# Proof of concept for a steady state optimization problem

Here we follow the classic "rugby ball" example by [Borrvall & Petersson 2003](https://doi.org/10.1002/fld.426)



This is example is primarily used to summarize the list of remaining tasks.

However, the source code is littered with TODO comments and further details.


## Todo (in Neko):


### Simcomp_handler 
this would be like what Tim made for source terms. I feel like Tim is the man for this job.
Importantly, simulation components are initialized in `neko.f90` by passing the case through. It would be nice to have functionality to append simcomps, this way we can pass `C` or `adj`.

The reason we need this is because currently computing an objective function or sensitivity is done with a proceedure owned by the `objective_function_t`, and is computed on the last step. When we move to unsteady calculations these will be accumulated during the run, so a `simcomp_t` will be the perfect tool. So then we can start appending simcomps that compute time averaged objective functions and sensitivity.

### Steady states
I would preffer if this is done in Neko, because the simulation component we have right now is a bit hacky. I suppose have a `am_I_done` function that returns a boolean of y/n is a good idea!
- Then in most cases this is just checking `t < t_final`
- But then we could check norms etc
- Other people could check converged statistics etc
- could even have a usr defined stopping condition too.


## Todo (in Neko-top):

### Write the `problem` module
You can see in `driver.f90` I've put a fair amount of the skeleton for the `problem` module, which would hold everything. But really this driver should only be calling:
```fortran
	problem%init()
	problem%solve()
	problem%free()
```
But we should have a discussion about how we want to lay out the case file first for a neko-top case.

Honestly, this bit ***really*** confuses me. Because in theory, a `problem` can exist without a fluid, and it's only when we start defining certain objective functions and constraints that we require a fluid simulation. This level of abstraction/objected oriented thinking is beyond me, but I'm sure we'll require some refactoring to make this right.

### Masks
This one is rather simple, just including masks on objective functions and masks on the optimization domain. We could then compute volume based on the optimization domain, not the computation domain.

We should discuss though, how we define a mask. Could be point zones like you suggested Tim. In `Nek5000` we filled up arrays the size of the mesh with booleans. That's rather light, and we could still use all the point zone functionality to initialize these masks.

### incorporate the PDE filters
We already have a branch in Neko where we did the PDE filters. Either we PR that in neko, so it can be used with the standard Brinkman term there, or we just migrate everything to neko-top. Then we can start including it in the `design_t`.

### Handle mapping to multiple coefficients 
This is a BIG one, and we need to discuss the best way forward before implementing anything. 

Right now, a `design_t` only has the Brinkman term to map to. But in principle, we may want to map to many coefficients, eg conductivity $\kappa$ etc for CHT. So I'm thinking it would be better for us to have a `field_list_t` for all the coefficients, and then each one has a unique subroutine for its mapping. 

So when we call `design%map_forward()`, it loops through and maps all the coefficients. I feel like this is going to be tricky to keep track of, but we could use a `get_by_name` so that's nice. It's like our own little self contained registry of sorts.

Similarly, the `objective_function_t` type has a field to hold the sensitivity wrt to the coefficient (for us, it's just the brinkman term $\chi$) ie, $\frac{\partial F}{\partial \chi}$. This gets passed to the `design_t` to chain rule backwards, ie, $\frac{\partial F}{\partial \rho}$. 

But with more coefficients we'll need to pass a list back, ie $\frac{\partial F}{\partial \chi}$, $\frac{\partial F}{\partial C}$, $\frac{\partial F}{\partial \kappa}$ etc, and some kind of instructions of how to assemble the sensitivity. 

All in all, this is a bit tricky in my mind.

### Handle multiple constraints (such that we have multiple sensitivities)
In a similar vibe, we would want a design to have a `field_list_t` for the sensitivities coming from the objective function $F$ and all the constraints $C_i$. ie, $\frac{\partial F}{\partial \rho}$, $\frac{\partial C_1}{\partial \rho}$, $\frac{\partial C_2}{\partial \rho}$  etc. 

All the initialization and execution of all of these need to be wrapped up nicely, similar to `source_terms` I guess.


### log file to write with KKT, constraint values etc
Of course we would still have the neko log file, but we should have a neko-top log file that prints out all the usefull information at each optimization iteration.

### merge in the passive scalars
We have most of the adjoint passive scalars written on a different branch, but we haven't been able to test it since we don't have a driver.


### ALL the JSON stuff
You can see everything is hardcoded and we need to include a bunch of JSON stuff. This can be done after we discuss the neko-top case file.

## Todo (in general):
- There are a ton more minor comments written throughout the source code here, so I'm sure some have been missed in this `readme.md`
- Double check all non-dimensionalization and scalings
- Implement the `adjoint_minimum_dissipation_source_term_t` correctly (ie, in weak form)
- Restart correctly (or maybe not! in steady calculations it can often be faster if your initial condition is the final condition from the previous iteration)
- Nothing is freed correctly
- Particularly regarding steady calculations, there's no reason to recalculate the adjoint forcing each time, it will always be the same.
 
