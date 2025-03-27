# Proof of concept for a steady state optimization problem {#rugby-ball}

Here we follow the classic "rugby ball" example by 
[Borrvall & Petersson 2003](https://doi.org/10.1002/fld.426)

Here we introduce the following new classes:


## `design_variable`
This could be any design variable in the optimization problem. eg :
- coefficients
- splines
- levelsets
- etc

Currently,
we only consider topology optimization and hence we just a type `topopt_design`.

Hopefully, this can be a derived type of the more abstract `design_variable`.

### `topopt_design`
This contains the design field \f$\rho\f$ as well as the following key proceedures
- `map_forward` maps the design varaible to coeffients (currently only the 
Brinkman amplitude \f$\rho \mapsto \chi \f$)

In the future, this should also contain filters \f$\rho \mapsto \tilde{\rho}
\mapsto \chi\f$

- `map_backwards` uses chain rule to map the sensitivity of the coefficents
to the sensitivity to the design variable \f$\frac{\partial F}{\partial \chi}
 \mapsto \frac{\partial F}{\partial \rho}\f$

## Optimizer

A generic type, optimizer, has been added. Currently, only the `mma_optimizer` 
is supported, and all required parameters for the MMA solver are specified in 
the JSON case file.

### New Parameters
1. **Scaling**  
   Based on the [Svanberg MMA Method](https://people.kth.se/~krille/mmagcmma.pdf), MMA performs optimally
   when function values are between 1 and 100. To achieve this:
   - If `auto_scale` is set to `true` in the case file:
     - The objective value is scaled at each iteration to 
	   match the specified `scale` parameter value. The derivatives are also 
	   scaled accordingly. (variable scaling factor)
   - If `auto_scale` is set to `false`:
     - The objective value and its derivatives are multiplied by the `scale` 
	   parameter value before being passed to the optimizer. (constant scaling
	   factor)
2. **Impact of Scaling**  
   Note that these scaling operations affect only the optimizer's internal 
   computations. The original values of the objective function and its 
   derivatives remain unchanged.


## `problem`
This contains a way to evaluate an optimization problem, ie, all objective
functions and constraints. Currently, we only
use gradient descent algorithms using first order gradient information, so
the key proceedures are
- `init` to initialize the type of optimization problem
- `compute` this is to evalute all objective functions \f$F\f$ and constraints \f$C_i\f$ 
- `compute_sensitivity` this is to evalute all the sensitivities
\f$\frac{\partial F}{\partial \rho}\f$, \f$\frac{\partial C_1}{\partial \rho}\f$,
\f$\frac{\partial C_2}{\partial \rho}\f$ etc.
- in princple only the `compute` is required. (eg, gradient free algorithms) or
more may be required (eg, `compute_hessian`)

Within a `problem` we also contain fluid and adjoint schemes as well as instance
 of the following classes.

Currently we have only implemented a steady topologogy optimization case.



### `objective_function`
this could be either an objective function or a constraint.
The key proceedures are
- `init` allowing various source terms to be appended to the adjoint.
- `compute` to compute the objective/constraint function value.
- `compute_sensitivity` to compute the sensitivity of the objective/constraint
with respect to the design.
- `free` to free.

Here we have made 2 derived types,

`minimum_dissipation_objective function` for objectives
\f$F = \int_\Omega |\nabla \mathbf{u}|^2 d\Omega + K \int_\Omega \chi |u|^2 d\Omega\f$

`volume_constraint` for constraints either \f$V < V_{max}\f$ or \f$V > V_{min}\f$ where
\f$V = \int_\Omega \tilde{\rho} d \Omega \f$.

Note, it should be over \f$\Omega_O\f$ after masks are included.

This is example is primarily used to summarize the list of remaining tasks.

However, the source code is littered with TODO comments and further details.


## Todo (in Neko):


### Simcomp_handler 
this would be like what Tim made for source terms. I feel like Tim is the man 
for this job. Importantly, simulation components are initialized in `neko.f90` 
by passing the case through. It would be nice to have functionality to append 
simcomps, this way we can pass `C` or `adj`.

The reason we need this is because currently computing an objective function or 
sensitivity is done with a proceedure owned by the `objective_function_t`, and 
is computed on the last step. When we move to unsteady calculations these will 
be accumulated during the run, so a `simcomp_t` will be the perfect tool. 
So then we can start appending simcomps that compute time averaged objective 
functions and sensitivity.

### Steady states
I would preffer if this is done in Neko, because the simulation component 
we have right now is a bit hacky. I suppose have a `am_I_done` function that 
returns a boolean of y/n is a good idea!
- Then in most cases this is just checking `t < t_final`
- But then we could check norms etc
- Other people could check converged statistics etc
- could even have a usr defined stopping condition too.


## Todo (in Neko-top):

### Write the `problem` module
You can see in `driver.f90` I've put a fair amount of the skeleton for the 
`problem` module, which would hold everything. 
But really this driver should only be calling:
```fortran
	problem%init()
	problem%solve()
	problem%free()
```
But we should have a discussion about how we want to lay out the case file 
first for a neko-top case.

Honestly, this bit ***really*** confuses me. Because in theory, a `problem` 
can exist without a fluid, and it's only when we start defining certain 
objective functions and constraints that we require a fluid simulation. 
This level of abstraction/objected oriented thinking is beyond me, but I'm 
sure we'll require some refactoring to make this right.

### Masks
This one is rather simple, just including masks on objective functions and 
masks on the optimization domain. We could then compute volume based on the 
optimization domain, not the computation domain.

We should discuss though, how we define a mask. Could be point zones like you 
suggested Tim. In `Nek5000` we filled up arrays the size of the mesh with 
booleans. That's rather light, and we could still use all the point zone 
functionality to initialize these masks.

### incorporate the PDE filters
We already have a branch in Neko where we did the PDE filters. Either we PR 
that in neko, so it can be used with the standard Brinkman term there, or we 
just migrate everything to neko-top. 
Then we can start including it in the `design_t`.

### Handle mapping to multiple coefficients 
This is a BIG one, and we need to discuss the best way forward before 
implementing anything. 

Right now, a `design_t` only has the Brinkman term to map to. But in principle, 
we may want to map to many coefficients, eg conductivity \f$\kappa\f$ etc for CHT. 
So I'm thinking it would be better for us to have a `field_list_t` for all the 
coefficients, and then each one has a unique subroutine for its mapping. 

So when we call `design%map_forward()`, it loops through and maps all the 
coefficients. I feel like this is going to be tricky to keep track of, but we 
could use a `get_by_name` so that's nice. It's like our own little self 
contained registry of sorts.

Similarly, the `objective_function_t` type has a field to hold the sensitivity 
wrt to the coefficient (for us, it's just the brinkman term \f$\chi\f$) 
ie, \f$\frac{\partial F}{\partial \chi}\f$. This gets passed to the `design_t` 
to chain rule backwards, ie, \f$\frac{\partial F}{\partial \rho}\f$. 

But with more coefficients we'll need to pass a list back, 
ie \f$\frac{\partial F}{\partial \chi}\f$, \f$\frac{\partial F}{\partial C}\f$, 
\f$\frac{\partial F}{\partial \kappa}\f$ etc, 
and some kind of instructions of how to assemble the sensitivity. 

All in all, this is a bit tricky in my mind.

### Handle multiple constraints (such that we have multiple sensitivities)
In a similar vibe, we would want a design to have a `field_list_t` for the 
sensitivities coming from the objective function \f$F\f$ and all the constraints 
\f$C_i\f$. ie, \f$\frac{\partial F}{\partial \rho}\f$, 
\f$\frac{\partial C_1}{\partial \rho}\f$, \f$\frac{\partial C_2}{\partial \rho}\f$  etc. 

All the initialization and execution of all of these need to be wrapped up 
nicely, similar to `source_terms` I guess.


### log file to write with KKT, constraint values etc
Of course we would still have the neko log file, but we should have a neko-top 
log file that prints out all the usefull information at each optimization 
iteration.

### merge in the passive scalars
We have most of the adjoint passive scalars written on a different branch, 
but we haven't been able to test it since we don't have a driver.


### ALL the JSON stuff
You can see everything is hardcoded and we need to include a bunch of JSON stuff. 
This can be done after we discuss the neko-top case file.

### ALL the GPU backends!
This can come at the end...

### Adjoint solver
There is still ALOT to do on the adjoint solver regarding the treatment of BCs 
etc! Don't forget to come back to this!

Also, this point should NEVER be removed from the readme, as there are frequent 
updates in the `fluid_scheme.f90`, `fluid_pnpn.f90` etc, and these changes 
should be reflected in `adjoint_scheme.f90` and `adjoint_pnpn.f90` etc.

## Todo (in general):
- There are a ton more minor comments written throughout the source code here, 
so I'm sure some have been missed in this `readme.md`
- Double check all non-dimensionalization and scalings
- Implement the `adjoint_minimum_dissipation_source_term_t` correctly 
(ie, in weak form)
- Restart correctly (or maybe not! in steady calculations it can often be 
faster if your initial condition is the final condition from the previous 
iteration)
- Nothing is freed correctly
- Particularly regarding steady calculations, there's no reason to recalculate 
the adjoint forcing each time, it will always be the same.
- Obviously have MMA based ok KKT not number of itereations 
- Having our own sampler (the one here is hardcoded, but as we introduce new 
eqns and coeffients it will be different)
- Update the `design_t` to include all the cool features Tim put in for 
initializing a design... instead of this hard coded circle 
- We have RAMP currently implemented, but we should also include SIMP, linear  
- Format everything to `Neko` coding standard
- there are ***MANY*** times that I've used `type` instead of `class`, 
particularly regarding the `topopt_design_t`. I suppose these could also be 
calculated differently if we have different `design_variable_t` coming. 
This will need some untangling.
- The mapping functions, eg `linear_mapping` or `ramp_mapping` all use 
`field_math` but it's inefficient and reads poorly. I (Harry) call dibs on 
writing the backend kernels if we decide to switch, because I'm hoping they're 
easy.
 For now, these get called so infrequently it doesn't matter...
