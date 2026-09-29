# Objectives and constraints {#objectives_and_constraints}

\tableofcontents

`Neko-top` allows the user to solve constrained optimization problems,
often times but not limited to, topology optimization problems involving fluid 
mechanics.

Objectives and constraints enter the `.case` file in `neko-top` as lists

```
{
    "version": 1.0
    "case": {},
    "optimization": {
        "objectives": [],
        "constraints": []
    }
}
```

In `neko-top` multiple objectives can be prescribed in a list, resulting in
a multi-objective optimization problem which is handled by a weighted sum of all
prescribed objectives
\f[
\mathcal{F} = \sum_i w_i \mathcal{F}_i,
\f]
where \f$\mathcal{F}i\f$ is an objective value \f$w_i\f$ is a prescribed weight.

Multiple constraints on the other hand are also entered in a list but are
handled through the MMA functionality discussed in 
[The Method of Moving Asymptotes (MMA)](@ref MMA).

The following objectives

1. [Viscous dissipation](@ref objective_dissipation)
2. [Brinkman dissipation](@ref objective_velocity_penalty)
3. [Scalar mixing](@ref objective_scalar_mixing)

and constraints

1. [Volume constraint](@ref constraint_volume)

have currently been implemented in `neko-top`·

## Objectives {#objectives}

### Time integration for unsteady objectives {#objectives_time_integration}
For unsteady simulations, the instantaneous objective is averaged over time,
\f[
\mathcal{F} = \frac{1}{|W|} \int_W f(t) \, dt,
\f]
where \f$W\f$ is the objective's own time window and \f$|W|\f$ its length.
Averaging rather than integrating keeps the objective on the same scale as the
instantaneous quantity it is built from.

The window can be restricted with the optional input parameters `start_time`
and `end_time`. Their defaults are `0.0` and `+\infty`, so the full simulated
horizon is used unless a smaller window is prescribed.

The objective is sampled once per completed timestep, and \f$|W|\f$ is the
total length of time actually sampled. A window is therefore truncated to the
part of it that was simulated, and the reported number is the mean over what
was simulated. In particular the average is normalised by the window, not by
the length of the run, so **a windowed objective does not change when the run
is made longer**: an objective windowed to \f$[2.5, 6.0]\f$ reports the same
value whether the simulation stops at \f$t = 6\f$ or continues to \f$t = 20\f$.

A window that never overlaps the simulated interval accumulates nothing. The
objective then reports zero and a warning is written to the log.

\warning Setting `end_time` to exactly the simulation's own `end_time` is not
the same as leaving it at its `+\infty` default. The adjoint forcing terms that
objectives register are gated by `source_term_t`'s time window, which compares
against the simulation time without a tolerance, so the final step -- where the
accumulated time lands a few ULP above the prescribed `end_time` -- silently
loses its forcing and the sensitivity is wrong. Leave `end_time` unset unless a
genuinely shorter window is wanted.

The window only applies to unsteady problems. A steady problem evaluates its
objectives once, on the converged field, so `start_time` and `end_time` have no
effect there.

The two paths agree wherever they are asking the same question. Once a problem
has reached a steady state, an unsteady objective averaged over a window inside
the converged part of the run reports what the steady path reports for the same
problem: the field no longer changes, so its average over the window is its
converged value. `tests/unit/objectives/steady_unsteady_converged` is exactly
this check, and it holds to the iterative solvers' own round-off. Note that
`steady_simcomp` freezes the fluid on convergence but not the scalar, so set
`scalar_coupled` when a scalar objective is involved, or the objective may be
evaluated on a scalar that has not settled.

They also agree trivially at the end of any run, converged or not: an unsteady
objective windowed to the final timestep alone reports exactly what the steady
path reports, because a one-sample average is the sample. Selecting that single
step is easier than it looks. The time loop stops at the *first* step reaching
the simulation's `end_time`, so unless `end_time` falls exactly on a step, the
run overshoots it by up to one `dt` -- with `end_time` \f$0.0475\f$ and `dt`
\f$0.005\f$ it takes ten steps and finishes at \f$t = 0.05\f$. A window whose
`start_time` is the simulation's own `end_time` therefore captures precisely
that last step and nothing before it. Keeping `end_time` off a step boundary
also removes any dependence on which side of the comparison round-off in the
accumulated time falls, which otherwise decides whether the run takes ten steps
or eleven.

For example,

```json
{
    "type": "scalar_mixing",
    "weight": 1.0,
    "start_time": 2.5,
    "end_time": 6.0
}
```

accumulates the scalar-mixing objective only over the interval
\f$2.5 \le t \le 6.0\f$, and reports its mean over that interval.

For the underlying adjoint formulation and how objective functions generate
adjoint forcing terms, please refer to
[Adjoint sensitivity analysis](@ref adjoint).

### Viscous dissipation {#objective_dissipation}

This objective is used to either minimize or maximize the viscous dissipation.
It takes the form
\f[
\mathcal{F} = \frac{1}{|\Omega_\text{obj}|}\int_{\Omega_\text{obj}}
\frac{\mu}{2} |\nabla \mathbf{u}|^2 d\Omega,
\f]
where \f$\mathbf{u}\f$ is the fluid velocity, \f$\Omega_\text{obj}\f$ is the
objective domain, \f$\mu\f$ is the dynamic viscosity, and
\f$|\nabla \mathbf{u}|^2\f$ denotes the Frobenius norm of the velocity
gradient. Here \f$|\Omega_\text{obj}|\f$ denotes the volume of the objective
domain.

The objective can be selected by prescribing `"type": "viscous_dissipation"`
and has the following input parameters:


| Name | Description  | Admissible values | Default value |
|------|--------------|-------------------|---------------|
| `weight`| The weight used in the objective. | Real | `1.0` |
| `mask_name` | The name of the `point_zone` indicating \f$\Omega_\text{obj}\f$. | String | `""`|
| `name`| The name that will appear in `objective_data.csv` | String | `Dissipation`|
| `start_time` | Start of the active time window for unsteady accumulation. | Real | `0.0` |
| `end_time` | End of the active time window for unsteady accumulation. | Real | `+\infty` |

### Brinkman dissipation {#objective_velocity_penalty}
In the works of [A. Gersborg-Hansen et al. (2005)](https://link.springer.com/article/10.1007/s00158-004-0508-7)
an objective function of the form
\f[
\mathcal{F} = \frac{1}{|\Omega_\text{obj}|}\int \frac{1}{2} 
\left[
\underset{I}{\nabla \mathbf{u} \cdot \left(\nabla \mathbf{u} + (\nabla \mathbf{u})^T \right)}
+ \underset{II}{\chi \mathbf{u} \cdot \mathbf{u} }
\right] d\Omega,
\f]

where \f$\mathbf{u}\f$ is the fluid velocity and \f$\chi\f$ the Brinkman amplitude
was used, claiming:
> Term I is half of the part of the dissipation function that is
> associated with the in-plane components of the stretching
> tensor (cf. Currie (2003)), while term II is half the part associated 
> with the out-of-plane components. The latter part
> arises from the parabolic velocity profile in the lubrication
> theory (2). From an optimization perspective, term II links
> the cost function directly to the two-dimensional velocity
> field.

[Later works](https://doi.org/10.1016/j.compfluid.2022.105387) have argued that
this second term can also be used to penalize intermediate values of the
design indicator and promote binary designs. Hence, this second term is considered
a "velocity penalty" in `neko-top` and takes the form
\f[
\mathcal{F} = \frac{1}{|\Omega_\text{obj}|}\int_{\Omega_\text{obj}} 
\frac{1}{2} \chi \mathbf{u}^2 d\Omega.
\f]

The objective can be selected by prescribing `"type": "brinkman_dissipation"`
and has the following input parameters:

\note the naming convention of `"brinkman_dissipation"` comes from the original claim
based on lubrication theory written by Gersborg-Hansen et al.


| Name | Description  | Admissible values | Default value |
|------|--------------|-------------------|---------------|
| `weight`| The weight used in the objective. | Real | `1.0` |
| `mask_name` | The name of the `point_zone` indicating \f$\Omega_\text{obj}\f$. | String | `""`|
| `name`| The name that will appear in `objective_data.csv` | String | `Out of plane stresses`|
| `dealias_forcing`| If dealiasing should be applied to adjoint forcing term | logical | `.true.`|
| `dealias_sensitivity`| If dealiasing should be applied to sensitivity contribution | logical | `.true.`|
| `start_time` | Start of the active time window for unsteady accumulation. | Real | `0.0` |
| `end_time` | End of the active time window for unsteady accumulation. | Real | `+\infty` |

### Scalar mixing {#objective_scalar_mixing}

This objective is used to either minimize or maximize the mixing of a passive
scalar.
It takes the form

\f[
\mathcal{F} = \frac{1}{|\Omega_\text{obj}|}\int_{\Omega_\text{obj}} 
\frac{1}{2} (\phi - \phi_\text{ref})^2 d\Omega,
\f]

where \f$\phi\f$ is the scalar field, \f$\Omega_\text{obj}\f$ is the
objective domain, \f$|\Omega_\text{obj}|\f$
denotes the volume of the objective domain and \f$ \phi_\text{ref}\f$ is
a target concentration.

The objective can be selected by prescribing `"type": "scalar_mixing"` 
and has the following input parameters:


| Name | Description  | Admissible values | Default value |
|------|--------------|-------------------|---------------|
| `weight`| The weight used in the objective. | Real | `1.0` |
| `mask_name` | The name of the `point_zone` indicating \f$\Omega_\text{obj}\f$. | String | `""`|
| `target_concentration` | \f$\phi_\text{ref}\f$ in the above equation. | Real | `0.5`|
| `name`| The name that will appear in `objective_data.csv` | String | `Scalar Mixing`|
| `start_time` | Start of the active time window for unsteady accumulation. | Real | `0.0` |
| `end_time` | End of the active time window for unsteady accumulation. | Real | `+\infty` |

## Constraints {#constraints}

### Volume constraint {#constraint_volume}

This constraint is used to constrain the volume of the design in the domain.
It takes the form

\f[
\mathcal{C} = \frac{1}{|\Omega_\text{opt}|}\int_{\Omega_\text{opt}} \rho d\Omega,
\f]

where \f$\rho\f$ is the material indicator, \f$\Omega_\text{opt}\f$ is the
optimization domain and \f$|\Omega_\text{opt}|\f$
denotes the volume of the optimization domain.

The constraint can be used to enforce either a minimum or maximum volume, i.e.
\f$ \mathcal{C} > \mathcal{C}_\text{min} \f$ or \f$ \mathcal{C} < \mathcal{C}_\text{max} \f$.

\note Currently the volume constraint can only be applied to the unfiltered
material indicator function, but in the future we aim to allow it to be prescribed
to intermediate stages of the mapping cascade.

The constraint can be selected by prescribing `"type": "volume"` 
and has the following input parameters:


| Name | Description  | Admissible values | Default value |
|------|--------------|-------------------|---------------|
| `limit` | \f$ \mathcal{C}_\text{min} \f$ or \f$ \mathcal{C}_\text{max} \f$  in the above equation | Real | - |
| `is_max` | Indicate whether a minimum or maximum volume constraint should be applied. | `.true.` or `.false.` | `.false.` |
| `mask_name` | The name of the `point_zone` indicating \f$\Omega_\text{obj}\f$. | String | `""`|
| `name`| The name that will appear in `objective_data.csv` | String | `Volume constraint`|
| `mapping`| A potential to, for instance, compute the volume based on a filtered design. For more information please refer to @ref mixer_mapping| Json | `""`|
