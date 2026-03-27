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
