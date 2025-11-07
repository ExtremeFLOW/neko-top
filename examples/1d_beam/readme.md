# A Simple 1D Beam Optimization Test Case for MMA {#beam}

This example demonstrates a simple test case for the Method of Moving Asymptotes (MMA) as implemented in **NEKO-TOP**.
The problem is a **1D cantilever beam**, discretized into multiple elements, where the design variable is the **beam height** `h` at each element. The objective is to **maximize stiffness (minimize tip deflection)** while **minimizing beam weight** and enforcing **stress constraints**.
The number of design variables (1D beam elements) matches the number of Gauss–Lobatto–Legendre (GLL) points in the 3D spectral element mesh provided in the `box.nmsh` file from the examples repository. This creates a one-to-one mapping between the GLL points in the spectral element mesh and the 1D beam elements.
This setup allows to use the **Neko’s mesh partitioning** for parallel execution: when running on multiple MPI nodes, Neko distributes the GLL points across nodes, which directly translates to distributing the 1D beam elements in the parallel system.


## Problem Description

In the `driver.f90` file, the beam is divided into `n` elements, each of length $L_e$. The **design variable** is the height $h(k)$ of element $k$.

> Note: Each element has **one design variable**, corresponding to its height.

### Beam Properties

- **Structure**: Cantilever beam of length `L = 2.0 m` and width `b = 0.02 m`.
- **Material**: Density `ρ = 7800 kg/m³`, Young’s modulus `E = 210 GPa`.
- **Loading**: Tip load `P = 1000 N`.
- **Design Variables**: Element-wise beam height `h` ∈ [`h_min = 0.005 m`, `h_max = 0.05 m`]. Moment of inertia (the beam is assumed to be rectangular with depth b and hight h):  

  \f[
  I(k) = \frac{b \cdot h(k)^3}{12}
  \f]

- **Objectives**:
  1. Minimize **tip deflection** (normalized by maximum allowed `u_tip_max = 0.25 m`).
  2. Minimize **beam weight**.
- **Constraints**:
  - Element-wise stress must satisfy `σ ≤ σ_max = 250 MPa`.


### Design Variable Mapping

The optimization variable \f$ x(k) \f$ is **normalized** between 0 and 1. The **physical element height** \f$ h(k) \f$ is obtained by linear projection:

\f[
h(k) = h_{\min} + x(k) \cdot (h_{\max} - h_{\min})
\f]
---

### Tip Deflection

For a cantilever beam under a **point load** \f$ P \f$ at the free end, the **tip deflection** \f$ u_\text{tip} \f$ is computed analytically as:

\f[
u_\text{tip} = \sum_{k=1}^{n} \frac{\Delta_k}{E \cdot I(k)}
\f]

with

\f[
\Delta_k = \frac{(L_\text{total} - x_{k-1})^3 - (L_\text{total} - x_k)^3}{3}
\f]

and

\f[
x_k = k \cdot L_e
\f]

where \f$ x_k \f$ is the end coordinate of element `k`.

The **sensitivity** of the tip deflection with respect to the element height \f$ h(k) \f$ is computed analytically as:

\f[
\frac{\partial J}{\partial h(k)} = - \frac{\Delta_k \cdot 3 b h(k)^2}{12 \cdot E \cdot I(k)^2} = - \frac{\Delta_k \cdot 3}{E \cdot b \cdot h(k)^4}
\f]

In terms of the normalized variable \f$ x(k) \f$:

\f[
\frac{\partial J}{\partial x(k)} = \frac{\partial J}{\partial h(k)} \cdot (h_{\max} - h_{\min})
\f]

---

## Stress Constraints

The **stress constraint** ensures that the maximum bending stress in each element does not exceed the allowable limit:

\f[
\sigma_e = \frac{M_e \cdot c_e}{I_e} \leq \sigma_\text{max}(k)
\f]

where:

- \f$ M_e = P \cdot (L_\text{total} - x_e) \f$ is the **bending moment** at element `e`
- \f$ c_e = \frac{h(k)}{2} \f$ is the distance from the **neutral axis**
- \f$ I_e = \frac{b \cdot h(k)^3}{12} \f$ is the **moment of inertia**

The **constraint function** used in MMA is:

\f[
g_\sigma(k) = \frac{\sigma_e}{\sigma_\text{max}(k)} - 1 \leq 0
\f]

The **sensitivity** of the stress constraint with respect to the element height \f$ h(k) \f$ is:

\f[
\frac{\partial g_\sigma}{\partial h(k)} =
\frac{\partial}{\partial h(k)} \left( \frac{M_e \cdot c_e}{I_e \cdot \sigma_\text{max}(k)} - 1 \right)
= \frac{M_e}{\sigma_\text{max}(k)} \left( \frac{1}{2 I_e} - \frac{c_e \cdot 3 b h(k)^2 / 12}{I_e^2} \right)
\f]

In terms of the normalized variable \f$ x(k) \f$:

\f[
\frac{\partial g_\sigma}{\partial x(k)} = \frac{\partial g_\sigma}{\partial h(k)} \cdot (h_{\max} - h_{\min})
\f]

## Problem Formulation

**Minimize:**

\f[
f(\mathbf{x}) = w_1 \cdot \frac{u_{\text{tip}}(\mathbf{x})}{u_{\text{tip}}^{\max}} + w_2 \cdot \frac{m(\mathbf{x})}{m_{\max}}
\f]

**Subject to:**

\f[
g_j(\mathbf{x}) = \frac{\sigma_j(\mathbf{x})}{\sigma_{\max}} - 1 \le 0, \quad j = 1, \dots, N_\text{constraints}
\f]

\f[
0 \le x_i \le 1, \quad i = 1, \dots, N_\text{elements}
\f]

**Where:**


- \f$ \mathbf{x} \f$ are the normalized design variables.
- \f$ u_{\text{tip}} \f$ is the tip deflection.
- \f$ m \f$ is the total beam mass.
- \f$ \sigma_j \f$ is the bending stress at the `j`-th constraint location.
- \f$ N_\text{constraints} \f$ can be set in `driver.f90` using `num_constraints = 10` as default.
- The weights \f$ w_1, w_2 \f$ for the objective function can be set in `driver.f90` using
`call deflection%init_from_components("tip_deflection", des, 1.0_rp)` and 
`call beamweight%init_from_components("beamweight", des, 1.0_rp)`
with `1.0_rp` as the default weight for both objectives.

## Notes

- This example is **1D** and uses **analytical formulas** for both the objective and constraints.
- It is intended as a **minimal, testable example** for verifying your MMA implementation performance.
- Beam parameters can be updated at the beginning of `example_problem.f90`:
```fortran
  ! ========================================================================== !
  ! Global beam parameters
  real(rp), public :: L_total   = 2.0_rp      ! total length (m)
  real(rp), public :: b         = 0.02_rp     ! width (m)
  real(rp), public :: rho       = 7800.0_rp   ! density (kg/m3)
  real(rp), public :: E         = 210.0e9_rp  ! Young's modulus (Pa)
  real(rp), public :: P         = 1000.0_rp   ! tip load (N)
  real(rp), public :: h_min     = 0.005_rp    ! min height (m)
  real(rp), public :: h_max     = 0.05_rp     ! max height (m)
  real(rp), public :: u_tip_max = 0.25_rp     ! max tip deflection (m)
  ! ========================================================================== !
```
- Number of constraints and their distribution can be altered in `driver.f90` file:
```fortran
  ! Set up distributed stress constraints
  num_constraints = 10
  num_constraint_partitions=10
```
where the beam is divided into `num_constraint_partitions` sections, and then the `num_constraints` constraints are distributed evenly across the first elements of each section.
