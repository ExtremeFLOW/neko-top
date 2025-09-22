# A Simple 1D Beam Optimization Test Case for MMA

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
 $I(k) = \frac{b \cdot h(k)^3}{12}$
- **Objectives**:
  1. Minimize **tip deflection** (normalized by maximum allowed `u_tip_max = 0.25 m`).
  2. Minimize **beam weight**.
- **Constraints**:
  - Element-wise stress must satisfy `σ ≤ σ_max = 250 MPa`.
  - Stress constraints are distributed across MPI processes.


### Design Variable Mapping

The optimization variable $x(k)$ is **normalized** between 0 and 1. The **physical element height** $h(k)$ is obtained by linear projection:

$$
h(k) = h_{min} + x(k) \cdot (h_{max} - h_{min})
$$

---

### Tip Deflection

For a cantilever beam under a **point load** $P$ at the free end, the **tip deflection** $u_\text{tip}$ is computed analytically as:

$$
u_\text{tip} = \sum_{k=1}^{n} \frac{\Delta_k}{E \cdot I(k)}
$$

with

$$
\Delta_k = \frac{(L_\text{total} - x_{k-1})^3 - (L_\text{total} - x_k)^3}{3}
$$

and

$$
x_k = k \cdot L_e
$$

where $x_k$ is the end coordinate of element $k$.

The **sensitivity** of the tip deflection with respect to the element height $h(k)$ is computed analytically as:

$$
\frac{\partial J}{\partial h(k)} = - \frac{\Delta_k \cdot 3 b h(k)^2}{12 \cdot E \cdot I(k)^2} = - \frac{\Delta_k \cdot 3}{E \cdot b \cdot h(k)^4}
$$

In terms of the normalized variable $x(k)$:

$$
\frac{\partial J}{\partial x(k)} = \frac{\partial J}{\partial h(k)} \cdot (h_{max} - h_{min})
$$

---

## Stress Constraints

The **stress constraint** ensures that the maximum bending stress in each element does not exceed the allowable limit:

$$
\sigma_e = \frac{M_e \cdot c_e}{I_e} \leq \sigma_\text{max}(k)
$$

where:

- $M_e = P \cdot (L_\text{total} - x_e)$ is the **bending moment** at element $e$
- $c_e = \frac{h(k)}{2}$ is the distance from the **neutral axis**
- $I_e = \frac{b \cdot h(k)^3}{12}$ is the **moment of inertia**

The **constraint function** used in MMA is:

$$
g_\sigma(k) = \frac{\sigma_e}{\sigma_\text{max}(k)} - 1 \leq 0
$$

The **sensitivity** of the stress constraint with respect to the element height $h(k)$ is:

$$
\frac{\partial g_\sigma}{\partial h(k)} = \frac{\partial}{\partial h(k)} \left( \frac{M_e \cdot c_e}{I_e \cdot \sigma_\text{max}(k)} - 1 \right)
= \frac{M_e}{\sigma_\text{max}(k)} \left( \frac{1}{2 I_e} - \frac{c_e \cdot 3 b h(k)^2 / 12}{I_e^2} \right)
$$

In terms of the normalized variable $x(k)$:

$$
\frac{\partial g_\sigma}{\partial x(k)} = \frac{\partial g_\sigma}{\partial h(k)} \cdot (h_{max} - h_{min})
$$


## Problem Formulation

**Minimize:**

$$
f(\mathbf{x}) = w_1 \cdot \frac{u_{tip}(\mathbf{x})}{u_{tip}^{max}} + w_2 \cdot \frac{m(\mathbf{x})}{m_{max}}
$$

**Subject to:**

$$
g_j(\mathbf{x}) = \frac{\sigma_j(\mathbf{x})}{\sigma_{max}} - 1 \le 0, \quad j = 1, \dots, N_\text{ constraints}
$$

$$
0 \le x_i \le 1, \quad i = 1, \dots, N_\text{ elements}
$$

**Where:**

- $\mathbf{x}$ are the normalized design variables.
- $u_{tip}$ is the tip deflection.
- $m$ is the total beam mass.
- $\sigma_j$ is the bending stress at the $j$-th constraint location.
- $N_\text{ constraints} can be set in the driver.f90 in lines:
```
  ! Set up distributed stress constraints
  num_constraints = 10
```
- the weights $$w_1, w_2$$ for the objective function can be set in the driver.f90 in lines:
```
  call deflection%init_from_components("tip_deflection", des, 1.0_rp)
  call beamweight%init_from_components("beamweight", des, 1.0_rp)
```
where 1.0_rp is the weight for each objective.


## Notes

- This example is **1D** and uses **analytical formulas** for both the objective and constraints.
- It is intended as a **minimal, testable example** for verifying your MMA implementation performance.
