# A Simple 1D Beam Optimization Test Case for MMA

This example demonstrates a simple test case for the Method of Moving Asymptotes (MMA), applied to a **1D cantilever beam** discretized into elements. The design variable is the **beam height** `h` at each element, and the objective is to **maximize stiffness (minimize tip deflection)** while satisfying **stress constraints**.

The goal is to test the MMA implementation on a small, analytically tractable structural optimization problem, without requiring a full 3D FEM solver.

---

## Problem Description

In the `driver.f90` file, the beam is divided into `n` elements, each of length $L_e$. The **design variable** is the height $h(k)$ of element $k$.

> Note: Each element has **one design variable**, corresponding to its height.

### Beam Properties

- Width $b$ (constant)
- Young’s modulus $E$ (constant)
- Element height $h(k)$ (design variable)
- Moment of inertia:  
  $$
  I(k) = \frac{b \cdot h(k)^3}{12}
  $$
- Maximum allowable stress: $\sigma_\text{max}(k)$
- **Bounds** on height: $h_{min} \leq h(k) \leq h_{max}$

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





## Notes

- This example is **1D** and uses **analytical formulas** for both the objective and constraints.
- It is intended as a **minimal, testable example** for verifying your MMA implementation performance.
- All variables and sensitivities are scaled to work with the MMA optimizer directly.
