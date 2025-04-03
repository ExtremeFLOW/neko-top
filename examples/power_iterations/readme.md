# PowerIterations {#power-iterations}

We hav a split field, being driven by a baseflow and described by a field of
pertubations.

My goal ios to implement a simulation component which can compute the spectral
values of the flow. This way we should be able to describe the stability of the
system.
The spectral values are computed by the power iterations method.

## Power Iterations

The power iterations method is a method for computing the largest eigenvalue of
a matrix. The method is based on the fact that the largest eigenvalue of a
matrix is the limit of the ratio of the norm of the matrix to the norm of a
vector iteratively multiplied by the matrix.

For a given time step `t` we compute \f$\lambda_t\f$ based on the perturbation
field `u`:

\f[
    \lambda_{t} = \frac{ \sum_{GLL} u_{t} \cdot \bar{u}_{t-1} }
        {\sum_{GLL} \bar{u}_{t-1} \cdot \bar{u}_{t-1}} \\
    \bar{u}_{t} = \frac{ u_{t} }{ ||u_{t}||_2 }\\
    ||u||_2 = \sqrt{ \frac{c}{V} \sum_{GLL} m (u_{t} \cdot u_{t}) }
\f]

Where:
- \f$u_{t}\f$ is the perturbation field at time `t`.
- \f$\bar{u}_{t}\f$ is the normalized perturbation field at time `t`.
- \f$m\f$ is the mass matrix component for the gll points.
- \f$c\f$ is a constant which in Nek5000 is set to 0.5.
- \f$V\f$ is the volume of the domain.



## Implementation

This is already implemented for Nek5000 and we are going to use that as the
baseline for our work here.

- https://github.com/KTH-Nek5000/KTH_Framework
  - Examples/ext_cyl_ARN
  - Examples/ext_cyl_PWI
  - Toolbox/tools/tstpr/pwit/pwit.f
  - Toolbox/tools/tstpr/tstpr.f