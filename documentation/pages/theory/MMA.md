# The Method of Moving Asymptotes (MMA) {#MMA}

\tableofcontents

The Method of Moving Asymptotes (MMA) is the optimization algorithm currently
used in `neko-top`. In a case file, MMA is configured in `optimization.solver`,
with MMA-specific settings under `optimization.solver.mma`.

```json
{
  "optimization": {
    "solver": {
      "type": "mma",
      "mma": {
        // MMA options go here
      }
    }
  }
}
```

## MMA Options overview

| Json entry   | Description                                            | Admissible Values     | Default value                                      |
| ------------ | ------------------------------------------------------ | --------------------- | -------------------------------------------------- |
| `max_iter`   | Maximum number of iterations in the MMA subsolver.     | Integer, `> 0`        | `100`                                              |
| `epsimin`    | Minimum barrier parameter used by MMA.                 | Real, `> 0`           | `1.0e-9 * sqrt(m + n_global)`                      |
| `asyinit`    | Initial asymptote distance factor.                     | Real, `> 0`           | `0.2`                                              |
| `asyincr`    | Asymptote increase factor.                             | Real, `> 0`           | `1.05`                                             |
| `asydecr`    | Asymptote decrease factor.                             | Real, `> 0`           | `0.65`                                             |
| `backend`    | MMA execution backend.                                 | `"cpu"` or `"device"` | `"cpu"` (or `"device"` when `NEKO_BCKND_DEVICE=1`) |
| `subsolver`  | Interior-point variant used in MMA subproblems.        | `"dip"` or `"dpip"`   | `"dip"`                                            |
| `xmin`       | Global lower bound applied to all design variables.    | Real                  | `0.0`                                              |
| `xmax`       | Global upper bound applied to all design variables.    | Real                  | `1.0`                                              |
| `a0`         | MMA scalar coefficient \f$a_0\f$.                      | Real                  | `1.0`                                              |
| `a`          | MMA vector coefficient \f$a_i\f$ (applied uniformly).  | Real                  | `0.0`                                              |
| `c`          | MMA vector coefficient \f$c_i\f$ (applied uniformly).  | Real                  | `100.0`                                            |
| `d`          | MMA vector coefficient \f$d_i\f$ (applied uniformly).  | Real                  | `0.0`                                              |
| `scale`      | Constraint scaling target value.                       | Real, `> 0`           | `10.0`                                             |
| `auto_scale` | Enable adaptive scaling of constraints each iteration. | `.true.` or `.false.` | `.false.`                                          |

### MMA subproblem and asymptote options
- `max_iter` and `epsimin` define limits for the MMA subsolver iterations and
  barrier parameter.
- `asyinit`, `asyincr`, and `asydecr` control moving-asymptote updates.
- `subsolver` selects the interior-point variant used inside MMA (`dip`/`dpip`).

### Design bounds and MMA coefficients
- `xmin` and `xmax` set global lower/upper bounds for all design variables.
- `a0`, `a`, `c`, and `d` are coefficients in the standard MMA subproblem
  formulation.

### Backend and scaling options
- `backend` selects whether MMA update/KKT operations run on CPU or device
  routines.
- `scale` sets the constraint scaling magnitude.
- `auto_scale` enables adaptive scaling from the current constraint value at
  each iteration.
