# Design Types {#design_types}

\tableofcontents

This page documents the JSON options for `optimization.design`.

## Overview {#design_types-overview}

The design section always lives under:

```json
"optimization": {
  "design": {
    "type": "..."
  }
}
```

Currently supported design types are:

- `"brinkman"`
- `"simple"`

## Brinkman Design {#design_types-brinkman}

`"type": "brinkman"` is the main design used with the mapping cascade and
Brinkman source-term coupling.

### JSON keys {#design_types-brinkman-keys}

| Key | Type | Required | Default | Notes |
|---|---|---|---|---|
| `type` | string | Yes | - | Must be `"brinkman"` |
| `name` | string | No | `"Brinkman Design"` | Design name |
| `domain.type` | string | No | `"full"` | `"full"` or `"point_zone"` |
| `domain.zone_name` | string | Cond. | - | Required when `domain.type = "point_zone"` |
| `dealias` | bool | No | `true` | Controls dealiasing for Brinkman source terms |
| `verbose_design` | bool | No | `false` | Write all forward mapping stages |
| `verbose_sensitivity` | bool | No | `false` | Write all backward sensitivity stages |
| `output_precision` | string | No | `"sp"` | `"sp"` or `"dp"` for design/sensitivity fld output |
| `mapping` | array | No | omitted | Mapping cascade definition, see [Mapping cascade](@ref mapping_cascade) |
| `initial_distribution` | object | No | omitted | Initial design field setup (details below) |

### `initial_distribution` options {#design_types-brinkman-initial}

If omitted, the design is initialized to zero.

#### Uniform

```json
"initial_distribution": {
  "type": "uniform",
  "value": 0.0
}
```

| Key | Type | Required | Default |
|---|---|---|---|
| `type` | string | Yes | - |
| `value` | real | Yes | - |

#### Point zone

```json
"initial_distribution": {
  "type": "point_zone",
  "base_value": 0.0,
  "zone_name": "my_zone",
  "zone_value": 1.0
}
```

| Key | Type | Required | Default |
|---|---|---|---|
| `type` | string | Yes | - |
| `base_value` | real | Yes | - |
| `zone_name` | string | Yes | - |
| `zone_value` | real | Yes | - |

#### Field file

```json
"initial_distribution": {
  "type": "field",
  "file_name": "design0.f00001",
  "interpolate": false,
  "tolerance": 1.0e-6,
  "mesh_file_name": "none"
}
```

| Key | Type | Required | Default | Notes |
|---|---|---|---|---|
| `type` | string | Yes | - | Must be `"field"` |
| `file_name` | string | Yes | - | Field sample file (e.g. `design0.f00001`) |
| `interpolate` | bool | No | `false` | Enable interpolation from another mesh |
| `tolerance` | real | No | `1.0e-6` | Interpolation tolerance |
| `mesh_file_name` | string | No | `"none"` | Optional mesh/source field file for interpolation |

## Simple Design {#design_types-simple}

`"type": "simple"` is a lightweight design representation.

\note `simple` design currently supports only single-rank execution.
\warning In the current implementation, the `box` sizing keys are read with the
fully-qualified paths `optimization.design.domain.nx/ny/nz/limits`.

### JSON keys {#design_types-simple-keys}

| Key | Type | Required | Default | Notes |
|---|---|---|---|---|
| `type` | string | Yes | - | Must be `"simple"` |
| `name` | string | No | `"Simple Design"` | Design name |
| `domain.type` | string | Yes | - | Currently supports `"box"` |
| `domain.nx` | integer | Yes | - | Number of cells in x-direction |
| `domain.ny` | integer | Yes | - | Number of cells in y-direction |
| `domain.nz` | integer | Yes | - | Number of cells in z-direction |
| `domain.limits` | real[6] | Yes | - | `[xmin, xmax, ymin, ymax, zmin, zmax]` |

Example:

```json
"design": {
  "type": "simple",
  "name": "Simple Design",
  "domain": {
    "type": "box",
    "nx": 16,
    "ny": 16,
    "nz": 16,
    "limits": [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
  }
}
```
