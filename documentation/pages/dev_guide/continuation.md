# Continuation {#continuation}

This document describes the continuation framework introduced in Neko-TOP to
support gradual variation of optimization parameters during the optimization
loop.

The implementation is based on the module `continuation_scheduler`.


## 1. Overview

The continuation framework allows selected parameters to change over
optimization iterations according to a predefined schedule.

Instead of keeping parameters fixed, users can define:

- A sequence of values for a parameter
- How many iterations each value should persist

This adds the flexibility to guide the optimizer using different values for
various parameters through the optimization loop.


## 2. Core Concept

Each continuation parameter is defined as:

- A **name**
- A **target pointer** (link to actual runtime variable)
- A **list of values**
- A number of iterations per value

At runtime at the end of each optimization iteration:

- The scheduler updates all registered parameters
- Each parameter selects a value based on the current iteration and assigned
  list of values


## 3. Module: `continuation_scheduler`

### 3.1 Main Types

#### continuation_parameter_t

Represents a single continuation-controlled variable.

```fortran
type continuation_parameter_t
   character(len=:), allocatable :: name
   real(rp), pointer :: target => null()
   real(rp), allocatable :: values(:)
   integer :: iterations_per_value = 0
end type
```


### 3.2 Continuation Scheduler
This is a `continuation_scheduler_t`, which allows registration of multiple
continuation parameters and provides an update functionality that updates all
registered parameters.
```fortran
type continuation_scheduler_t
     type(continuation_parameter_t), allocatable :: params(:)
     integer :: default_iterations = 1
   contains
     procedure :: init => init_scheduler
     procedure :: free => free_scheduler
     procedure :: json_get_or_register
     procedure :: register_parameter
     procedure :: update
     procedure :: get_param_name
     procedure :: get_n_params
end type
```

### 3.3 Global Instance
A global instance of the continuation scheduler is provided as
`nekotop_continuation`. It is initialized in the driver at the beginning of the
simulation. Any other module that requires a continuation strategy for one or
more of its parameters can register them through this global object. Note that
`nekotop_continuation%update(iter)` is called at each iteration in the
optimization loop.

```fortran
type(continuation_scheduler_t), public, target :: nekotop_continuation
```

## 4. Initialization

```fortran
call nekotop_continuation%init(json)
```
Reads `optimization.solver.continuation_iterations` in the case file as the
default frequency of updating all parameters registered in the scheduler. Note
that individual frequency for can be registered in the json file as
`{parameters_name}_iterations` in the case file.


## 5. Registering Parameters
```fortran
call nekotop_continuation%register_parameter(name, target, values, iterations)
```
It is possible to use `json_get_or_register` function in the scheduler to read
the parameter form the json file. Note that if the parameter in the json is an
array, then it will automatically register in the scheduler as well. If a single
scalar is given in the json, then this function will work as a single call of
`json_get_or_default`.
```fortran
call nekotop_continuation%json_get_or_register(json, name, target, default_value)
```

## 6. Update Mechanism
```fortran
call nekotop_continuation%update(iter)
```
`idx = (iter - 1) / iterations_per_value + 1` `target = values(idx)`


## 7. JSON Integration
The following utility function helps to read a parameter with continuation
strategy form the case file:
```fortran
call json_get_with_continuation(json, name, values, default_value, iter)
```
If the parameter is a single scalar, this will be a constant parameter in the
optimization loop and will not be registered to `nekotop_continuation`. If it is
an array, then it will be registered with these values with
`optimization.solver.continuation_iterations` as the default value for its
frequency, unless the individual frequency is defined for that particular
parameter as `{parameters_name}_iterations`.


## 8. Integration Points

### Optimizer loop
```fortran
call nekotop_continuation%update(current_iteration)
```
### Logging
Continuation parameters are appended to log headers and values in the
`optimizer`.


## 9. Free()
The scheduler should be freed at the end of the simulation in the driver.
```fortran
call nekotop_continuation%free()
```

## 10. Example
```json
"beta": [1.0, 5.0, 10.0],
"beta_iterations": 50
```

## 11. Summary

Continuation enables gradual parameter variation during optimization to improve
robustness and convergence.
