# Thermal dissipation from a rod {#easy-E}

This showcases a very crude example of handling heat transfer from a rod to the
surrounding fluid. The rod have an assigned temperature and the fluid has an
inflow temperature. 

Since Neko only support a single medium for now, the fluid is represented by a
cylindrical region where the temperature is assigned manually and the flow is 
restricted through a Brinkman term.

The Peclet number is set to 2000 and the Reynolds number is set to 2000.

The rod temperature is set to 3 different values, 1, 10 and 20.

Please run the prepare script which will setup the domain.

```bash
./prepare
```

Then run the simulation with the following command:

```bash
neko easy-E_[N].case
```

In order to avoid numerical issues, we ramp the temperature to the target value
linearly over the first 1% of the total simulation time.

The fluid is outputted 250 times, which support an animation of the flow with 25
frames per second without time dilation.

## Problem setup:

Domain: 3D, 4x1x1 box.
Discretization: 40x10x10 cells.
Time: 10.0 seconds.
Rod temperature: 1, 10 and 20.
Rod location: 1.0, 0.5, 0.0.
Rod geometry:
- Constructed from cylinder with a sphere at the internal end.
- Cylinder: 0.1 radius, 0.5 length, placed vertically.
- Sphere: 0.1 radius.

```
----------------------------------------
->                                    ->
->                                    ->
->        O                           ->
->       | |                          ->
->       | |                          ->
----------------------------------------
```

