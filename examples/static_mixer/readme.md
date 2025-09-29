# Static mixer example used for benchmarking {#static-mixers}

This example contain the static mixer examples of our code. So far only the
design by Casper Schousboe Andreasen is here for reference.

## PetSC Static mixer design

Details:

- Domain: Duct of dimension 4 x 1 x 1
- BCs:
    - Inlet (x=0): Parabolic velocity profile, u_max = 1, v = w = 0
      (Currently constant)
    - Two inlets for temperature:
        - Bottom half (z < 0.5): Dirichlet, T = 0
        - Top half (z >= 0.5): Dirichlet, T = 1
    - Walls (y=0, y=1, z=0, z=1): No-slip, Neumann for thermal.
    - Outlet (x=4): Pressure outflow.
    - Initial condition: Fluid at constant velocity u = 1, v = w = 0
- Re = 1000
- Brinkman IBM for solid obstacles.
    - Max restriction = 100
    - Min restriction = 0
