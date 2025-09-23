# Static mixer design.

Please note this example depends on a (as of 20th of February 2024) unreleased
feature of Neko. Please checkout the `feature/immersed_boundary` pull request
branch in the Neko library to run this example.


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
