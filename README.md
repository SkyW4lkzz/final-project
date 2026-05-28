# LBM Final Project

This project contains a D2Q9 lattice Boltzmann solver with a BGK collision
operator and optional OpenMP parallelism. It currently supports two cases:
plane Poiseuille flow and flow past a circular cylinder.

- periodic boundary condition in the streamwise direction
- bounce-back wall boundary condition at the top and bottom walls
- Guo forcing term for the pressure-gradient/body-force source
- numerical profile compared with the analytical parabolic velocity profile
- runtime reported in MLUPS for later performance comparison

## Structure

```text
include/lbm.hpp      Shared D2Q9 constants, config, data types, declarations
src/main.cpp         Command-line parsing and case dispatch
src/lbm.cpp          Common lattice utilities, initialization, output helpers
src/poiseuille.cpp   Poiseuille update loop, analytical profile, error metrics
src/cylinder.cpp     Cylinder solid mask usage, Zou-He inlet/outlet, run loop
script/plot.py       Poiseuille profile plotting
script/animate.py    2D velocity-field animation
```

## Build

```sh
cmake -S . -B build
cmake --build build
```

or

```sh
./script/compile.sh
```

## Run

```sh
./build/lbm_solver
```

Useful options:

```sh
./build/lbm_solver --case poiseuille --grid 32 --steps 20000 --nu 0.05 --boundary on-grid
./build/lbm_solver --case cylinder --grid 32 --steps 20000 --nu 0.05 --re 5
./build/lbm_solver --case cylinder --grid 32 --steps 20000 --nu 0.05 --re 40
```

OpenMP thread count can be selected at runtime:

```sh
./build/lbm_solver --threads 1
./build/lbm_solver --threads 4
```

For Poiseuille runs, the solver writes an averaged velocity profile CSV and
appends run-level metrics to `test/poiseuille_summary.csv`. For cylinder runs,
it writes the final 2D velocity field CSV.

To write 2D velocity-field snapshots and create an animation:

```sh
./build/lbm_solver --case cylinder --grid 32 --steps 20000 --nu 0.05 --re 40 --snapshot-interval 100
.venv/bin/python script/animate.py
```

For the cylinder case, the left boundary is a Zou-He velocity inlet, the right
boundary is a Zou-He pressure outlet, top and bottom are no-slip walls, and the
cylinder is no-slip bounce-back. The cylinder case is configured by Reynolds
number, and the solver computes the inlet velocity as
`inlet_ux = Re * nu / (2 * cylinder_radius)`.
By default, the cylinder radius is `0.125 * min(nx, ny)`. For example,
`--grid 32 --nu 0.05 --re 5` gives a `128x32` domain, radius `4`, and
`inlet_ux = 0.03125`.

## Verification Target

For Poiseuille flow, the code reports the relative L2 error between the
streamwise velocity averaged over `x` and the analytical solution

```text
u(y) = g y (H - y) / (2 nu),    tau = 0.5 + 3 nu
```

where the wall locations depend on the selected boundary condition. The default
is on-grid bounce-back. The previous mid-grid setup is still available with
`--boundary mid-grid`.

For Poiseuille flow, `--grid 32` is interpreted as `nx = 32`, `ny = 32`, so the
on-grid channel height is `H = ny - 1 = 31`. The default body force
`0.000390625` is inferred from a pressure drop with magnitude `0.0125` across
the streamwise length `L = nx = 32`.
