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
./build/lbm_solver --case poiseuille --grid 32x33 --steps 10000 --tau 0.65 --force-x 0.000390625 --boundary on-grid
./build/lbm_solver --case cylinder --grid 400x100 --steps 20000 --tau 0.6 --re 40 --cylinder-radius 10 --output test/cylinder_field.csv
./build/lbm_solver --case cylinder --grid 400x100 --steps 20000 --tau 0.6 --re 150 --cylinder-radius 10 --output test/cylinder_field.csv
```

OpenMP thread count can be selected at runtime:

```sh
./build/lbm_solver --threads 1
./build/lbm_solver --threads 4
```

For Poiseuille runs, the solver writes an averaged velocity profile CSV. For
cylinder runs, it writes the final 2D velocity field CSV.

To write 2D velocity-field snapshots and create an animation:

```sh
./build/lbm_solver --case cylinder --grid 400x100 --steps 5000 --tau 0.6 --re 40 --cylinder-radius 10 --snapshot-interval 100 --snapshot-dir test/snapshots
.venv/bin/python script/animate.py --input-dir test/snapshots --output test/velocity_field.gif
```

For the cylinder case, the left boundary is a Zou-He velocity inlet, the right
boundary is a Zou-He pressure outlet, top and bottom are no-slip walls, and the
cylinder is no-slip bounce-back. The cylinder case is configured by Reynolds
number, and the solver computes the inlet velocity as
`inlet_ux = Re * nu / (2 * cylinder_radius)`.
For example, `--grid 400x100 --tau 0.6 --cylinder-radius 10 --re 40`
gives `inlet_ux = 0.0666666666666667`.

## Verification Target

For Poiseuille flow, the code reports the relative L2 error between the
streamwise velocity averaged over `x` and the analytical solution

```text
u(y) = g y (H - y) / (2 nu),    nu = cs^2 (tau - 0.5)
```

where the wall locations depend on the selected boundary condition. The default
is on-grid bounce-back. The previous mid-grid setup is still available with
`--boundary mid-grid`.

The default `32x33` grid corresponds to an on-grid channel height
`H = ny - 1 = 32`. The default body force `0.000390625` is the pressure-gradient
equivalent of a pressure drop with magnitude `0.0125` across `L = 32`.
