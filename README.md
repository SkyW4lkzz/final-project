# LBM Final Project

This project currently contains a single-thread D2Q9 lattice Boltzmann solver
with a BGK collision operator. The implemented verification case is plane
Poiseuille flow driven by a constant body force:

- periodic boundary condition in the streamwise direction
- bounce-back wall boundary condition at the top and bottom walls
- Guo forcing term for the pressure-gradient/body-force source
- numerical profile compared with the analytical parabolic velocity profile
- runtime reported in MLUPS for later performance comparison

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
./build/lbm_solver --grid 40x33 --steps 10000 --tau 1.0 --force-x 0.00125 --boundary on-grid
./build/lbm_solver --grid 128x64 --steps 20000 --tau 0.8 --force-x 1e-6
```

The solver writes the averaged velocity profile to
`test/poiseuille_profile.csv` by default.

To write 2D velocity-field snapshots and create an animation:

```sh
./build/lbm_solver --grid 64 --steps 5000 --snapshot-interval 100 --snapshot-dir test/field_snapshots
.venv/bin/python script/animate.py --input-dir test/field_snapshots --output test/velocity_field.gif
```

## Verification Target

For Poiseuille flow, the code reports the relative L2 error between the
streamwise velocity averaged over `x` and the analytical solution

```text
u(y) = g y (H - y) / (2 nu),    nu = cs^2 (tau - 0.5)
```

where the wall locations depend on the selected boundary condition. The default
is on-grid bounce-back. The previous mid-grid setup is still available with
`--boundary mid-grid`.

The default `40x33` grid corresponds to a streamwise length `L = 40` and an
on-grid channel height `H = ny - 1 = 32`. The default body force
`0.00125` is the pressure-gradient equivalent of a pressure drop with magnitude
`0.05` across `L = 40`.
