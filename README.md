# LBM Final Project

This project contains a D2Q9 lattice Boltzmann method solver with a BGK
collision operator and optional OpenMP parallelism. It currently supports:

- plane Poiseuille flow for correctness verification
- flow past a circular cylinder for 2D flow visualization
- flow past a symmetric airfoil for 2D flow visualization
- on-grid and mid-grid wall bounce-back for Poiseuille comparison
- Zou-He velocity inlet and density outlet for obstacle-flow cases
- CSV outputs for plotting, snapshots, timing, and MLUPS comparison

## Structure

```text
include/lbm.hpp      Shared D2Q9 constants, config, data types, declarations
src/main.cpp         Command-line parsing, validation, and case dispatch
src/lbm.cpp          Common lattice utilities, initialization, output helpers
src/poiseuille.cpp   Poiseuille update loop, analytical profile, error metrics
src/cylinder.cpp     Cylinder update loop, Zou-He inlet/outlet, run metrics
src/airfoil.cpp      Airfoil update loop, Zou-He inlet/outlet, run metrics
script/compile.sh    CMake configure/build helper
script/plot.py       Poiseuille profile plot
script/animate.py    2D velocity-field GIF animation
```

Generated files are written under `test/`.

## Build

```sh
./script/compile.sh
```

or manually:

```sh
cmake -S . -B build
cmake --build build
```

## Poiseuille Flow

The Poiseuille case is driven by a constant body force equivalent to a pressure
gradient. The streamwise direction is periodic, and the top/bottom walls use
no-slip bounce-back.

Default reference setup:

```text
Pressure drop = 0.0125
Channel height H = 32
Viscosity nu = 0.05
tau = 0.5 + 3 nu = 0.65
```

For Poiseuille, `--grid N` means physical channel height `H = N`.

```text
on-grid:  ny = H + 1
mid-grid: ny = H + 2
```

Example:

```sh
./build/lbm_solver --case poiseuille --grid 32 --steps 20000 --nu 0.05 --boundary on-grid
./build/lbm_solver --case poiseuille --grid 32 --steps 20000 --nu 0.05 --boundary mid-grid
```

Outputs:

```text
test/poiseuille_profile.csv
test/poiseuille_summary.csv
test/poiseuille_snapshots/field_*.csv
test/poiseuille_snapshots/profile_*.csv
```

The profile CSV contains:

```text
y, ux_average, ux_analytical, error
```

The summary CSV appends run-level metrics including:

```text
nx, ny, nu, tau, force_x, boundary, steps, threads, computing_time, mlups, relative_l2
```

The analytical solution used for verification is:

```text
u(y) = g y (H - y) / (2 nu)
```

where `g` is the applied body force.

## Cylinder Flow

The cylinder case uses:

- no-slip bounce-back on the top/bottom walls
- no-slip bounce-back on the circular cylinder
- Zou-He velocity boundary condition at the inlet
- Zou-He density boundary condition at the outlet

For cylinder flow, `--grid N` creates a rectangular domain:

```text
nx = 4N
ny = N
```

The default cylinder location and radius are:

```text
center = (nx / 4, ny / 2)
radius = round(0.125 * min(nx, ny))
```

The cylinder case is configured by Reynolds number and a safe inlet velocity.
The solver infers viscosity from:

```text
nu = inlet_ux * D / Re
tau = 0.5 + 3 nu
```

where `D = 2 * radius`.

Current safety validation:

```text
0.01 <= inlet_ux <= 0.08
0.56 <= tau <= 1.0
```

This means higher Reynolds number cases require larger grids. Recommended
starting cases:

```text
Re = 5    grid = 32
Re = 40   grid = 64
Re = 150  grid = 256
```

Examples:

```sh
./build/lbm_solver --case cylinder --grid 32 --steps 20000 --re 5 --inlet-ux 0.05
./build/lbm_solver --case cylinder --grid 64 --steps 40000 --re 40 --inlet-ux 0.05
./build/lbm_solver --case cylinder --grid 256 --steps 80000 --re 150 --inlet-ux 0.05
```

Do not pass `--nu` or `--tau` for the cylinder case. The program rejects them
because cylinder viscosity is inferred from `Re` and `inlet_ux`.

Outputs:

```text
test/cylinder_field.csv
test/cylinder_snapshots/field_*.csv
```

The field CSV contains:

```text
x, y, rho, ux, uy, velocity, solid
```

## Airfoil Flow

The airfoil case uses the same inlet, outlet, and wall treatment as the
cylinder case, but replaces the circular cylinder mask with a symmetric
NACA 0012-style airfoil mask.

For airfoil flow, `--grid N` creates:

```text
nx = 4N
ny = N
```

The default airfoil geometry is:

```text
leading edge = (nx / 4, ny / 2)
chord = round(0.5 * min(nx, ny))
thickness ratio = 0.12
```

The Reynolds number is based on chord length:

```text
nu = inlet_ux * chord / Re
tau = 0.5 + 3 nu
```

Example:

```sh
./build/lbm_solver --case airfoil --grid 64 --steps 40000 --re 40 --inlet-ux 0.05
```

Outputs:

```text
test/airfoil_field.csv
test/airfoil_snapshots/field_*.csv
```

## OpenMP

Use `--threads` to control the OpenMP thread count:

```sh
./build/lbm_solver --case poiseuille --grid 64 --steps 80000 --threads 1
./build/lbm_solver --case poiseuille --grid 64 --steps 80000 --threads 4
```

Use MLUPS from the program output or `test/poiseuille_summary.csv` to compare
parallel performance.

## Plotting

Create the Poiseuille profile plot:

```sh
.venv/bin/python script/plot.py
```

Input:

```text
test/poiseuille_profile.csv
test/poiseuille_summary.csv
```

Output:

```text
test/poiseuille_profile.png
```

## Animation

Generate field snapshots by running with a positive snapshot interval:

```sh
./build/lbm_solver --case poiseuille --grid 32 --steps 20000 --snapshot-interval 100
./build/lbm_solver --case cylinder --grid 64 --steps 40000 --re 40 --inlet-ux 0.05 --snapshot-interval 100
./build/lbm_solver --case airfoil --grid 64 --steps 40000 --re 40 --inlet-ux 0.05 --snapshot-interval 100
```

Then create the GIF:

```sh
.venv/bin/python script/animate.py
```

In `script/animate.py`, choose the case by editing:

```python
ANIMATE_CASE = "poiseuille"  # "poiseuille", "cylinder", "airfoil", or "all"
```

Outputs:

```text
test/poiseuille_velocity_field.gif
test/cylinder_velocity_field.gif
test/airfoil_velocity_field.gif
```

`FRAME_INTERVAL` controls how many snapshots are used in the animation. For
example, if snapshots are written every 100 steps and `FRAME_INTERVAL = 500`,
the GIF uses every fifth snapshot.

## VS Code

The `.vscode` folder contains launch/task configurations for:

- `Run LBM: Poiseuille`
- `Run LBM: Cylinder`
- `Run LBM: Airfoil`

The Poiseuille configuration asks for grid size, steps, boundary condition, and
threads. The cylinder and airfoil configurations ask for grid size, steps,
Reynolds number, snapshot interval, and threads. Obstacle-flow viscosity is
inferred automatically.
