from __future__ import annotations

import re
from pathlib import Path

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scienceplots
from mpl_toolkits.axes_grid1 import make_axes_locatable
from tqdm import tqdm

plt.style.use(["science", "grid"])
plt.rcParams.update(
    {
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "grid.alpha": 0.,
        "axes.labelsize": 8,
        "axes.titlesize": 10,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 6,
    }
)

ANIMATE_CASE = "cylinder"  # Choose "all", "poiseuille", or "cylinder"

FPS = 10
FRAME_INTERVAL = 100

FIGURE_HEIGHT = 3.5
MAX_FIGURE_WIDTH = 10.5
POISEUILLE_DISPLAY_RATIO = 4.0

REQUIRED_COLUMNS = {"x", "y", "rho", "ux", "uy", "velocity", "solid"}
CASE_PATHS = {
    "poiseuille": (
        Path("test/poiseuille_snapshots"),
        Path("test/poiseuille_velocity_field.gif"),
    ),
    "cylinder": (
        Path("test/cylinder_snapshots"),
        Path("test/cylinder_velocity_field.gif"),
    ),
}


def step_from_path(path: Path) -> int:
    match = re.search(r"(\d+)", path.stem)
    return int(match.group(1)) if match else 0


def read_snapshots(input_dir: Path, case_name: str) -> list[tuple[int, pd.DataFrame]]:
    paths = sorted(input_dir.glob("field_*.csv"), key=step_from_path)
    paths = [path for path in paths if step_from_path(path) % FRAME_INTERVAL == 0]
    if not paths:
        raise FileNotFoundError(
            f"No field_*.csv snapshots found in {input_dir}. "
            f"Generate them with: ./build/lbm_solver --case {case_name} --snapshot-interval 100"
        )

    snapshots = []
    for path in paths:
        data = pd.read_csv(path)
        if "velocity" not in data.columns and "speed" in data.columns:
            data = data.rename(columns={"speed": "velocity"})

        missing = REQUIRED_COLUMNS.difference(data.columns)
        if missing:
            missing_list = ", ".join(sorted(missing))
            raise ValueError(f"{path} is missing required columns: {missing_list}")

        snapshots.append(
            (step_from_path(path), data.sort_values(["y", "x"]).reset_index(drop=True))
        )
    return snapshots


def field_grid(data: pd.DataFrame) -> np.ndarray:
    return (
        data.pivot(index="y", columns="x", values="velocity")
        .sort_index(ascending=True)
        .to_numpy()
    )


def solid_grid(data: pd.DataFrame) -> np.ndarray:
    return (
        data.pivot(index="y", columns="x", values="solid")
        .sort_index(ascending=True)
        .to_numpy()
    )


def build_animation(
    snapshots: list[tuple[int, pd.DataFrame]], output_path: Path, case_name: str
) -> None:
    first = snapshots[0][1]
    first_field = field_grid(first)
    first_solid = solid_grid(first)

    vmin = min(field_grid(data).min() for _, data in snapshots)
    vmax = max(field_grid(data).max() for _, data in snapshots)
    if np.isclose(vmin, vmax):
        vmax = vmin + 1.0

    ny, nx = first_field.shape
    display_ratio = POISEUILLE_DISPLAY_RATIO if case_name == "poiseuille" else nx / ny
    image_aspect = "auto" if case_name == "poiseuille" else "equal"
    figure_width = min(MAX_FIGURE_WIDTH, FIGURE_HEIGHT * display_ratio)
    fig, ax = plt.subplots(figsize=(figure_width, FIGURE_HEIGHT))
    image = ax.imshow(
        first_field,
        origin="lower",
        cmap="inferno",
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
        aspect=image_aspect,
    )
    solid_overlay = ax.imshow(
        np.ma.masked_where(first_solid == 0, first_solid),
        origin="lower",
        cmap="gray_r",
        interpolation="nearest",
        aspect=image_aspect,
        vmin=0,
        vmax=1,
    )
    divider = make_axes_locatable(ax)
    colorbar_axis = divider.append_axes("bottom", size="5%", pad=0.35)
    cbar = fig.colorbar(image, cax=colorbar_axis, orientation="horizontal")
    cbar.set_label("Velocity")
    colorbar_axis.xaxis.set_ticks_position("bottom")
    colorbar_axis.xaxis.set_label_position("bottom")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    ax.set_title(f"{case_name.title()} Velocity Field, Step {snapshots[0][0]}", pad=8)
    fig.tight_layout()

    def update(frame: int):
        step, data = snapshots[frame]
        image.set_data(field_grid(data))
        solid = solid_grid(data)
        solid_overlay.set_data(np.ma.masked_where(solid == 0, solid))
        ax.set_title(f"{case_name.title()} Velocity Field, Step {step}", pad=8)
        return [image, solid_overlay]

    anim = animation.FuncAnimation(
        fig, update, frames=len(snapshots), interval=100, blit=False
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)

    def update_progress(frame: int, total: int) -> None:
        progress.update(1)

    with tqdm(total=len(snapshots), desc=case_name.title(), unit=" frames") as progress:
        anim.save(
            output_path,
            writer=animation.PillowWriter(fps=FPS),
            progress_callback=update_progress,
        )
    plt.close(fig)


def selected_cases() -> dict[str, tuple[Path, Path]]:
    if ANIMATE_CASE == "all":
        return CASE_PATHS
    if ANIMATE_CASE in CASE_PATHS:
        return {ANIMATE_CASE: CASE_PATHS[ANIMATE_CASE]}
    raise ValueError('ANIMATE_CASE must be "all", "poiseuille", or "cylinder".')


def main() -> int:
    animated_cases = 0
    for case_name, (snapshot_dir, output_path) in selected_cases().items():
        if not list(snapshot_dir.glob("field_*.csv")):
            print(f"Skip {case_name}: no snapshots found in {snapshot_dir}")
            continue

        snapshots = read_snapshots(snapshot_dir, case_name)
        build_animation(snapshots, output_path, case_name)
        animated_cases += 1
        print(f"Case: {case_name}")
        print(f"Frames: {len(snapshots)}")
        print(f"Output: {output_path}")

    if animated_cases == 0:
        raise FileNotFoundError(
            "No snapshots found in test/poiseuille_snapshots or test/cylinder_snapshots."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
