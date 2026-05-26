from __future__ import annotations
import argparse
import re
from pathlib import Path
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scienceplots

plt.style.use(["science", "grid"])
plt.rcParams.update(
    {
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "grid.alpha": 0,
        "axes.labelsize": 8,
        "axes.titlesize": 10,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 6,
    }
)

REQUIRED_COLUMNS = {"x", "y", "rho", "ux", "uy", "speed", "solid"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Animate 2D LBM velocity-field snapshots."
    )
    parser.add_argument(
        "-i",
        "--input-dir",
        type=Path,
        default=Path("test/snapshots"),
        help="Directory containing field_*.csv snapshots.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("test/velocity_field.gif"),
        help="Output animation path. GIF is the simplest option.",
    )
    parser.add_argument(
        "-q",
        "--quantity",
        choices=["speed", "ux", "uy", "rho"],
        default="speed",
        help="Field quantity to render as a heatmap.",
    )
    parser.add_argument(
        "--fps", type=int, default=12, help="Frames per second for saved animation."
    )
    parser.add_argument(
        "--quiver-stride",
        type=int,
        default=4,
        help="Arrow spacing for velocity vectors. Use 0 to disable arrows.",
    )
    parser.add_argument(
        "--show", action="store_true", help="Show the animation interactively."
    )
    return parser.parse_args()


def step_from_path(path: Path) -> int:
    match = re.search(r"(\d+)", path.stem)
    return int(match.group(1)) if match else 0


def read_snapshots(input_dir: Path) -> list[tuple[int, pd.DataFrame]]:
    paths = sorted(input_dir.glob("field_*.csv"), key=step_from_path)
    if not paths:
        raise FileNotFoundError(
            f"no field_*.csv snapshots found in {input_dir}. "
            "Generate them with: ./build/lbm_solver --snapshot-interval 100"
        )

    snapshots = []
    for path in paths:
        data = pd.read_csv(path)
        missing = REQUIRED_COLUMNS.difference(data.columns)
        if missing:
            missing_list = ", ".join(sorted(missing))
            raise ValueError(f"{path} is missing required columns: {missing_list}")
        snapshots.append(
            (step_from_path(path), data.sort_values(["y", "x"]).reset_index(drop=True))
        )
    return snapshots


def field_grid(data: pd.DataFrame, column: str) -> np.ndarray:
    return (
        data.pivot(index="y", columns="x", values=column)
        .sort_index(ascending=True)
        .to_numpy()
    )


def coordinate_grids(data: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    x_values = np.sort(data["x"].unique())
    y_values = np.sort(data["y"].unique())
    return np.meshgrid(x_values, y_values)


def build_animation(
    snapshots: list[tuple[int, pd.DataFrame]],
    output_path: Path,
    quantity: str,
    fps: int,
    quiver_stride: int,
    show: bool,
) -> None:
    first = snapshots[0][1]
    first_field = field_grid(first, quantity)
    x_grid, y_grid = coordinate_grids(first)

    vmin = min(field_grid(data, quantity).min() for _, data in snapshots)
    vmax = max(field_grid(data, quantity).max() for _, data in snapshots)
    if np.isclose(vmin, vmax):
        vmax = vmin + 1.0

    fig, ax = plt.subplots(figsize=(4, 4), constrained_layout=True)
    image = ax.imshow(
        first_field,
        origin="lower",
        cmap="coolwarm",
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
        extent=[
            first["x"].min() - 0.5,
            first["x"].max() + 0.5,
            first["y"].min() - 0.5,
            first["y"].max() + 0.5,
        ],
        aspect="auto",
    )
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label("Velocity")
    """
    quiver = None
    if quiver_stride > 0:
        stride = max(1, quiver_stride)
        ux = field_grid(first, "ux")
        uy = field_grid(first, "uy")
        quiver = ax.quiver(
            x_grid[::stride, ::stride],
            y_grid[::stride, ::stride],
            ux[::stride, ::stride],
            uy[::stride, ::stride],
            color="white",
            alpha=0.75,
            scale_units="xy",
            scale=None,
            width=0.003,
        )
    """
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")

    def update(frame: int):
        step, data = snapshots[frame]
        image.set_data(field_grid(data, quantity))
        artists: list[object] = [image]
        """
        if quiver is not None:
            ux = field_grid(data, "ux")
            uy = field_grid(data, "uy")
            stride = max(1, quiver_stride)
            quiver.set_UVC(ux[::stride, ::stride], uy[::stride, ::stride])
            artists.append(quiver)
        """
        ax.set_title(f"2D Velocity Field, Step {step}")
        return artists

    anim = animation.FuncAnimation(
        fig, update, frames=len(snapshots), interval=100, blit=False
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    anim.save(output_path, writer=animation.PillowWriter(fps=fps))

    if show:
        plt.show()
    plt.close(fig)


def main() -> int:
    args = parse_args()
    snapshots = read_snapshots(args.input_dir)
    build_animation(
        snapshots=snapshots,
        output_path=args.output,
        quantity=args.quantity,
        fps=args.fps,
        quiver_stride=args.quiver_stride,
        show=args.show,
    )
    print(f"Frames: {len(snapshots)}")
    print(f"Quantity: {args.quantity}")
    print(f"Output: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
