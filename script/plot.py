from __future__ import annotations
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scienceplots

plt.style.use(["science", "grid"])
plt.rcParams.update(
    {
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "grid.alpha": 0.25,
        "axes.labelsize": 8,
        "axes.titlesize": 10,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 6,
    }
)


REQUIRED_COLUMNS = {"y", "ux_average", "ux_analytical", "error"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot the LBM Poiseuille velocity profile.")
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        default=Path("test/profile.csv"),
        help="Input CSV written by lbm_solver.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("test/profile.png"),
        help="Output figure path.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Show the figure interactively after saving.",
    )
    return parser.parse_args()


def read_profile(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"input CSV not found: {path}")

    data = pd.read_csv(path)
    missing = REQUIRED_COLUMNS.difference(data.columns)
    if missing:
        missing_list = ", ".join(sorted(missing))
        raise ValueError(f"{path} is missing required columns: {missing_list}")

    return data.sort_values("y").reset_index(drop=True)


def compute_metrics(data: pd.DataFrame) -> tuple[float, float]:
    numerical = data["ux_average"].to_numpy(dtype=float)
    analytical = data["ux_analytical"].to_numpy(dtype=float)
    error = numerical - analytical

    denominator = np.sum(analytical**2)
    relative_l2 = (
        np.sqrt(np.sum(error**2) / denominator)
        if denominator > 0.0
        else np.sqrt(np.sum(error**2))
    )
    max_absolute = np.max(np.abs(error))
    return relative_l2, max_absolute


def plot_profile(data: pd.DataFrame, output_path: Path, input_path: Path) -> tuple[float, float]:
    y = data["y"].to_numpy(dtype=float)
    numerical = data["ux_average"].to_numpy(dtype=float)
    analytical = data["ux_analytical"].to_numpy(dtype=float)
    relative_l2, max_absolute = compute_metrics(data)

    fig, ax = plt.subplots(figsize=(3, 3), constrained_layout=True)
    ax.plot(numerical, y, "o", markersize=3, label="LBM")
    ax.plot(analytical, y, "-", label="Analytical")
    ax.set_xlabel(r"$u_x$")
    ax.set_ylabel(r"$y$")
    ax.set_title("Poiseuille Velocity Profile")
    ax.legend()
    ax.text(
        0.03,
        0.97,
        f"relative L2 = {relative_l2:.3e}\nmax abs = {max_absolute:.3e}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=7,
        bbox={"facecolor": "white", "edgecolor": "0.8", "alpha": 0.85},
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight")
    return relative_l2, max_absolute


def main() -> int:
    args = parse_args()
    data = read_profile(args.input)
    relative_l2, max_absolute = plot_profile(data, args.output, args.input)

    print(f"Input: {args.input}")
    print(f"Output: {args.output}")
    print(f"Relative_l2: {relative_l2:.6e}")
    print(f"max_absolute: {max_absolute:.6e}")

    if args.show:
        plt.show()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
