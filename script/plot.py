from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scienceplots

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

REQUIRED_COLUMNS = {"y", "ux_average", "ux_analytical", "error"}
PROFILE_PATH = Path("test/poiseuille_profile.csv")
SUMMARY_PATH = Path("test/poiseuille_summary.csv")
OUTPUT_PATH = Path("test/poiseuille_profile.png")


def read_profile(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Input CSV not found: {path}")

    data = pd.read_csv(path)
    missing = REQUIRED_COLUMNS.difference(data.columns)
    if missing:
        missing_list = ", ".join(sorted(missing))
        if {"x", "rho", "ux", "uy", "velocity", "solid"}.issubset(data.columns):
            raise ValueError(
                f"{path} is a 2D field CSV, not a Poiseuille profile CSV. "
                "Run the poiseuille case or pass test/poiseuille_profile.csv with --input."
            )
        raise ValueError(f"{path} is missing required columns: {missing_list}")

    return data.sort_values("y").reset_index(drop=True)


def read_relative_l2(path: Path) -> float | None:
    if not path.exists():
        return None

    data = pd.read_csv(path)
    if data.empty or "relative_l2" not in data.columns:
        return None
    return float(data.iloc[-1]["relative_l2"])


def plot_profile(
    data: pd.DataFrame, output_path: Path, relative_l2: float | None
) -> None:
    y = data["y"].to_numpy(dtype=float)
    numerical = data["ux_average"].to_numpy(dtype=float)
    analytical = data["ux_analytical"].to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(4, 4), constrained_layout=True)
    ax.plot(numerical, y, "o", markersize=3, label="LBM")
    ax.plot(analytical, y, "-", label="Analytical")
    ax.set_xlabel(r"$u_x$")
    ax.set_ylabel(r"$y$")
    title = "Poiseuille Velocity Profile"
    if relative_l2 is not None:
        title += f", Error = {relative_l2:.3e}"
    ax.set_title(title)
    ax.legend()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight")


def main() -> int:
    data = read_profile(PROFILE_PATH)
    relative_l2 = read_relative_l2(SUMMARY_PATH)
    plot_profile(data, OUTPUT_PATH, relative_l2)
    print(f"Input : {PROFILE_PATH}")
    print(f"Output: {OUTPUT_PATH}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
