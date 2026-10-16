from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

RESULTS_FOLDERS = [
    # "GAP_p1_80x20",
    # "GAP_p2_160x40",
    "GAP_p2_100x100"
]

GID_REACTION_FILES = [
    "REACTION_size_0.5.txt"#,
    # "REACTION_size_1.txt",
    # "REACTION_size_2.txt",
]

COLORS = [
    "tab:blue",
    "tab:orange",
    "tab:grey",
    "tab:red",
    "tab:purple",
    "tab:brown",
]


def read_reaction_file(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    times = data[:, 0]
    reaction_x = abs(data[:, 1])
    return times, reaction_x


def load_results_history(results_folder: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    history_path = results_folder / "field_histories.npz"
    if not history_path.exists():
        raise FileNotFoundError(f"{history_path} non trovato.")
    with np.load(history_path) as data:
        time_history = data["time_history"]
        reaction_history = data["reaction_history"]
        imposed_history = data.get("imposed_displacement_history")
    return time_history, reaction_history, imposed_history


def get_stored_reaction_curve(results_folder: Path) -> tuple[np.ndarray, np.ndarray]:
    time_history, reaction_history, imposed_history = load_results_history(results_folder)
    if imposed_history is None or imposed_history.size != reaction_history.size:
        imposed = 250.0 * time_history
    else:
        imposed = imposed_history
    return imposed, reaction_history


def main() -> None:
    results_to_plot = RESULTS_FOLDERS
    gid_files_to_plot = GID_REACTION_FILES

    base_dir = Path(__file__).resolve().parent
    gid_dir = base_dir / "GID_case" / "Reaction_results"

    fig, ax = plt.subplots(figsize=(8, 5))
    color_idx = 0

    for folder in results_to_plot:
        results_dir = base_dir / "Results" / folder
        try:
            stored_imposed, stored_reaction = get_stored_reaction_curve(results_dir)
        except FileNotFoundError as err:
            print(err)
            continue
        color = COLORS[color_idx % len(COLORS)]
        color_idx += 1
        ax.plot(
            stored_imposed,
            stored_reaction,
            marker="o",
            markersize = 2,
            linestyle="-",
            linewidth=1.5,
            color=color,
            # label=f"Results/{folder}",
            label=f"Gap-SBM_size_1.8"
        )

    for filename in gid_files_to_plot:
        reaction_path = gid_dir / filename
        try:
            times, reaction_x = read_reaction_file(reaction_path)
        except FileNotFoundError:
            print(f"{reaction_path} non trovato.")
            continue
        imposed = 250.0 * times
        color = COLORS[color_idx % len(COLORS)]
        color_idx += 1
        ax.plot(
            imposed,
            reaction_x,
            marker="^",
            markersize = 2,
            linestyle="--",
            linewidth=1.2,
            color=color,
            # label=f"GID_case/{filename}",
            label=f"FEM_size_0.5",
        )

    ax.set_xlabel("Imposed displacement")
    ax.set_ylabel("Reaction X")
    ax.set_title("Reaction X vs imposed displacement")
    ax.grid(True, linestyle="--", alpha=0.6)
    ax.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
