from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DEFAULT_FOLDER_NAMES = [
    "GAP_p2_160x40",
    "GAP_p2_320x80",
    "GAP_p2_100x100",
    "GAP_p1_80x20",
    "GAP_p1_320x80",
    
    # "GAP_p1_80x100"
]
DEFAULT_FOLDER_COLORS = [
    "tab:orange",
    "tab:blue",
    "tab:green",
    "tab:red"]


def load_results_folder(folder_name: str) -> tuple[np.ndarray, np.ndarray]:
    results_dir = Path(__file__).resolve().parent / "Results" / folder_name
    history_path = results_dir / "field_histories.npz"
    if not history_path.exists():
        raise FileNotFoundError(f"{history_path} non trovato per la cartella {folder_name}.")
    with np.load(history_path) as data:
        time_history = data["time_history"]
        reaction_history = data["reaction_history"]
        imposed_history = data.get("imposed_displacement_history")
    if imposed_history is None or imposed_history.size != reaction_history.size:
        imposed = 250.0 * time_history
    else:
        imposed = imposed_history
    return imposed, reaction_history


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot Reaction vs imposed for multiple Results folders.")
    parser.add_argument(
        "folders",
        nargs="*",
        help="Nome delle sottocartelle dentro Results/ da plottare.",
    )
    args = parser.parse_args()

    if args.folders:
        folder_names = args.folders
    else:
        folder_names = DEFAULT_FOLDER_NAMES

    fig, ax = plt.subplots(figsize=(8, 5))
    for idx, folder in enumerate(folder_names):
        try:
            imposed, reaction = load_results_folder(folder)
        except FileNotFoundError as err:
            print(err)
            continue
        color = DEFAULT_FOLDER_COLORS[idx % len(DEFAULT_FOLDER_COLORS)]
        ax.plot(imposed, reaction, label=folder, marker="o", linewidth=1.5, color=color)

    ax.set_xlabel("Imposed displacement")
    ax.set_ylabel("Reaction force")
    ax.set_title("Reaction vs imposed (Results)")
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
