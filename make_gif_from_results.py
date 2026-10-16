from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

DEFAULT_POINT_RADIUS_DATA = 0.025


def _sanitize_folder_name(raw_name: str) -> str:
    sanitized = "".join(ch if ch.isalnum() or ch in ("-", "_") else "_" for ch in raw_name).strip("_")
    return sanitized


def request_existing_results_folder(base_folder: str = "Results") -> Path:
    base_dir = Path(__file__).resolve().parent / base_folder
    base_dir.mkdir(parents=True, exist_ok=True)

    candidate: Path | None = None
    while candidate is None:
        raw_input = input(f"Inserisci il nome della sottocartella dentro {base_folder}/: ").strip()
        sanitized = _sanitize_folder_name(raw_input)
        if not sanitized:
            print("Il nome indicato non è valido. Include caratteri alfanumerici, '-' o '_'.")
            continue
        potential = base_dir / sanitized
        if not potential.exists():
            print(f"{potential} non esiste. Assicurati di averlo scritto correttamente.")
            continue
        candidate = potential

    print(f"Leggo i dati da: {candidate}")
    return candidate


def _build_levels(step_min: float, step_max: float, levels_count: int) -> tuple[np.ndarray, np.ndarray]:
    if np.isclose(step_min, step_max, rtol=1e-9, atol=1e-12):
        lower = np.nextafter(step_min, -np.inf)
        upper = np.nextafter(step_max, np.inf)
        if not np.isfinite(lower):
            lower = step_min - 1.0
        if not np.isfinite(upper):
            upper = step_max + 1.0
        if lower >= step_min:
            lower = step_min - max(abs(step_min) * 1e-6, 1e-9)
        if upper <= step_max:
            upper = step_max + max(abs(step_max) * 1e-6, 1e-9)
        if lower >= upper:
            upper = lower + max(abs(lower) * 1e-6, 1e-9)
        levels = np.array([lower, upper])
        ticks = np.array([step_min])
    else:
        levels = np.linspace(step_min, step_max, levels_count)
        levels = np.unique(levels)
        if len(levels) < 2:
            next_level = np.nextafter(step_min, np.inf)
            if not np.isfinite(next_level):
                next_level = step_min + 1.0
            levels = np.array([step_min, next_level])
        ticks_count = min(10, len(levels))
        ticks = np.linspace(levels[0], levels[-1], ticks_count)
    return levels, ticks


def create_gif(
    base_x: np.ndarray,
    base_y: np.ndarray,
    values_history: np.ndarray,
    disp_x_history: np.ndarray,
    disp_y_history: np.ndarray,
    times: np.ndarray,
    output_path: Path,
    title: str,
    x_limits: tuple[float, float] | None = None,
    y_limits: tuple[float, float] | None = None,
) -> None:
    if values_history.size == 0 or times.size == 0:
        raise ValueError("History arrays are empty.")
    levels_count = 100
    fig, ax = plt.subplots(figsize=(6, 6))
    plt.subplots_adjust(bottom=0.2)
    x_def = base_x + disp_x_history[0]
    y_def = base_y + disp_y_history[0]
    circles = [
        plt.Circle((x_def[i], y_def[i]), DEFAULT_POINT_RADIUS_DATA)
        for i in range(len(x_def))
    ]
    collection = PatchCollection(circles, cmap="viridis")
    collection.set_array(values_history[0])
    collection.set_edgecolor("none")
    ax.add_collection(collection)
    step_min = float(np.nanmin(values_history[0]))
    step_max = float(np.nanmax(values_history[0]))
    levels, ticks = _build_levels(step_min, step_max, levels_count)
    collection.set_clim(levels[0], levels[-1])
    cbar = fig.colorbar(collection, ax=ax, orientation="vertical")
    cbar.set_ticks(ticks)
    ax.set_title(f"{title} (t = {times[0]:.6f})")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")

    if x_limits is not None:
        ax.set_xlim(x_limits)
    else:
        ax.set_xlim(float(np.min(base_x)), float(np.max(base_x)))
    if y_limits is not None:
        ax.set_ylim(y_limits)
    else:
        ax.set_ylim(float(np.min(base_y)), float(np.max(base_y)))

    def animate(frame: int) -> None:
        x_def = base_x + disp_x_history[frame]
        y_def = base_y + disp_y_history[frame]
        for circle, x, y in zip(circles, x_def, y_def):
            circle.center = (x, y)
        step_values = values_history[frame]
        collection.set_array(step_values)
        step_min = float(np.nanmin(step_values))
        step_max = float(np.nanmax(step_values))
        levels, ticks = _build_levels(step_min, step_max, levels_count)
        collection.set_clim(levels[0], levels[-1])
        cbar.update_normal(collection)
        cbar.set_ticks(ticks)
        ax.set_title(f"{title} (t = {times[frame]:.6f})")

    anim = FuncAnimation(fig, animate, frames=len(times), repeat=False)
    writer = PillowWriter(fps=4)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    anim.save(output_path, writer=writer)
    plt.close(fig)
    print(f"In registrata la GIF in {output_path}")


def main() -> None:
    results_folder = request_existing_results_folder()
    history_path = results_folder / "field_histories.npz"
    if not history_path.exists():
        raise FileNotFoundError(f"{history_path} non trovato.")

    with np.load(history_path) as data:
        history_keys = list(data.keys())
        time_history = np.array(data["time_history"])
        gauss_x = np.array(data["gauss_x"])
        gauss_y = np.array(data["gauss_y"])
        disp_x_history = np.array(data["disp_x_history"])
        disp_y_history = np.array(data["disp_y_history"])

        variable_options = {
            "displacement_x": "disp_x_history",
            "displacement_y": "disp_y_history",
            "stress_xx": "stress_xx_history",
            "stress_xy": "stress_xy_history",
            "stress_yy": "stress_yy_history",
            "plastic_strain": "plastic_strain_history",
        }
        if "uniaxial_stress_history" in history_keys:
            variable_options["uniaxial_stress"] = "uniaxial_stress_history"

        print("Variabili disponibili:")
        for name in variable_options:
            print(f"  {name}")
        chosen = input(f"Digitare il nome della variabile da animare (default 'stress_xx'): ").strip()
        data_key = variable_options.get(chosen, "stress_xx_history")
        if data_key not in data:
            raise KeyError(f"{data_key} non presente nei dati.")

        values_history = np.array(data[data_key])

    if ag_x := gauss_x.size:
        base_x = np.array(gauss_x)
    else:
        raise ValueError("Nessun punto gaussiano disponibile.")
    base_y = np.array(gauss_y)

    auto_x_min = float(np.min(base_x + np.min(disp_x_history, axis=0)))
    auto_x_max = float(np.max(base_x + np.max(disp_x_history, axis=0)))
    auto_y_min = float(np.min(base_y + np.min(disp_y_history, axis=0)))
    auto_y_max = float(np.max(base_y + np.max(disp_y_history, axis=0)))
    def _zoom_limits(min_val: float, max_val: float) -> tuple[float, float]:
        span = max_val - min_val
        margin = span * 0.02 if span > 0 else 0.1
        return (min_val - margin, max_val + margin)
    x_limits = _zoom_limits(auto_x_min, auto_x_max)
    y_limits = _zoom_limits(auto_y_min, auto_y_max)

    if values_history.ndim == 1:
        values_history = values_history.reshape(values_history.shape[0], 1)

    min_length = min(
        values_history.shape[0],
        disp_x_history.shape[0],
        disp_y_history.shape[0],
        time_history.shape[0],
    )
    if min_length == 0:
        raise ValueError("Nessuno step disponibile per la GIF.")

    output_path = results_folder / f"{data_key}_animation.gif"
    create_gif(
        base_x=np.array(base_x),
        base_y=np.array(base_y),
        values_history=values_history[:min_length],
        disp_x_history=disp_x_history[:min_length],
        disp_y_history=disp_y_history[:min_length],
        times=time_history[:min_length],
        output_path=output_path,
        title=f"Animazione {chosen or 'stress_xx'}",
        x_limits=x_limits,
        y_limits=y_limits,
    )


if __name__ == "__main__":
    main()
