from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.widgets import Slider


def _sanitize_folder_name(raw_name: str) -> str:
    sanitized = "".join(ch if ch.isalnum() or ch in ("-", "_") else "_" for ch in raw_name).strip("_")
    return sanitized


def request_existing_results_folder(base_folder: str = "Results") -> Path:
    """
    Prompt the user for the subfolder stored under WORKSPACE/Results and return it.
    """
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

    print(f"Ricarico i dati da: {candidate}")
    return candidate


def build_filtered_triangulation(x_coord: Iterable[float], y_coord: Iterable[float]) -> mtri.Triangulation:
    triangulation = mtri.Triangulation(np.array(x_coord), np.array(y_coord))
    filtered = []
    for triangle in triangulation.triangles:
        x1, y1 = triangulation.x[triangle[0]], triangulation.y[triangle[0]]
        x2, y2 = triangulation.x[triangle[1]], triangulation.y[triangle[1]]
        x3, y3 = triangulation.x[triangle[2]], triangulation.y[triangle[2]]

        # Lunghezze dei lati del triangolo
        edge_lengths = (
            np.hypot(x1 - x2, y1 - y2),
            np.hypot(x2 - x3, y2 - y3),
            np.hypot(x3 - x1, y3 - y1),
        )
        if max(edge_lengths) <= 2:
            filtered.append(triangle)
    if not filtered:
        # Evita triangolazioni vuote: usa quella completa se il filtro rimuove tutto
        filtered = triangulation.triangles
    return mtri.Triangulation(triangulation.x, triangulation.y, filtered)


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


def create_time_slider_figure(triangulation: mtri.Triangulation, values_history: np.ndarray, times: np.ndarray, title: str) -> None:
    if values_history.size == 0 or times.size == 0:
        return

    n_steps = values_history.shape[0]
    steps_array = [np.asarray(values_history[step_index]) for step_index in range(n_steps)]
    levels_count = 100
    global_min = float(np.nanmin(values_history))
    global_max = float(np.nanmax(values_history))
    levels, ticks = _build_levels(global_min, global_max, levels_count)
    cmap = "viridis"

    fig, ax = plt.subplots(figsize=(6, 6))
    plt.subplots_adjust(bottom=0.2)
    contour = ax.tricontourf(triangulation, steps_array[0], levels=levels, cmap=cmap)
    cbar = fig.colorbar(contour, ax=ax, orientation="vertical")
    cbar.set_ticks(ticks)
    ax.set_title(f"{title} (t = {times[0]:.6f})")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")

    slider_ax = fig.add_axes([0.15, 0.08, 0.7, 0.03])
    step_slider = Slider(slider_ax, "Step", 0, n_steps - 1, valinit=0, valstep=1, valfmt="%0.0f")

    def update(_val: float) -> None:
        nonlocal contour
        index = int(step_slider.val)
        for coll in contour.collections:
            coll.remove()
        contour = ax.tricontourf(triangulation, steps_array[index], levels=levels, cmap=cmap)
        cbar.update_normal(contour)
        ax.set_title(f"{title} (t = {times[index]:.6f})")
        fig.canvas.draw_idle()

    step_slider.on_changed(update)
    fig._step_slider = step_slider
    fig._colorbar_ref = cbar


def create_deformed_time_slider_figure(
    triangulation: mtri.Triangulation,
    values_history: np.ndarray,
    disp_x_history: np.ndarray,
    disp_y_history: np.ndarray,
    times: np.ndarray,
    title: str,
) -> None:
    if (
        values_history.size == 0
        or disp_x_history.size == 0
        or disp_y_history.size == 0
        or times.size == 0
    ):
        return

    min_length = min(
        values_history.shape[0],
        disp_x_history.shape[0],
        disp_y_history.shape[0],
        times.shape[0],
    )
    if min_length == 0:
        return

    values_history = values_history[:min_length]
    disp_x_history = disp_x_history[:min_length]
    disp_y_history = disp_y_history[:min_length]
    times = times[:min_length]

    base_x = np.array(triangulation.x)
    base_y = np.array(triangulation.y)

    steps_array = [np.asarray(values_history[step]) for step in range(min_length)]
    levels_count = 100
    initial_values = steps_array[0]
    step_min = float(np.nanmin(initial_values))
    step_max = float(np.nanmax(initial_values))
    levels, ticks = _build_levels(step_min, step_max, levels_count)

    fig, ax = plt.subplots(figsize=(6, 6))
    plt.subplots_adjust(bottom=0.2)
    x_def = base_x + np.asarray(disp_x_history[0])
    y_def = base_y + np.asarray(disp_y_history[0])
    global_min = float(np.nanmin(values_history))
    global_max = float(np.nanmax(values_history))
    levels, ticks = _build_levels(global_min, global_max, levels_count)
    cmap = "viridis"
    contour = ax.tricontourf(
        mtri.Triangulation(x_def, y_def, triangulation.triangles),
        initial_values,
        levels=levels,
        cmap=cmap,
    )

    cbar = fig.colorbar(contour, ax=ax, orientation="vertical")
    cbar.set_ticks(ticks)
    ax.set_title(f"{title} (t = {times[0]:.6f})")
    ax.set_aspect("equal", adjustable="datalim")

    slider_ax = fig.add_axes([0.15, 0.08, 0.7, 0.03])
    step_slider = Slider(slider_ax, "Step", 0, min_length - 1, valinit=0, valstep=1, valfmt="%0.0f")

    def update(_val: float) -> None:
        nonlocal contour
        index = int(step_slider.val)
        x_def = base_x + np.asarray(disp_x_history[index])
        y_def = base_y + np.asarray(disp_y_history[index])
        step_values = steps_array[index]
        for coll in contour.collections:
            coll.remove()
        contour = ax.tricontourf(
            mtri.Triangulation(x_def, y_def, triangulation.triangles),
            step_values,
            levels=levels,
            cmap=cmap,
        )
        cbar.update_normal(contour)
        ax.set_title(f"{title} (t = {times[index]:.6f})")
        fig.canvas.draw_idle()

    step_slider.on_changed(update)
    fig._step_slider = step_slider
    fig._colorbar_ref = cbar

def print_summary(results_folder: Path) -> None:
    summary_path = results_folder / "summary.json"
    if not summary_path.exists():
        return
    with open(summary_path, "r", encoding="utf-8") as summary_file:
        summary = json.load(summary_file)
    print("\nRiepilogo dei risultati salvati:")
    for key, value in summary.items():
        print(f"  {key}: {value}")


def main() -> None:
    results_folder = request_existing_results_folder()
    print_summary(results_folder)

    history_path = results_folder / "field_histories.npz"
    if not history_path.exists():
        raise FileNotFoundError(f"{history_path} non trovato.")

    with np.load(history_path) as data:
        time_history = data["time_history"]
        gauss_x = data["gauss_x"]
        gauss_y = data["gauss_y"]
        disp_x_history = data["disp_x_history"]
        disp_y_history = data["disp_y_history"]
        stress_xx_history = data["stress_xx_history"]
        stress_xy_history = data["stress_xy_history"]
        stress_yy_history = data["stress_yy_history"]
        uniaxial_stress_history = data["uniaxial_stress_history"]
        plastic_strain_history = data["plastic_strain_history"]
        reaction_history = data["reaction_history"]
        imposed_history = data.get("imposed_displacement_history")

    triangulation = build_filtered_triangulation(gauss_x, gauss_y)
    response = input("Vuoi visualizzare le mappe sulla configurazione deformata? [s/N]: ").strip().lower()
    use_deformed = response.startswith("s")

    if use_deformed and (
        disp_x_history.size == 0 or disp_y_history.size == 0
    ):
        print("I dati di spostamento non sono sufficienti; mostro la configurazione originale.")
        use_deformed = False

    if use_deformed:
        def plot(title: str, data: np.ndarray) -> None:
            create_deformed_time_slider_figure(
                triangulation,
                data,
                disp_x_history,
                disp_y_history,
                time_history,
                title,
            )
    else:
        def plot(title: str, data: np.ndarray) -> None:
            create_time_slider_figure(triangulation, data, time_history, title)

    plot("Displacement X", disp_x_history)
    plot("Displacement Y", disp_y_history)
    plot("Stress XX", stress_xx_history)
    plot("Stress XY", stress_xy_history)
    plot("Stress YY", stress_yy_history)
    plot("Uniaxial Stress", uniaxial_stress_history)
    plot("Plastic Strain", plastic_strain_history)

    if imposed_history is not None and reaction_history.size == imposed_history.size:
        plt.figure()
        plt.plot(imposed_history, reaction_history, marker="o")
        plt.xlabel("Imposed Displacement")
        plt.ylabel("Reaction Force")
        plt.title("Reaction vs imposed displacement")
        plt.grid(True)
    elif reaction_history.size:
        plt.figure()
        plt.plot(time_history, reaction_history, marker="o")
        plt.xlabel("Time")
        plt.ylabel("Reaction Force")
        plt.title("Reaction vs time")
        plt.grid(True)

    plt.show()


if __name__ == "__main__":
    main()
