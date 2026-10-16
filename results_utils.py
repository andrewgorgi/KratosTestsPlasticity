import json
from datetime import datetime
from pathlib import Path
from typing import Iterable, Optional

import numpy as np


def _sanitize_folder_name(raw_name: str) -> str:
    """Return a filesystem-safe version of the requested folder name."""
    sanitized = "".join(
        ch if ch.isalnum() or ch in ("-", "_") else "_" for ch in raw_name
    ).strip("_")
    return sanitized


def request_results_directory(script_file: str, base_folder: str = "Results") -> Path:
    """
    Ask the user for the name of the results folder and ensure it exists.

    The folder is created inside the provided base_folder relative to the script.
    """
    script_dir = Path(script_file).resolve().parent
    base_dir = script_dir / base_folder
    base_dir.mkdir(parents=True, exist_ok=True)
    sanitized = ""
    while not sanitized:
        raw_input = input(f"Inserisci il nome della sottocartella da creare in {base_folder}/: ").strip()
        sanitized = _sanitize_folder_name(raw_input)
        if not sanitized:
            print("Il nome indicato non è valido. Inserisci almeno un carattere alfanumerico, '-' o '_'.")

    candidate = base_dir / sanitized
    suffix = 1
    while candidate.exists():
        candidate = base_dir / f"{sanitized}_{suffix:02d}"
        suffix += 1

    candidate.mkdir(parents=True, exist_ok=False)
    print(f"I risultati verranno salvati in: {candidate}")
    return candidate


def _as_float_array(sequence: Iterable) -> np.ndarray:
    seq_list = list(sequence)
    if not seq_list:
        return np.array([], dtype=float)
    return np.array(seq_list, dtype=float)


def save_simulation_results(
    simulation,
    output_dir: Path,
    imposed_displacements: Optional[Iterable[float]] = None,
) -> None:
    """
    Persist the data gathered during the simulation run inside output_dir.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    history_data = {
        "time_history": _as_float_array(simulation.time_history),
        "gauss_x": _as_float_array(simulation.gauss_x),
        "gauss_y": _as_float_array(simulation.gauss_y),
        "disp_x_history": _as_float_array(simulation.disp_x_history),
        "disp_y_history": _as_float_array(simulation.disp_y_history),
        "stress_xx_history": _as_float_array(simulation.stress_xx_history),
        "stress_xy_history": _as_float_array(simulation.stress_xy_history),
        "stress_yy_history": _as_float_array(simulation.stress_yy_history),
        "uniaxial_stress_history": _as_float_array(simulation.uniaxial_stress_history),
        "strain_xx_history": _as_float_array(simulation.strain_xx_history),
        "strain_xy_history": _as_float_array(simulation.strain_xy_history),
        "strain_yy_history": _as_float_array(simulation.strain_yy_history),
        "plastic_strain_history": _as_float_array(simulation.plastic_strain_history),
        "reaction_history": _as_float_array(simulation.reaction_history),
        "total_area_history": _as_float_array(simulation.total_area_history),
        "time_step_samples": _as_float_array(simulation.time_step),
    }

    if imposed_displacements is not None:
        history_data["imposed_displacement_history"] = _as_float_array(imposed_displacements)

    np.savez(output_dir / "field_histories.npz", **history_data)

    metadata = {
        "steps_recorded": int(len(simulation.time_history)),
        "gauss_points": int(len(simulation.gauss_x)),
        "reference_gauss_index": simulation.reference_gauss_index,
        "results_folder": str(output_dir),
    }
    if hasattr(simulation, "initialize_time"):
        metadata["initialize_time_sec"] = float(simulation.initialize_time)
    if hasattr(simulation, "solve_time"):
        metadata["solve_time_sec"] = float(simulation.solve_time)
    with open(output_dir / "summary.json", "w", encoding="utf-8") as metadata_file:
        json.dump(metadata, metadata_file, indent=2)

    time_len = len(simulation.time_history)
    reaction_len = len(simulation.reaction_history)

    if time_len and time_len == reaction_len:
        columns = [np.array(simulation.time_history, dtype=float)]
        headers = ["time"]

        if imposed_displacements is not None:
            imposed_array = np.array(list(imposed_displacements), dtype=float)
            if len(imposed_array) == time_len:
                columns.append(imposed_array)
                headers.append("imposed_displacement")

        columns.append(np.array(simulation.reaction_history, dtype=float))
        headers.append("reaction_force")

        reaction_table = np.column_stack(columns)
        header_line = ",".join(headers)
        np.savetxt(
            output_dir / "reaction_history.csv",
            reaction_table,
            delimiter=",",
            header=header_line,
            comments="",
        )
