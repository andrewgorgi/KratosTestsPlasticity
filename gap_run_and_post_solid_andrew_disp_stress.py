import KratosMultiphysics
import KratosMultiphysics.IgaApplication
import KratosMultiphysics.ConstitutiveLawsApplication as CLapp
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

import sympy as sp

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.widgets import Slider

from results_utils import request_results_directory, save_simulation_results

class CustomAnalysisStage(StructuralMechanicsAnalysis):
    def __init__(self, model, parameters):
        super().__init__(model, parameters)
        self._model_ref = model
        self.time_history = []
        self.gauss_x = []
        self.gauss_y = []
        self.disp_x_history = []
        self.disp_y_history = []
        self.stress_xx_history = []
        self.stress_xy_history = []
        self.stress_yy_history = []
        
        self.strain_xx_history = []
        self.strain_xy_history = []
        self.strain_yy_history = []
        
        self.plastic_strain_history = []
        self.reaction_history = []
        self.total_area_history = []
        self.time_step = []
        self.nSTEP = 0
        self.reference_gauss_index = None
    
    def InitializeSolutionStep(self):
        current_time = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
        model["skin_model_part"].ProcessInfo[KratosMultiphysics.TIME] = current_time
        self.time_step.append(current_time)
        self.nSTEP += 1
        super().InitializeSolutionStep()

    def OutputSolutionStep(self):
        super().OutputSolutionStep()
        computing_model_part = self._GetSolver().GetComputingModelPart()
        current_time = computing_model_part.ProcessInfo[KratosMultiphysics.TIME]
        if self.time_history and abs(current_time - self.time_history[-1]) < 1e-12:
            return

        mp = computing_model_part.GetModel().GetModelPart("IgaModelPart.StructuralAnalysisDomain")

        x_coord = []
        y_coord = []
        disp_x = []
        disp_y = []
        stress_xx = []
        stress_xy = []
        stress_yy = []
        strain_xx = []
        strain_xy = []
        strain_yy = []
        
        plastic_strain_values = []
        total_area = 0.0

        for elem in mp.Elements:
            geom = elem.GetGeometry()
            N = geom.ShapeFunctionsValues()
            center = geom.Center()

            weight = elem.GetValue(KratosMultiphysics.INTEGRATION_WEIGHT)
            total_area += weight
            
            plastic_strain = elem.CalculateOnIntegrationPoints(CLapp.EQUIVALENT_PLASTIC_STRAIN, self._GetSolver().GetComputingModelPart().ProcessInfo)

            equiv_plastic_value = plastic_strain[0] if plastic_strain else 0.0
            plastic_strain_values.append(equiv_plastic_value)
        
            x_coord.append(center.X)
            y_coord.append(center.Y)

            curr_disp_x = 0.0
            curr_disp_y = 0.0
            for i, node in enumerate(geom):
                curr_disp_x += N[0, i] * node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0)
                curr_disp_y += N[0, i] * node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0)

            disp_x.append(curr_disp_x)
            disp_y.append(curr_disp_y)

            cauchy_stress = elem.GetValue(KratosMultiphysics.IgaApplication.CAUCHY_STRESS)
            stress_xx.append(cauchy_stress[0])
            stress_yy.append(cauchy_stress[1])
            stress_xy.append(cauchy_stress[2])
            
            GL_strain_tensor = elem.GetValue(KratosMultiphysics.GREEN_LAGRANGE_STRAIN_TENSOR)
            strain_xx.append(GL_strain_tensor[0, 0])
            strain_yy.append(GL_strain_tensor[1, 0])
            strain_xy.append(GL_strain_tensor[1, 1])

        if not self.gauss_x:
            self.gauss_x = x_coord
            self.gauss_y = y_coord
        if self.reference_gauss_index is None and self.gauss_x:
            target_x, target_y = -4.8, 3.8
            target_idx = 0
            min_dist = float("inf")
            for idx, (gx, gy) in enumerate(zip(self.gauss_x, self.gauss_y)):
                dx = gx - target_x
                dy = gy - target_y
                dist = dx * dx + dy * dy
                if dist < min_dist:
                    min_dist = dist
                    target_idx = idx
            self.reference_gauss_index = target_idx

        self.time_history.append(current_time)
        self.disp_x_history.append(disp_x)
        self.disp_y_history.append(disp_y)
        self.stress_xx_history.append(stress_xx)
        self.stress_xy_history.append(stress_xy)
        self.stress_yy_history.append(stress_yy)
        self.strain_xx_history.append(strain_xx)
        self.strain_xy_history.append(strain_xy)
        self.strain_yy_history.append(strain_yy)
        self.plastic_strain_history.append(plastic_strain_values)
        reaction_vals = compute_boundary_stresses_error(self._model_ref)
        reaction_value = reaction_vals[0] if reaction_vals else 0.0
        self.reaction_history.append(reaction_value)
        self.total_area_history.append(total_area)


def build_filtered_triangulation(x_coord, y_coord):
    """Create a triangulation and drop oversized triangles to improve plots."""
    triangulation = mtri.Triangulation(x_coord, y_coord)
    areas = []
    for triangle in triangulation.triangles:
        x1, y1 = x_coord[triangle[0]], y_coord[triangle[0]]
        x2, y2 = x_coord[triangle[1]], y_coord[triangle[1]]
        x3, y3 = x_coord[triangle[2]], y_coord[triangle[2]]
        area = 0.5 * abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
        areas.append(area)

    average_area = sum(areas) / len(areas)
    filtered = [tri for i, tri in enumerate(triangulation.triangles) if areas[i] <= average_area * 1.5]
    return mtri.Triangulation(x_coord, y_coord, filtered)


def _build_levels(step_min, step_max, levels_count):
    """Return strictly increasing contour levels and colorbar ticks."""
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


def create_time_slider_figure(triangulation, values_history, times, title):
    if not values_history:
        return

    fig, ax = plt.subplots(figsize=(6, 6))
    plt.subplots_adjust(bottom=0.2)

    levels_count = 100
    contour_sets = []
    colorbar_ticks = []
    for step_values in values_history:
        step_min = min(step_values)
        step_max = max(step_values)
        levels, ticks = _build_levels(step_min, step_max, levels_count)
        colorbar_ticks.append(ticks)
        contour_set = ax.tricontourf(
            triangulation,
            step_values,
            levels=levels,
            cmap="jet"
        )
        contour_sets.append(contour_set)
        for collection in contour_set.collections:
            collection.set_visible(False)

    for collection in contour_sets[0].collections:
        collection.set_visible(True)

    cbar = fig.colorbar(contour_sets[0], ax=ax, orientation="vertical")
    cbar.set_ticks(colorbar_ticks[0])
    fig._colorbar_ref = cbar

    ax.set_title(f"{title} (t = {times[0]:.6f})")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")

    slider_ax = fig.add_axes([0.15, 0.08, 0.7, 0.03])
    step_slider = Slider(slider_ax, "Step", 0, len(times) - 1, valinit=0, valstep=1, valfmt="%0.0f")

    def update(val):
        nonlocal cbar
        index = int(step_slider.val)
        for step, contour_set in enumerate(contour_sets):
            is_visible = step == index
            for collection in contour_set.collections:
                collection.set_visible(is_visible)
        cbar.remove()
        cbar = fig.colorbar(contour_sets[index], ax=ax, orientation="vertical")
        cbar.set_ticks(colorbar_ticks[index])
        fig._colorbar_ref = cbar
        ax.set_title(f"{title} (t = {times[index]:.6f})")
        ax.figure.canvas.draw_idle()

    step_slider.on_changed(update)
    fig._step_slider = step_slider
    fig._colorbar_ref = cbar

def validate_time_history(time_history):
    if not time_history:
        raise RuntimeError("No time history recorded during simulation.")
    

def plot_reaction_vs_time(time_values, reaction_values):
    if not time_values or not reaction_values:
        return
    plt.figure()
    plt.plot(time_values, reaction_values, marker="o")
    plt.xlabel("Imposed Displacement")
    plt.ylabel("Reaction force")
    plt.title("Reaction vs imposed displacement")
    plt.grid(True)


def plot_sigma_xx_vs_strain_xx(time_values, stress_history, strain_history, gauss_index=0):
    if not time_values or not stress_history or not strain_history:
        return
    if gauss_index is None:
        gauss_index = 0
    sampled_times = []
    sigma_vals = []
    strain_vals = []
    for t, sigma_step, strain_step in zip(time_values, stress_history, strain_history):
        if gauss_index < len(sigma_step) and gauss_index < len(strain_step):
            sampled_times.append(t)
            sigma_vals.append(sigma_step[gauss_index])
            strain_vals.append(strain_step[gauss_index])
    if not sampled_times:
        return
    fig, ax = plt.subplots()
    ax.plot(strain_vals, sigma_vals, color="gray", alpha=0.4, linewidth=1.0)
    scatter = ax.scatter(strain_vals, sigma_vals, c=sampled_times, cmap="viridis", marker="o")
    ax.set_xlabel("Strain xx")
    ax.set_ylabel("Sigma xx")
    ax.set_title("Sigma_xx vs Strain_xx at selected Gauss point")
    cbar = fig.colorbar(scatter, ax=ax, orientation="vertical")
    cbar.set_label("Time")
    

def plot_reaction_distribution(y_coords, reaction_contributions, title):
    if not y_coords or not reaction_contributions:
        return
    plt.figure()
    plt.scatter(y_coords, reaction_contributions, c=reaction_contributions, cmap="coolwarm", marker="o")
    plt.xlabel("y coordinate")
    plt.ylabel("Reaction contribution")
    plt.title(title)
    plt.grid(True, linestyle="--", alpha=0.4)


def compute_boundary_stresses_error(model, plot_distribution=False):
    mp = model["IgaModelPart.left"]
    
    # Lista per memorizzare i dati    
    Reaction_vals = []
    Tangential_force_vals = []
    
    reaction_val = 0.0
    
    def get_normal_stress_vec(cond):
        return np.array(cond.GetValue(KratosMultiphysics.NORMAL_STRESS), dtype=float)

    # Funzione per estrarre i dati dalle condizioni e salvarli in un array
    def accumulate(mp, plot_distribution=False):
        nonlocal reaction_val
        y_coords = []
        reaction_contributions = []
        for cond in mp.Conditions:
            # peso aggregato della condition
            w = cond.GetValue(KratosMultiphysics.INTEGRATION_WEIGHT)
        
            # direzioni bordo
            n_hat = cond.GetValue(KratosMultiphysics.NORMAL)

            # stress numerico vettoriale e proiezioni
            s_vec = get_normal_stress_vec(cond)
            sigma_nn_num   = float(np.dot(s_vec, n_hat))
            
            contribution = sigma_nn_num*w
            reaction_val += contribution
            if plot_distribution:
                geom = cond.GetGeometry()
                center = geom.Center()
                y_coords.append(center[1])
                reaction_contributions.append(sigma_nn_num)

        if plot_distribution and y_coords:
            plot_reaction_distribution(y_coords, reaction_contributions, "Reaction contribution vs y (Gap left boundary)")

    accumulate(mp, plot_distribution)
    Reaction_vals.append(reaction_val)
    
    return Reaction_vals


if __name__ == "__main__":
    results_directory = request_results_directory(__file__)
    with open("ProjectParametersSolidGap.json", "r") as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = CustomAnalysisStage(model, parameters)
    simulation.Run()

    validate_time_history(simulation.time_history)

    filtered_triangulation = build_filtered_triangulation(simulation.gauss_x, simulation.gauss_y)

    create_time_slider_figure(filtered_triangulation, simulation.disp_x_history, simulation.time_history, "Displacement X")
    # create_time_slider_figure(filtered_triangulation, simulation.disp_y_history, simulation.time_history, "Displacement Y")
    # create_time_slider_figure(filtered_triangulation, simulation.stress_xx_history, simulation.time_history, "Cauchy Stress XX")
    # create_time_slider_figure(filtered_triangulation, simulation.stress_xy_history, simulation.time_history, "Cauchy Stress XY")
    # create_time_slider_figure(filtered_triangulation, simulation.stress_yy_history, simulation.time_history, "Cauchy Stress YY")
    create_time_slider_figure(filtered_triangulation, simulation.plastic_strain_history, simulation.time_history, "Equivalent Plastic Strain")
    
    
    t = sp.symbols('t')
    u_imposed = 0.1*t
    u_imposed_handle = sp.lambdify((t), u_imposed)
    
    
    time_step_vals = simulation.time_history
    Reaction_vals = simulation.reaction_history
    # Plot final reaction distribution along the constrained boundary
    compute_boundary_stresses_error(model, plot_distribution=True)
    imposed_displacement_vals = [u_imposed_handle(ti) for ti in simulation.time_history]
    
    plot_reaction_vs_time(imposed_displacement_vals, Reaction_vals)
    plot_sigma_xx_vs_strain_xx(
        simulation.time_history,
        simulation.stress_xx_history,
        simulation.strain_xx_history,
        simulation.reference_gauss_index
    )

    save_simulation_results(
        simulation,
        results_directory,
        imposed_displacements=imposed_displacement_vals,
    )

    plt.show()
