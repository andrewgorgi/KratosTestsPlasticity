import KratosMultiphysics
import KratosMultiphysics.IgaApplication
import KratosMultiphysics.ConstitutiveLawsApplication as CLapp
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

import sympy as sp

import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.widgets import Slider
from time import perf_counter
import math

PLOT_DEFORMED_CONFIGURATION = True
E_ANALYTIC = 67000
NU_ANALYTIC = 0.338


class CustomAnalysisStage(StructuralMechanicsAnalysis):
    def __init__(self, model, parameters):
        super().__init__(model, parameters)
        self._model_ref = model
        self.time_history = []
        self.gauss_x = []
        self.gauss_y = []
        self.gauss_weights = []
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
        self.uniaxial_stress_history = []
        self.total_area_history = []
        self.time_step = []
        self.nSTEP = 0
        self.reference_gauss_index = None
        self.initialize_time = 0.0
        self.solve_time = 0.0

    def Initialize(self):
        start = perf_counter()
        super().Initialize()
        self.initialize_time += perf_counter() - start
    
    def InitializeSolutionStep(self):
        current_time = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
        # model["skin_model_part"].ProcessInfo[KratosMultiphysics.TIME] = current_time
        self.time_step.append(current_time)
        self.nSTEP += 1
        super().InitializeSolutionStep()

    def SolveSolutionStep(self):
        start = perf_counter()
        super().SolveSolutionStep()
        self.solve_time += perf_counter() - start

    def OutputSolutionStep(self):
        super().OutputSolutionStep()
        computing_model_part = self._GetSolver().GetComputingModelPart()
        current_time = computing_model_part.ProcessInfo[KratosMultiphysics.TIME]
        if self.time_history and abs(current_time - self.time_history[-1]) < 1e-12:
            return
        
    
        # mp = computing_model_part.GetModel().GetModelPart("IgaModelPart.SBM_Support_outer")
        # length = 0
        # for cond in mp.Conditions:
        #     weight = cond.GetValue(KratosMultiphysics.INTEGRATION_WEIGHT)
        #     length +=weight
        
        # print("total length: ", length)
        # ef

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

        uniaxial_values = []
        for elem in mp.Elements:
            geom = elem.GetGeometry()
            N = geom.ShapeFunctionsValues()
            center = geom.Center()
            center.X = 0.0
            center.Y = 0.0
            for i, node in enumerate(geom):
                center.X += N[0, i] * node.X0
                center.Y += N[0, i] * node.Y0

            weight = elem.GetValue(KratosMultiphysics.INTEGRATION_WEIGHT)
            total_area += weight
                        
            plastic_strain = elem.CalculateOnIntegrationPoints(CLapp.EQUIVALENT_PLASTIC_STRAIN, self._GetSolver().GetComputingModelPart().ProcessInfo)

            equiv_plastic_value = plastic_strain[0] if plastic_strain else 0.0
            plastic_strain_values.append(equiv_plastic_value)
        
            x_coord.append(center.X)
            y_coord.append(center.Y)
            if len(self.gauss_weights) < len(x_coord):
                self.gauss_weights.append(weight)

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
            uniaxial_stress = elem.CalculateOnIntegrationPoints(
                CLapp.UNIAXIAL_STRESS, self._GetSolver().GetComputingModelPart().ProcessInfo
            )
            uniaxial_values.append(float(uniaxial_stress[0]) if uniaxial_stress else 0.0)

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
        self.uniaxial_stress_history.append(uniaxial_values)
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


def create_deformed_time_slider_figure(triangulation, values_history, disp_x_history, disp_y_history, times, title):
    """Same as create_time_slider_figure but warps the coordinates by the nodal displacements."""
    if not values_history or not disp_x_history or not disp_y_history:
        return
    min_length = min(len(values_history), len(disp_x_history), len(disp_y_history), len(times))
    if min_length == 0:
        return
    # enforce matching lengths
    values_history = values_history[:min_length]
    disp_x_history = disp_x_history[:min_length]
    disp_y_history = disp_y_history[:min_length]
    times = times[:min_length]

    base_x = np.array(triangulation.x)
    base_y = np.array(triangulation.y)

    steps_array = [np.asarray(values_history[step]) for step in range(min_length)]
    levels_count = 100
    initial_values = steps_array[0]
    step_min = float(np.min(initial_values))
    step_max = float(np.max(initial_values))
    levels, ticks = _build_levels(step_min, step_max, levels_count)

    fig, ax = plt.subplots(figsize=(6, 6))
    plt.subplots_adjust(bottom=0.2)

    x_def = base_x + np.asarray(disp_x_history[0])
    y_def = base_y + np.asarray(disp_y_history[0])
    scatter = ax.scatter(
        x_def,
        y_def,
        c=initial_values,
        cmap="viridis",
        marker="o",
        s=15,
    )
    scatter.set_clim(levels[0], levels[-1])

    cbar = fig.colorbar(scatter, ax=ax, orientation="vertical")
    cbar.set_ticks(ticks)
    ax.set_title(f"{title} (t = {times[0]:.6f})")
    ax.set_aspect("equal", adjustable="datalim")

    axcolor = "lightgoldenrodyellow"
    slider_ax = fig.add_axes([0.15, 0.05, 0.7, 0.03], facecolor=axcolor)
    step_slider = Slider(slider_ax, "Step", 0, min_length - 1, valinit=0, valstep=1)

    def update(val):
        index = int(step_slider.val)
        x_def = base_x + np.asarray(disp_x_history[index])
        y_def = base_y + np.asarray(disp_y_history[index])
        scatter.set_offsets(np.column_stack((x_def, y_def)))
        step_values = steps_array[index]
        scatter.set_array(step_values)
        step_min = float(np.min(step_values))
        step_max = float(np.max(step_values))
        levels, ticks = _build_levels(step_min, step_max, levels_count)
        scatter.set_clim(levels[0], levels[-1])
        cbar.update_normal(scatter)
        cbar.set_ticks(ticks)
        ax.set_title(f"{title} (t = {times[index]:.6f})")
        fig._colorbar_ref = cbar
        fig.canvas.draw_idle()

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
    mp = model["IgaModelPart.SBM_Support_outer"]
    
    # Lista per memorizzare i dati
    normal_stress_vals = []
    shear_stress_vals = []
    
    Reaction_vals = []
    Tangential_force_vals = []

    # Funzione per estrarre i dati dalle condizioni e salvarli in un array
    def extract_data_from_conditions(mp, plot_distribution=False):
        ## compute normal stress L2 error
        Reaction_force = 0.0
        Tangential_force = 0.0
        y_coordinates = []
        reaction_contributions = []
        
        total_length = 0.0
        for cond in mp.Conditions:
            # Estrai i dati da RESULTS_ON_TRUE_BOUNDARY per ogni condizione
            values_on_true_boundary = cond.GetValue(KratosMultiphysics.IgaApplication.RESULTS_ON_TRUE_BOUNDARY)
            # print(f"Condition ID: {cond.Id}, Values on True Boundary: {values_on_true_boundary}")
            for i in range(np.size(values_on_true_boundary,0)):
                weight = values_on_true_boundary[i, 0]  # Il peso è nella prima colonna (indice 0)
                x = values_on_true_boundary[i, 1]  # La coordinata x è nella seconda colonna (indice 1)
                y = values_on_true_boundary[i, 2]
                normal_stress = values_on_true_boundary[i, 4]  # Stress normale (indice 4)
                shear_stress = values_on_true_boundary[i, 5]  # Shear stress (indice 5)

                # Aggiungi i dati estratti alle rispettive liste
                normal_stress_vals.append(normal_stress)
                shear_stress_vals.append(shear_stress)
                
                Reaction_force += normal_stress*weight
                Tangential_force += shear_stress*weight
                if plot_distribution:
                    y_coordinates.append(y)
                    reaction_contributions.append(normal_stress * weight)

                total_length += weight
        
        Reaction_vals.append(Reaction_force)
        Tangential_force_vals.append(Tangential_force)
        print("total reaction: ", Reaction_force)
        if plot_distribution and y_coordinates:
            plot_reaction_distribution(y_coordinates, reaction_contributions, "Reaction contribution vs y (SBM support)")
        
        
            

    # Estrai i dati per entrambe le condizioni: Top e Top_free
    extract_data_from_conditions(mp, plot_distribution)
    
    return Reaction_vals


def _plot_disp_x_numeric_vs_analytic(triangulation, numeric_disp, gauss_x, gauss_y, time_value):
    analytic_disp = np.cos(gauss_x) * np.sin(gauss_y)*time_value
    # analytic_disp = np.cos(gauss_x/80) * np.sin(gauss_y/20)*time_value
    # analytic_disp = gauss_x+gauss_y
    global_min = float(min(np.nanmin(numeric_disp), np.nanmin(analytic_disp)))
    global_max = float(max(np.nanmax(numeric_disp), np.nanmax(analytic_disp)))
    levels, ticks = _build_levels(global_min, global_max, 80)
    cmap = "coolwarm"

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    cf_num = axes[0].tricontourf(triangulation, numeric_disp, levels=levels, cmap=cmap)
    axes[0].set_title(f"DISPLACEMENT_X numerico (t = {time_value:.4f})")
    axes[0].set_aspect("equal")
    cf_an = axes[1].tricontourf(triangulation, analytic_disp, levels=levels, cmap=cmap)
    axes[1].set_title("DISPLACEMENT_X analitico")
    axes[1].set_aspect("equal")
    fig.colorbar(cf_num, ax=axes, orientation="vertical", fraction=0.04, pad=0.02, ticks=ticks)

    return analytic_disp


def _plot_disp_x_error(triangulation, numeric_disp, analytic_disp):
    error = numeric_disp - analytic_disp
    max_abs = float(np.nanmax(np.abs(error))) if error.size else 0.0
    if max_abs <= 0.0:
        max_abs = 1e-9
    levels = np.linspace(-max_abs, max_abs, 80)
    levels = np.unique(levels)
    cmap = "coolwarm"

    fig, ax = plt.subplots(figsize=(5, 5), constrained_layout=True)
    cf_err = ax.tricontourf(triangulation, error, levels=levels, cmap=cmap)
    ax.set_title("Errore DISPLACEMENT_X (num - analitico)")
    ax.set_aspect("equal")
    fig.colorbar(cf_err, ax=ax, orientation="vertical")


def _compute_analytic_stress(gauss_x, gauss_y, time_value, E=E_ANALYTIC, nu=NU_ANALYTIC):
    x = np.asarray(gauss_x)
    y = np.asarray(gauss_y)
    t = time_value
    
    sigma_xx = ((-33500000/331*t*math.sin(x)*math.sin(y)) if (math.sqrt(447561*math.sin(x)**2*math.sin(y)**2 + 328683*math.cos(x)**2*math.cos(y)**2)*abs(t) <= 664317/500000) else (-67/993*t*(1000000*abs(t) + 221439*math.sqrt(3)/math.sqrt(258748*math.sin(x)**2*math.sin(y)**2 - 109561*math.sin(x)**2 - 109561*math.sin(y)**2 + 109561))*math.sin(x)*math.sin(y)/abs(t)))
    sigma_xy = (((33500000/669)*t*math.cos(x)*math.cos(y)) if (math.sqrt(447561*math.sin(x)**2*math.sin(y)**2 + 328683*math.cos(x)**2*math.cos(y)**2)*abs(t) <= 664317/500000) else (22177*math.sqrt(3)*t*math.cos(x)*math.cos(y)/(math.sqrt(149187*math.sin(x)**2*math.sin(y)**2 + 109561*math.cos(x)**2*math.cos(y)**2)*abs(t))))
    sigma_yy = ((-33500000/331*t*math.sin(x)*math.sin(y)) if (math.sqrt(447561*math.sin(x)**2*math.sin(y)**2 + 328683*math.cos(x)**2*math.cos(y)**2)*abs(t) <= 664317/500000) else (-67/993*t*(1000000*abs(t) + 221439*math.sqrt(3)/math.sqrt(258748*math.sin(x)**2*math.sin(y)**2 - 109561*math.sin(x)**2 - 109561*math.sin(y)**2 + 109561))*math.sin(x)*math.sin(y)/abs(t)))

    return sigma_xx, sigma_xy, sigma_yy


def _plot_stress_numeric_vs_analytic(triangulation, numeric, analytic, title, time_value):
    global_min = float(np.nanmin([np.nanmin(numeric), np.nanmin(analytic)]))
    global_max = float(np.nanmax([np.nanmax(numeric), np.nanmax(analytic)]))
    levels = np.linspace(global_min, global_max, 80)
    levels = np.unique(levels)
    cmap = "coolwarm"

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    cf_num = axes[0].tricontourf(triangulation, numeric, levels=levels, cmap=cmap)
    axes[0].set_title(f"{title} numerico (t = {time_value:.4f})")
    axes[0].set_aspect("equal")
    cf_an = axes[1].tricontourf(triangulation, analytic, levels=levels, cmap=cmap)
    axes[1].set_title(f"{title} analitico")
    axes[1].set_aspect("equal")
    fig.colorbar(cf_num, ax=axes, orientation="vertical", fraction=0.04, pad=0.02)


def _compute_von_mises(sxx: np.ndarray, syy: np.ndarray, sxy: np.ndarray) -> np.ndarray:
    sxx = np.asarray(sxx, dtype=float)
    syy = np.asarray(syy, dtype=float)
    sxy = np.asarray(sxy, dtype=float)
    return np.sqrt(np.maximum(0.0, sxx**2 - sxx * syy + syy**2 + 3.0 * sxy**2))


def _plot_von_mises_vs_uniaxial(triangulation, vm: np.ndarray, uniaxial: np.ndarray, time_value: float):
    vm = np.asarray(vm, dtype=float)
    uniaxial = np.asarray(uniaxial, dtype=float)
    global_min = float(np.nanmin([np.nanmin(vm), np.nanmin(uniaxial)]))
    global_max = float(np.nanmax([np.nanmax(vm), np.nanmax(uniaxial)]))
    levels = np.linspace(global_min, global_max, 80)
    levels = np.unique(levels)
    cmap = "coolwarm"

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    cf_vm = axes[0].tricontourf(triangulation, vm, levels=levels, cmap=cmap)
    axes[0].set_title(f"Von Mises numerico (t = {time_value:.4f})")
    axes[0].set_aspect("equal")
    cf_uni = axes[1].tricontourf(triangulation, uniaxial, levels=levels, cmap=cmap)
    axes[1].set_title("Uniaxial stress (numerico)")
    axes[1].set_aspect("equal")
    fig.colorbar(cf_vm, ax=axes, orientation="vertical", fraction=0.04, pad=0.02)


def _plot_vm_difference(triangulation, vm_analytic: np.ndarray, uniaxial: np.ndarray, time_value: float):
    diff = vm_analytic - np.asarray(uniaxial, dtype=float)
    max_abs = float(np.nanmax(np.abs(diff))) if diff.size else 0.0
    if max_abs <= 0.0:
        max_abs = 1e-9
    levels = np.linspace(-max_abs, max_abs, 80)
    levels = np.unique(levels)
    cmap = "coolwarm"
    fig, ax = plt.subplots(figsize=(5, 5), constrained_layout=True)
    cf = ax.tricontourf(triangulation, diff, levels=levels, cmap=cmap)
    ax.set_title(f"Diff Von Mises analitico - uniaxial (t = {time_value:.4f})")
    ax.set_aspect("equal")
    fig.colorbar(cf, ax=ax, orientation="vertical")


def _l2_relative(numeric: np.ndarray, analytic: np.ndarray, weights: np.ndarray) -> float:
    if numeric.size == 0 or analytic.size == 0 or weights.size == 0:
        return float("nan")
    numeric = np.asarray(numeric, dtype=float)
    analytic = np.asarray(analytic, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if numeric.shape != analytic.shape or numeric.shape != weights.shape:
        return float("nan")
    num = np.sum(weights * (numeric - analytic) ** 2)
    den = np.sum(weights * analytic**2)
    if den <= 0:
        return float("nan")
    return float(np.sqrt(num / den))
    
    
if __name__ == "__main__":
    if os.path.exists("txt_files/Id_active_control_points_condition.txt"):
        os.remove("txt_files/Id_active_control_points_condition.txt")
    if os.path.exists("txt_files/Id_active_control_points.txt"):
        os.remove("txt_files/Id_active_control_points.txt")
    if os.path.exists("txt_files/boundary_GPs.txt"):
        os.remove("txt_files/boundary_GPs.txt")
    if os.path.exists("txt_files/results_on_true_boundary.txt"):
        os.remove("txt_files/results_on_true_boundary.txt")  
        
        
    with open("ProjectParametersSolid2.json", "r") as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = CustomAnalysisStage(model, parameters)
    simulation.Run()
    print(f"Tempo Initialize: {simulation.initialize_time:.3f} s")
    print(f"Tempo Solve (solo solver): {simulation.solve_time:.3f} s")

    validate_time_history(simulation.time_history)

    filtered_triangulation = build_filtered_triangulation(simulation.gauss_x, simulation.gauss_y)

    if PLOT_DEFORMED_CONFIGURATION:
        def plot_on_deformed(values_history, title):
            create_deformed_time_slider_figure(
                filtered_triangulation,
                values_history,
                simulation.disp_x_history,
                simulation.disp_y_history,
                simulation.time_history,
                title
            )

        # plot_on_deformed(simulation.disp_x_history, "Displacement X on deformed geometry")
        # plot_on_deformed(simulation.disp_y_history, "Displacement Y on deformed geometry")
        # plot_on_deformed(simulation.stress_xx_history, "Cauchy Stress XX on deformed geometry")
        # plot_on_deformed(simulation.stress_xy_history, "Cauchy Stress XY on deformed geometry")
        # plot_on_deformed(simulation.stress_yy_history, "Cauchy Stress YY on deformed geometry")
        # plot_on_deformed(simulation.plastic_strain_history, "Equivalent Plastic Strain on deformed geometry")
        # plot_on_deformed(simulation.uniaxial_stress_history, "Uniaxial stress on deformed geometry")
    else:
        create_time_slider_figure(filtered_triangulation, simulation.disp_x_history, simulation.time_history, "Displacement X")
        create_time_slider_figure(filtered_triangulation, simulation.disp_y_history, simulation.time_history, "Displacement Y")
        create_time_slider_figure(filtered_triangulation, simulation.stress_xx_history, simulation.time_history, "Cauchy Stress XX")
        create_time_slider_figure(filtered_triangulation, simulation.stress_xy_history, simulation.time_history, "Cauchy Stress XY")
        create_time_slider_figure(filtered_triangulation, simulation.stress_yy_history, simulation.time_history, "Cauchy Stress YY")
        create_time_slider_figure(filtered_triangulation, simulation.plastic_strain_history, simulation.time_history, "Equivalent Plastic Strain")
    
    
    t = sp.symbols('t')
    u_imposed = 250*t
    u_imposed_handle = sp.lambdify((t), u_imposed)
    
    
    time_step_vals = simulation.time_history
    Reaction_vals = simulation.reaction_history
    imposed_displacement_vals = [u_imposed_handle(ti) for ti in simulation.time_history]
    print(imposed_displacement_vals)
    
    # plot_reaction_vs_time(imposed_displacement_vals, Reaction_vals)
    # plot_sigma_xx_vs_strain_xx(
    #     simulation.time_history,
    #     simulation.stress_xx_history,
    #     simulation.strain_xx_history,
    #     simulation.reference_gauss_index
    # )
    # # Plot final reaction distribution along the boundary
    # compute_boundary_stresses_error(model, plot_distribution=True)

    # Confronto con soluzione analitica u = cosh(x) sinh(y), v = cosh(y) sinh(x)
    last_disp_x = np.asarray(simulation.disp_x_history[-1]) if simulation.disp_x_history else np.array([])
    gauss_x = np.asarray(simulation.gauss_x)
    gauss_y = np.asarray(simulation.gauss_y)
    if last_disp_x.size and gauss_x.size and gauss_y.size:
        analytic_disp = _plot_disp_x_numeric_vs_analytic(
            filtered_triangulation,
            last_disp_x,
            gauss_x,
            gauss_y,
            simulation.time_history[-1] if simulation.time_history else 0.0,
        )
        _plot_disp_x_error(filtered_triangulation, last_disp_x, analytic_disp)

        # Stress analytic vs numerico all'ultimo step
        last_stress_xx = np.asarray(simulation.stress_xx_history[-1]) if simulation.stress_xx_history else np.array([])
        last_stress_xy = np.asarray(simulation.stress_xy_history[-1]) if simulation.stress_xy_history else np.array([])
        last_stress_yy = np.asarray(simulation.stress_yy_history[-1]) if simulation.stress_yy_history else np.array([])
        last_uniaxial = np.asarray(simulation.uniaxial_stress_history[-1]) if simulation.uniaxial_stress_history else np.array([])
        if last_stress_xx.size and last_stress_xy.size and last_stress_yy.size:
            sigma_xx_an, sigma_xy_an, sigma_yy_an = _compute_analytic_stress(
                gauss_x, gauss_y, simulation.time_history[-1]
            )
            _plot_stress_numeric_vs_analytic(
                filtered_triangulation, last_stress_xx, sigma_xx_an, "Cauchy Stress XX", simulation.time_history[-1]
            )
            _plot_stress_numeric_vs_analytic(
                filtered_triangulation, last_stress_xy, sigma_xy_an, "Cauchy Stress XY", simulation.time_history[-1]
            )
            _plot_stress_numeric_vs_analytic(
                filtered_triangulation, last_stress_yy, sigma_yy_an, "Cauchy Stress YY", simulation.time_history[-1]
            )
            weights = np.asarray(simulation.gauss_weights, dtype=float)
            if weights.size == gauss_x.size:
                ana_disp = np.cos(gauss_x) * np.sin(gauss_y) * simulation.time_history[-1]
                err_disp = _l2_relative(last_disp_x, ana_disp, weights)
                err_sxx = _l2_relative(last_stress_xx, sigma_xx_an, weights)
                err_sxy = _l2_relative(last_stress_xy, sigma_xy_an, weights)
                err_syy = _l2_relative(last_stress_yy, sigma_yy_an, weights)
                print(f"L2 rel disp_x: {err_disp:.6e}")
                print(f"L2 rel sigma_xx: {err_sxx:.6e}")
                print(f"L2 rel sigma_xy: {err_sxy:.6e}")
                print(f"L2 rel sigma_yy: {err_syy:.6e}")
            if last_uniaxial.size:
                vm_analytic = _compute_von_mises(sigma_xx_an, sigma_yy_an, sigma_xy_an)
                print(vm_analytic)
                _plot_von_mises_vs_uniaxial(
                    filtered_triangulation, vm_analytic, last_uniaxial, simulation.time_history[-1]
                )
                _plot_vm_difference(
                    filtered_triangulation, vm_analytic, last_uniaxial, simulation.time_history[-1]
                )

    plt.show()
