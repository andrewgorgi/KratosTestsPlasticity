import math
import os
import time

import KratosMultiphysics
import KratosMultiphysics.IgaApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import (
    StructuralMechanicsAnalysis,
)
import numpy as np
import sympy as sp


class StructuralMechanicsAnalysisWithFlush(StructuralMechanicsAnalysis):
    def __init__(self, model, project_parameters, flush_frequency: int = 1):
        super().__init__(model, project_parameters)
        self.flush_frequency = flush_frequency
        self.last_flush = time.time()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()


if __name__ == "__main__":
    l2_error_rel = []
    computational_area_vec = []
    h = []

    output_dir = "discard"
    os.makedirs(output_dir, exist_ok=True)

    insertion = [7, 14, 28, 56, 112]
    degree = 1
    tot = len(insertion)

    # Analitico: u = cosh(x)sinh(y), v = cosh(y)sinh(x)
    x, y = sp.symbols("x y")
    u_exact = sp.cosh(x) * sp.sinh(y)
    v_exact = sp.cosh(y) * sp.sinh(x)
    u_exact_handle = sp.lambdify((x, y), u_exact)
    v_exact_handle = sp.lambdify((x, y), v_exact)

    for i in range(tot):
        insertions = insertion[i]
        print(f"insertions: {insertions}")
        print(("---\n\n----------------------------------------------") * 5)

        index_param_space = 1
        with open("ProjectParametersSolid2.json", "r", encoding="utf-8") as f:
            parameters = KratosMultiphysics.Parameters(f.read())
            parameters["modelers"][index_param_space]["Parameters"]["number_of_knot_spans"] = KratosMultiphysics.Parameters(
                f"[{insertions}, {insertions}]"
            )
            parameters["modelers"][index_param_space]["Parameters"]["polynomial_order"] = KratosMultiphysics.Parameters(
                f"[{degree}, {degree}]"
            )

        model = KratosMultiphysics.Model()
        simulation = StructuralMechanicsAnalysis(model, parameters)

        for fp in [
            "txt_files/output_results_GPs.txt",
            "txt_files/Gauss_Point_coordinates.txt",
            "txt_files/boundary_GPs.txt",
            "txt_files/results_on_true_boundary.txt",
        ]:
            if os.path.exists(fp):
                os.remove(fp)

        start_time = time.time()
        simulation.Run()
        print(f"-- > {time.time() - start_time:.3f} seconds for simulation.Run()")

        mp = model["IgaModelPart.StructuralAnalysisDomain"]

        L2_err_int = 0.0
        L2_norm_int = 0.0
        computational_area = 0.0

        for elem in mp.Elements:
            if elem.GetValue(KratosMultiphysics.ACTIVATION_LEVEL) == 2:
                continue

            geom = elem.GetGeometry()
            N = geom.ShapeFunctionsValues()
            weight = float(elem.GetValue(KratosMultiphysics.INTEGRATION_WEIGHT))

            xgp = 0.0
            ygp = 0.0
            uh = 0.0
            vh = 0.0
            for a, node in enumerate(geom):
                uh += N[0, a] * node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                vh += N[0, a] * node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                xgp += N[0, a] * node.X
                ygp += N[0, a] * node.Y

            u_ex = u_exact_handle(xgp, ygp)
            v_ex = v_exact_handle(xgp, ygp)

            L2_err_int += ((uh - u_ex) ** 2 + (vh - v_ex) ** 2) * weight
            L2_norm_int += (u_ex**2 + v_ex**2) * weight
            computational_area += weight

        L2_rel = math.sqrt(L2_err_int / L2_norm_int) if L2_norm_int > 0 else float("nan")

        print(f"L2 error (relativo, vettoriale): {L2_rel}")
        print("\n")

        l2_error_rel.append(L2_rel)
        computational_area_vec.append(computational_area)

        lower_x = parameters["modelers"][index_param_space]["Parameters"]["lower_point_xyz"][0].GetDouble()
        upper_x = parameters["modelers"][index_param_space]["Parameters"]["upper_point_xyz"][0].GetDouble()
        h_canditate = (upper_x - lower_x) / insertions
        h.append(h_canditate)

        simulation.Clear()

    print("h = ", h)
    average_slope = 0.0
    for i in range(tot - 1):
        slope = (math.log(l2_error_rel[i + 1]) - math.log(l2_error_rel[i])) / (
            math.log(h[i + 1]) - math.log(h[i])
        )
        print("slope L2 = ", slope)
        average_slope += slope
    if tot > 1:
        print("average slope L2 = ", average_slope / (tot - 1))

    print("\n L2 error (vect) =", l2_error_rel)

    # Salva su file
    output_file = os.path.join(output_dir, "Convergence_L2_data.txt")
    with open(output_file, "w", encoding="utf-8") as f:
        f.write("h = " + str(h) + "\n")
        f.write("L2_vect = " + str(l2_error_rel) + "\n")
        f.write("Area = " + str(computational_area_vec) + "\n")
