import KratosMultiphysics
# import KratosMultiphysics.ConvectionDiffusionApplication
import KratosMultiphysics.IgaApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

import time
import json
import matplotlib.pyplot as plt
import sympy as sp
import numpy as np
import math
import os
import scipy.io


class StructuralMechanicsAnalysisWithFlush(StructuralMechanicsAnalysis):

    def __init__(self,model,project_parameters,flush_frequency=1):
        super(StructuralMechanicsAnalysisWithFlush,self).__init__(model,project_parameters)
        self.flush_frequency = flush_frequency
        self.last_flush = time.time()

    def FinalizeSolutionStep(self):
        super(StructuralMechanicsAnalysisWithFlush,self).FinalizeSolutionStep()
        
    def InitializeSolutionStep(self):
        current_time = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
        # model["skin_model_part"].ProcessInfo[KratosMultiphysics.TIME] = current_time
        super().InitializeSolutionStep()

if __name__ == "__main__":

    # ---- NUOVE LISTE PER GLI ERRORI VETTORIALI ----
    L2error_rel = []   # L2 relativo su ||(u,v)||
    H1error_rel = []   # H1 relativo su u,v
    computational_area_vec = []
    h = []
    condition_numbers = []

    # __________________________________________________________________________________________ 
    output_dir = "discard"
    os.makedirs(output_dir, exist_ok=True)

    if os.path.exists("txt_files/Condition_numbers.txt"):
        os.remove("txt_files/Condition_numbers.txt")
    
    # insertion = np.arange(7, 150, 1)
    # insertion = np.arange(121, 250, 1)
    insertion = [10,20, 40, 80]#, 160]#, 320]#, 120]#, 120, 240]#, 320]#,56]#, 112]#, 224]
    degree = 1
    
    tot = len(insertion)
    for i in range(tot):
        insertions = insertion[i]
        print("insertions: ", insertions) 
        print(("---\n\n----------------------------------------------")*5)

        index_param_space = 0

        with open('ProjectParametersSolidLinear.json','r') as f:
            parameters = KratosMultiphysics.Parameters(f.read())
            parameters["modelers"][index_param_space]["Parameters"]["number_of_knot_spans"] = KratosMultiphysics.Parameters(f"[{insertions}, {insertions}]")
            parameters["modelers"][index_param_space]["Parameters"]["polynomial_order"] = KratosMultiphysics.Parameters(f"[{degree}, {degree}]")

        model = KratosMultiphysics.Model()
        simulation = StructuralMechanicsAnalysisWithFlush(model,parameters)
        start_time = time.time()

        for fp in [
            "txt_files/output_results_GPs.txt",
            "txt_files/Gauss_Point_coordinates.txt",
            "txt_files/boundary_GPs.txt",
            "txt_files/results_on_true_boundary.txt",
        ]:
            if os.path.exists(fp):
                os.remove(fp)

        simulation.Run()

        # ===================== ERROR CALCULATION (VECTOR L2 & H1) =====================
        # Analitico
        x, y = sp.symbols('x y')
        
        t = 0.0005
        # u_exact = sp.cos(x/80)*sp.sin(y/20)*t
        # v_exact = sp.sin(x/80)*sp.cos(y/20)*t
        
        # u_exact = (x*x +x*y)*t
        # v_exact = 0
        
        # u_exact = (x**2*(x-4)**2*y**2*(y-2)**2)*t
        # v_exact = 0
        
        # u_exact = (x*(x-4)*y*(y-2))*t
        # v_exact = 0
        
        ## MIDDLE
        # u_exact = -2*sp.cos(sp.pi*x/4)*sp.sin(sp.pi*y/2)*t
        # v_exact = 4*sp.sin(sp.pi*x/4)*sp.cos(sp.pi*y/2)*t
        
        #plasticity
        u_exact = -2*sp.cos(sp.pi*x/4)*sp.sin(sp.pi*(y+1)/2)*t
        v_exact = 4*sp.sin(sp.pi*x/4)*sp.cos(sp.pi* (y+1)/2)*t
        
        # u_exact = sp.cos(x)*sp.sin(y)
        # v_exact = sp.sin(x)*sp.cos(y)

        # u_exact = sp.cos(10*y)*sp.sin(10*x)
        # v_exact = sp.cos(10*x)*sp.sin(10*y)

        u_exact_handle = sp.lambdify((x, y), u_exact)
        v_exact_handle = sp.lambdify((x, y), v_exact)

        # Gradienti analitici
        ux_exact = sp.diff(u_exact, x)
        uy_exact = sp.diff(u_exact, y)
        vx_exact = sp.diff(v_exact, x)
        vy_exact = sp.diff(v_exact, y)

        ux_exact_handle = sp.lambdify((x, y), ux_exact)
        uy_exact_handle = sp.lambdify((x, y), uy_exact)
        vx_exact_handle = sp.lambdify((x, y), vx_exact)
        vy_exact_handle = sp.lambdify((x, y), vy_exact)

        mp = model["IgaModelPart.StructuralAnalysisDomain"]

        # Integrali assoluti
        L2_err_int  = 0.0   # ∫ (|u_h-u|^2 + |v_h-v|^2)
        L2_norm_int = 0.0   # ∫ (u^2 + v^2)
        G_err_int   = 0.0   # ∫ (|∇u_h-∇u|^2 + |∇v_h-∇v|^2)
        G_norm_int  = 0.0   # ∫ (|∇u|^2 + |∇v|^2)

        computational_area = 0.0

        for elem in mp.Elements:
            if elem.GetValue(KratosMultiphysics.ACTIVATION_LEVEL) == 2:
                continue

            geom = elem.GetGeometry()
            N = geom.ShapeFunctionsValues()
            weight = float(elem.GetValue(KratosMultiphysics.INTEGRATION_WEIGHT))

            # Coordinate GP (stima) + u_h(x_gp), v_h(x_gp)
            xgp = 0.0
            ygp = 0.0
            uh_x = 0.0
            vh_y = 0.0

            # valore al GP (prima riga di N)
            for a, n in enumerate(geom):
                uh_x += N[0, a] * n.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                vh_y += N[0, a] * n.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                xgp  += N[0, a] * n.X
                ygp  += N[0, a] * n.Y

            # Analitico ai GP
            u_ex  = u_exact_handle(xgp, ygp)
            v_ex  = v_exact_handle(xgp, ygp)
            dux_e = ux_exact_handle(xgp, ygp)
            duy_e = uy_exact_handle(xgp, ygp)
            dvx_e = vx_exact_handle(xgp, ygp)
            dvy_e = vy_exact_handle(xgp, ygp)

            # Contributi L2
            L2_err_int  += ((uh_x - u_ex)**2 + (vh_y - v_ex)**2) * weight
            L2_norm_int += ((u_ex)**2 + (v_ex)**2) * weight

            # Gradiente numerico al GP
            du_dx_h = 0.0; du_dy_h = 0.0
            dv_dx_h = 0.0; dv_dy_h = 0.0

            # Derivate delle shape al GP 0 (coerente col tuo uso attuale)
            dN = geom.ShapeFunctionDerivatives(1, 0)  # ordine=1, gp_index=0
            for a, n in enumerate(geom):
                ux_a = n.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                vy_a = n.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                du_dx_h += ux_a * dN[a, 0]
                du_dy_h += ux_a * dN[a, 1]
                dv_dx_h += vy_a * dN[a, 0]
                dv_dy_h += vy_a * dN[a, 1]

            # Contributi H1 (parte gradiente)
            G_err_int  += ((du_dx_h - dux_e)**2 + (du_dy_h - duy_e)**2
                           + (dv_dx_h - dvx_e)**2 + (dv_dy_h - dvy_e)**2) * weight
            G_norm_int += ((dux_e)**2 + (duy_e)**2 + (dvx_e)**2 + (dvy_e)**2) * weight

            computational_area += weight
        
        # Normalizzazioni relative
        L2_rel = np.sqrt(L2_err_int / L2_norm_int) if L2_norm_int > 0 else float('nan')
        H1_err_int  = L2_err_int + G_err_int
        H1_norm_int = L2_norm_int + G_norm_int
        H1_rel = np.sqrt(H1_err_int / H1_norm_int) if H1_norm_int > 0 else float('nan')

        print('L2 error (relativo, vettoriale): ', L2_rel)
        print('H1 error COMPLETO (relativo, vettoriale): ', H1_rel)
        print("\n")

        L2error_rel.append(L2_rel)
        H1error_rel.append(H1_rel)
        computational_area_vec.append(computational_area)
        print("--> %s seconds for simulation.Run()" % (time.time() - start_time))

        # h (adatta se preferisci altra metrica)
        lower_x = parameters["modelers"][index_param_space]["Parameters"]["lower_point_xyz"][0].GetDouble()
        upper_x = parameters["modelers"][index_param_space]["Parameters"]["upper_point_xyz"][0].GetDouble()
        h_canditate1 = (upper_x - lower_x) / (insertions)
        h.append(h_canditate1)
        
        condition_number_entries = []
        if os.path.exists("txt_files/Condition_numbers.txt"):
            with open("txt_files/Condition_numbers.txt", "r") as file:
                for line in file:
                    values = line.split()
                    if len(values) >= 2:
                        try:
                            condition_number_entries.append(float(values[1]))
                        except ValueError:
                            continue
        if condition_number_entries:
            condition_numbers.append(condition_number_entries[-1])
        else:
            condition_numbers.append(float('nan'))
        
        simulation.Clear()
    
    # -------- Pendenze L2 (opzionale) --------
    print("h = ", h )
    average_slope = 0.0
    for i in range(tot-1):
        slope = (math.log(L2error_rel[i+1]) - math.log(L2error_rel[i])) / (math.log(h[i+1]) - math.log(h[i]))
        print('slope L2 = ', slope)
        average_slope += slope
    if tot > 1:
        print('average slope L2 = ', average_slope/(tot-1))
    print(h)
    print('\n L2 error (vect) =', L2error_rel)
    print('\n H1 error COMPLETO (vect) =', H1error_rel)

    # -------- Plot --------
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True, which="both", linestyle='--')

    # reference slopes
    if len(h) >= 1:
        stored = np.zeros(6)
        for k in range(6):
            stored[k] = (-(k+1)*(np.log(h[0])) + (np.log(L2error_rel[0])))
        degrees = np.arange(1, 7)
        yDegrees = np.zeros((degrees.size, len(h)))
        for iDeg in range(degrees.size):
            for j in range(len(h)):
                yDegrees[iDeg][j] = np.exp(degrees[iDeg] * np.log(h[j]) + stored[iDeg])
            plt.plot(h, yDegrees[iDeg], "--", label='vel %d' % tuple([degrees[iDeg]]))

    output_file = os.path.join(output_dir, "Convergence_data.txt")
    with open(output_file, "w") as f:
        f.write("h = " + str(h) + "\n")
        f.write("L2_vect = " + str(L2error_rel) + "\n")
        f.write("H1_full_vect = " + str(H1error_rel) + "\n")
        f.write("Area = " + str(computational_area_vec) + "\n")

    plt.plot(h, L2error_rel, 's-', markersize=1.5, linewidth=2, label='L^2 error (rel., vect)')
    plt.plot(h, H1error_rel, 'o-', markersize=1.5, linewidth=2, label='H^1 error FULL (rel., vect)')
    plt.ylabel('Error', fontsize=14, color='blue')
    plt.xlabel('h', fontsize=14, color='blue')
    plt.legend(loc='lower right')
    plt.savefig(os.path.join(output_dir, "L2_H1_full_vector_convergence.png"))
    plt.show()
    
    valid_condition_numbers = [cn for cn in condition_numbers if not math.isnan(cn)]
    if valid_condition_numbers:
        fig2, ax2 = plt.subplots()
        ax2.set_xscale('log')
        ax2.set_yscale('log')
        ax2.grid(True, which="both", linestyle='--')
        ax2.plot(h, condition_numbers, 's-', markersize=1.5, linewidth=2, label='Condition number')
        ax2.set_ylabel('Condition number', fontsize=14, color='blue')
        ax2.set_xlabel('h', fontsize=14, color='blue')
        ax2.legend(loc='lower right')
        fig2.savefig(os.path.join(output_dir, "Condition_number.png"))

    # === Stampa comoda finale richiesta ===
    h_list = h
    L2_error = L2error_rel
    H1_error = H1error_rel
    print("\n\n\n")
    print(f"h_list = {h_list}")
    print(f"L2_error = {L2_error}")
    print(f"H1_error = {H1_error}")
    print(f"Condition_numbers = {condition_numbers}")
