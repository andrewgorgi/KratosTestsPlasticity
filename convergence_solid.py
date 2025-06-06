
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

if __name__ == "__main__":

    L2error = []
    L2error_new = []
    L_inf_error = []
    condition_number_SVD = []
    condition_number_TSVD = []

    computational_area_vec = []
    h = []
    #__________________________________________________________________________________________ 
    insertions= 30

    output_dir = "discard"

    os.makedirs(output_dir, exist_ok=True)

    tot = 50

    if os.path.exists("txt_files/Condition_numbers.txt"):
        os.remove("txt_files/Condition_numbers.txt")
    
    insertion = [14, 28,56]#, 112]#, 224]# 112]#~, 300] 
    # insertion = [ 18, 36, 72]#, 144]#, 288] # Bunny+internal
    
    # insertion = [16,32,64, 128]#, 256]#, 68, 116]#, 256]#, 256]
    # insertion = [82, 83, 84]#, 256]
    # insertion = [9, 18, 36, 72, 144]#, 288]#, 256]
    degree = 2

    insertions = 20

    # tot = 150

    tot = len(insertion)
    for i in range(0,tot) :
        # insertions = insertions+1
        insertions = insertion[i]
        # insertions = insertions*2
        print("insertions: ", insertions) 
        print("---\n\n----------------------------------------------")
        print("---\n\n----------------------------------------------")
        print("---\n\n----------------------------------------------")
        print("---\n\n----------------------------------------------")
        print("---\n\n----------------------------------------------")

        with open('ProjectParametersSolid.json','r') as f:
            parameters = KratosMultiphysics.Parameters(f.read())
        #     # Access and modify the "insert_nb_per_span_u" parameter
            parameters["modelers"][0]["Parameters"]["number_of_knot_spans"] = KratosMultiphysics.Parameters(f"[{insertions}, {insertions}]")
            parameters["modelers"][0]["Parameters"]["polynomial_order"] = KratosMultiphysics.Parameters(f"[{degree}, {degree}]")

        model = KratosMultiphysics.Model()
        simulation = StructuralMechanicsAnalysis(model,parameters)
        start_time = time.time()

        if os.path.exists("txt_files/output_results_GPs.txt"):
            os.remove("txt_files/output_results_GPs.txt")
        if os.path.exists("txt_files/Gauss_Point_coordinates.txt"):
            os.remove("txt_files/Gauss_Point_coordinates.txt")
        if os.path.exists("txt_files/boundary_GPs.txt"):
            os.remove("txt_files/boundary_GPs.txt")
        if os.path.exists("txt_files/results_on_true_boundary.txt"):
            os.remove("txt_files/results_on_true_boundary.txt")  


        # FIXME:
        # skin_model_part_inner_initial = model.CreateModelPart("skinModelPart_inner_initial")
        # skin_model_part_inner_initial.CreateNewProperties(1)
        # skin_model_part_inner_initial.CreateNewNode(1, 1.5, 0.5, 0.0)
        # skin_model_part_inner_initial.CreateNewNode(2, 3.0, 0.5, 0.0)
        # skin_model_part_inner_initial.CreateNewNode(3, 3.0, 3.5, 0.0)
        # skin_model_part_inner_initial.CreateNewNode(4, 1.5, 3.5, 0.0)
        # skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 1, [1, 2], skin_model_part_inner_initial.GetProperties()[1])
        # skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 2, [2, 3], skin_model_part_inner_initial.GetProperties()[1])
        # skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 3, [3, 4], skin_model_part_inner_initial.GetProperties()[1])
        # skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 4, [4, 1], skin_model_part_inner_initial.GetProperties()[1])

        # # Create Outer skin boundary
        # skin_model_part_outer_initial = model.CreateModelPart("skinModelPart_outer_initial")
        # skin_model_part_outer_initial.CreateNewProperties(1)
        # skin_model_part_outer_initial.CreateNewNode(1, 0.4, 0.0, 0.0)
        # skin_model_part_outer_initial.CreateNewNode(2, 3.5, 0.0, 0.0)
        # skin_model_part_outer_initial.CreateNewNode(3, 3.2, 5.0, 0.0)
        # skin_model_part_outer_initial.CreateNewNode(4, 1.0, 5.5, 0.0)
        # # skin_model_part_outer_initial.CreateNewNode(1, 0.0, 0.0, 0.0)
        # # skin_model_part_outer_initial.CreateNewNode(2, 4.0, 0.0, 0.0)
        # # skin_model_part_outer_initial.CreateNewNode(3, 4.0, 6.0, 0.0)
        # # skin_model_part_outer_initial.CreateNewNode(4, 0.0, 6.0, 0.0)
        # skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 1, [1, 2], skin_model_part_outer_initial.GetProperties()[1])
        # skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 2, [2, 3], skin_model_part_outer_initial.GetProperties()[1])
        # skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 3, [3, 4], skin_model_part_outer_initial.GetProperties()[1])
        # skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 4, [4, 1], skin_model_part_outer_initial.GetProperties()[1])
        # #FIXME:

        
        simulation.Run()

        # ERROR CALCULATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ## Compute the L^2 error
        # Exact solution as function handle:
        x = sp.symbols('x')
        y = sp.symbols('y')
        # u_exact = x

        # u_exact =  x**2*y + 0.5*x*x # u_sol

        # u_exact = x**3 + 0.2*y*x  # v_sol

        u_exact = -sp.cos(x)*sp.sinh(y)

        # u_exact = 3*x+y

        # u_exact =  x**2*y
        u_exact_handle = sp.lambdify((x, y), u_exact)

        # Computation of the error:
        output = []
        x_coord = []
        y_coord = []
        weight = []


        mp = model["IgaModelPart.StructuralAnalysisDomain"]
        L2_err_new = 0.0
        L2_norm = 0.0

        #-------------------------------------------

        computational_area = 0.0
        for elem in mp.Elements:
            geom = elem.GetGeometry()

            if (elem.GetValue(KratosMultiphysics.ACTIVATION_LEVEL) == 2): continue

            N = geom.ShapeFunctionsValues()
            weight = elem.GetValue(KratosMultiphysics.INTEGRATION_WEIGHT)
            
            # GET GAUSS POINT COORD
            stress_xx = elem.GetValue(KratosMultiphysics.IgaApplication.CAUCHY_STRESS_XX)
            stress_yy = elem.GetValue(KratosMultiphysics.IgaApplication.CAUCHY_STRESS_YY)
            stress_xy = elem.GetValue(KratosMultiphysics.IgaApplication.CAUCHY_STRESS_XY)
            # sigma_x_current_time_step.append(stress_xx);  sigma_y_current_time_step.append(stress_yy)  
            # sigma_VM_current_time_step.append(sigma_vm); sigma_xy_current_time_step.append(stress_xy)  
            
            x = 0.0
            y = 0.0
            i = 0
            curr_disp_x = 0
            curr_disp_y = 0

            for n in geom:
                curr_disp_x += N[0, i] * n.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                curr_disp_y += N[0, i] * n.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                x += N[0, i] * n.X
                y += N[0, i] * n.Y
                i += 1

            x_coord.append(x)
            y_coord.append(y)

            L2_err_new += (curr_disp_x - u_exact_handle(x,y))**2 * weight
            L2_norm += (u_exact_handle(x,y))**2 * weight

            computational_area += weight


        # sum_weight = 0
        # with open("txt_files/output_results_GPs.txt", "r") as file:
        #     for line in file:
        #         values = line.split()
        #         if len(values) >= 3:
        #             output.append(float(values[0]))
        #             x_coord.append(float(values[1]))
        #             y_coord.append(float(values[2]))
        #             weight.append(float(values[3]))
        #             sum_weight = sum_weight + float(values[3])
                    
            
            
        # print(sum_weight)


        # L2_err = 0
        # L2_err_new = 0
        # error =[] 
        # exact_solution = []
        # maximum_abs_error = 1e-14
        # count = 0
        # for i in range(len(x_coord)) :
        #     error.append(abs(output[i]-u_exact_handle(x_coord[i],y_coord[i])))
        #     exact_solution.append(u_exact_handle(x_coord[i],y_coord[i]))
        #     # if (x_coord[i]-1)**2 + (y_coord[i]-1)**2 > 0.2478**2 :
        #     L2_err = L2_err + (output[i] - u_exact_handle(x_coord[i],y_coord[i]))**2
        #     L2_err_new = L2_err_new + (output[i] - u_exact_handle(x_coord[i],y_coord[i]))**2 * weight[i]
        #     if abs(output[i] - u_exact_handle(x_coord[i],y_coord[i])) > maximum_abs_error :
        #         maximum_abs_error = abs(output[i] - u_exact_handle(x_coord[i],y_coord[i]))

        # L_inf_error.append(maximum_abs_error)

        # L2_err = np.sqrt(sum_weight * L2_err/len(x_coord))
        # L2_err = np.sqrt(L2_err/count)

        L2_err_new = np.sqrt(L2_err_new/L2_norm)

        # print('L2 error : ', L2_err)
        print('L2 error WEIGHT : ', L2_err_new)

        
        # L2error.append(L2_err)
        L2error_new.append(L2_err_new)
        computational_area_vec.append(computational_area)
        print("--> %s seconds for simulation.Run()" % (time.time() - start_time))

        # Condition number using singular values
        # matrix = scipy.io.mmread('A_1.mm')
        # largest_singular_value = scipy.sparse.linalg.svds(matrix, k=1, which='LM', return_singular_vectors=False)
        # smallest_singular_value = scipy.sparse.linalg.svds(matrix, k=1, which='SM', return_singular_vectors=False, solver='lobpcg')
        
        # # tolerance = 1e-2  # Aumenta questa tolleranza se necessario
        # # smallest_singular_value = scipy.sparse.linalg.svds(matrix, k=1, which='SM', return_singular_vectors=False, tol=tolerance)

        # # Compute the singular values using SVD
        # U, s, V = np.linalg.svd(matrix.toarray(), full_matrices=False)
        # print("Singular values:", s)
        # condition_number_SVD.append(max(s) / min(s))

        # condition_number_SVD.append(largest_singular_value / smallest_singular_value)

        # print("Condition number:", condition_number_SVD)

        h_canditate1 = 4/(insertions)

        h.append(h_canditate1)
        
        simulation.Clear()
    
    print("h = ", h )
    average_slope = 0 
    for i in range(tot-1) :
        # print(h[i], h[i+1])
        slope = (math.log(L2error_new[i+1])-math.log(L2error_new[i])) / (math.log(h[i+1])-math.log(h[i]))
        print('slope L2 = ', slope)
        average_slope = average_slope + slope
    print('average slope L2 = ', average_slope/(tot-1))
    print(h)
    # print('\n L2error', L2error)
    print('\n L2 error', L2error_new)


    plt.xscale('log')
    plt.yscale('log')
    # plt.grid(True)
    plt.grid(True, which="both", linestyle='--')

    stored = np.zeros(6)
    stored[0] = (-(1)*(np.log(h[0])) + (np.log(L2error_new[0])))
    stored[1] = (-(2)*(np.log(h[0])) + (np.log(L2error_new[0])))
    stored[2] = (-(3)*(np.log(h[0])) + (np.log(L2error_new[0])))
    stored[3] = (-(4)*(np.log(h[0])) + (np.log(L2error_new[0])))
    stored[4] = (-(5)*(np.log(h[0])) + (np.log(L2error_new[0])))
    stored[5] = (-(6)*(np.log(h[0])) + (np.log(L2error_new[0])))
    degrees = np.arange(1, 7)
    yDegrees = np.zeros((degrees.size, len(h)))
    for i in range(0, degrees.size):
        for jtest in range(0, len(h)):
            yDegrees[i][jtest] = np.exp(degrees[i] * np.log(h[jtest]) + stored[i])
        plt.plot(h, yDegrees[i], "--", label='vel %d' % tuple([degrees[i]]))

    
   


    output_file = os.path.join(output_dir, "Convergence_data.txt")
    with open(output_file, "w") as f:
        f.write("h = " + str(h) + "\n")
        f.write("L2 = " + str(L2error_new) + "\n")
        f.write("Area = " + str(computational_area_vec) + "\n")



    # plt.plot(h, L_inf_error, 's-', 'LineWidth', 2, label='L^inf error')
    # plt.plot(h, L2error, 'bs-', 'LineWidth', 2, label='L^2 error')
    plt.plot(h, L2error_new, 's-', markersize=1.5, linewidth=2, label='L^2 error with Weight')
    plt.ylabel('L2 err',fontsize=14, color='blue')
    plt.xlabel('h',fontsize=14, color='blue')
    plt.legend(loc='lower right')
    plt.savefig(os.path.join(output_dir, "L2_convergence.png"))


    fig2, ax2 = plt.subplots()

    plt.xscale('log')
    plt.yscale('log')
    # plt.grid(True)
    plt.grid(True, which="both", linestyle='--')
    ax2.plot(h, computational_area_vec, 'go-', markersize=1.5, linewidth=2, label='computational area VS refinement')

    plt.show()

