import KratosMultiphysics
import KratosMultiphysics.IgaApplication
from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis

import time
import matplotlib.pyplot as plt
import sympy as sp
import numpy as np
import math
import os


class CustomAnalysisStage(ConvectionDiffusionAnalysis):

    def __init__(self,model,project_parameters,flush_frequency=1):
        super(ConvectionDiffusionAnalysis,self).__init__(model,project_parameters)
        self.flush_frequency = flush_frequency
        self.last_flush = time.time()
    
# Define the analytical solution
def analytical_temperature(x, y, t):
    return np.sin(x)*np.sinh(y)
    # return x**2+y**2
    # return x**3+y**3
    # return x+y


if __name__ == "__main__":

    L2error = []
    L_inf_error = []

    computational_area_vec = []
    h = []
    #__________________________________________________________________________________________ 
    insertions= 15

    # output_dir = "ResultsPolynomialSol2/p2/non_linear_problem"
    output_dir = "discard"
    # Crea la cartella se non esiste
    os.makedirs(output_dir, exist_ok=True)

    # tot = 40

    if os.path.exists("txt_files/Condition_numbers.txt"):
        os.remove("txt_files/Condition_numbers.txt")
    
    insertion = [14, 28, 56, 112]#, 224] 
    # insertion = [ 18, 36, 72, 144]#, 288] # Bunny+internal
    
    # insertion = [15, 21, 23, 24, 26, 31, 62, 63, 64]#, 256]#, 256]
    # insertion = [82, 83, 84]#, 256]
    # insertion = [9, 18, 36, 72, 144]#, 256]
    degree = 2

    tot = 50

    tot = len(insertion)
    for i in range(0,tot) :
        # insertions = insertions+1
        insertions = insertion[i]
        # insertions = insertions*2
        print("insertions: ", insertions) 
        print("---\n\n----------------------------------------------")

        with open('ProjectParameters.json','r') as f:
            parameters = KratosMultiphysics.Parameters(f.read())
        #     # Access and modify the "insert_nb_per_span_u" parameter
            parameters["modelers"][2]["Parameters"]["number_of_knot_spans"] = KratosMultiphysics.Parameters(f"[{insertions}, {insertions}]")
            parameters["modelers"][2]["Parameters"]["polynomial_order"] = KratosMultiphysics.Parameters(f"[{degree}, {degree}]")

        model = KratosMultiphysics.Model()
        simulation = CustomAnalysisStage(model, parameters)
        simulation.Run()

        # Extract the computed solution at a specific time step
        mp = model["IgaModelPart"]

        x_coord = []
        y_coord = []
        computed_temperature = []
        weights = []

        # Loop over elements to gather computed solution
        for elem in mp.Elements:
            if (elem.Id == 1):
                continue
            geom = elem.GetGeometry()
            N = geom.ShapeFunctionsValues()
            center = geom.Center()
            weight = elem.GetValue(KratosMultiphysics.INTEGRATION_WEIGHT)
            weights.append(weight)
            # Extract Gauss point (center) coordinates
            x_coord.append(center.X)
            y_coord.append(center.Y)

            # Initialize solution values at the center
            curr_temperature = 0
            # Compute nodal contributions using shape functions
            for i, node in enumerate(geom):
                curr_temperature += N[0, i] * node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0)
            computed_temperature.append(curr_temperature)


        # Get the current time after simulation run
        current_time = simulation._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
        print("Current time after simulation:", current_time)
        
        # Calculate errors and analytical solution
        temperature_error = []
        temperature_analytical_values = []

        for i in range(len(x_coord)):
            x = x_coord[i]
            y = y_coord[i]

            # Analytical solution at the point
            temperature_analytical = analytical_temperature(x, y, current_time)

            # Append analytical values for plotting
            temperature_analytical_values.append(temperature_analytical)

            # Error calculations
            temperature_error.append(abs(computed_temperature[i] - temperature_analytical))
        
        # Compute the L2 norm of the error for temperature
        temperature_l2_norm = 0
        for i in range(len(weights)):
            temperature_l2_norm += weights[i] * temperature_error[i]**2
        # Take the square root to finalize the L2 norm
        temperature_l2_norm = (temperature_l2_norm)**0.5
        
        # Compute the L2 norm of the analytical solution
        analytical_l2_norm = 0
        for i in range(len(weights)):
            analytical_l2_norm += weights[i] * temperature_analytical_values[i]**2
        analytical_l2_norm = (analytical_l2_norm)**0.5

        # Compute the normalized L2 error
        normalized_l2_error = temperature_l2_norm / analytical_l2_norm

        # Print the results
        print("L2 norm of temperature error:", temperature_l2_norm)
        # print("L2 norm of analytical solution:", analytical_l2_norm)
        print("Normalized L2 error:", normalized_l2_error)

        
        L2error.append(normalized_l2_error)


        h_canditate1 = 2/(insertions)

        h.append(h_canditate1)
        
        simulation.Clear()
    
    print("h = ", h )
    average_slope = 0 
    for i in range(tot-1) :
        # print(h[i], h[i+1])
        slope = (math.log(L2error[i+1])-math.log(L2error[i])) / (math.log(h[i+1])-math.log(h[i]))
        print('slope L2 = ', slope)
        average_slope = average_slope + slope
    print('average slope L2 = ', average_slope/(tot-1))
    print(h)
    print('\n L2 error', L2error)


    plt.xscale('log')
    plt.yscale('log')
    # plt.grid(True)
    plt.grid(True, which="both", linestyle='--')

    stored = np.zeros(5)
    stored[0] = (-(1+1)*(np.log(h[0])) + (np.log(L2error[0])))
    stored[1] = (-(2+1)*(np.log(h[0])) + (np.log(L2error[0])))
    stored[2] = (-(3+1)*(np.log(h[0])) + (np.log(L2error[0])))
    stored[3] = (-(4+1)*(np.log(h[0])) + (np.log(L2error[0])))
    stored[4] = (-(5+1)*(np.log(h[0])) + (np.log(L2error[0])))
    degrees = np.arange(2, 7)
    yDegrees = np.zeros((degrees.size, len(h)))
    for i in range(0, degrees.size):
        for jtest in range(0, len(h)):
            yDegrees[i][jtest] = np.exp(degrees[i] * np.log(h[jtest]) + stored[i])
        plt.plot(h, yDegrees[i], "--", label='vel %d' % tuple([degrees[i]]))

    
   


    output_file = os.path.join(output_dir, "Convergence_data.txt")
    with open(output_file, "w") as f:
        f.write("h = " + str(h) + "\n")
        f.write("L2 = " + str(L2error) + "\n")
        f.write("Area = " + str(computational_area_vec) + "\n")



    # plt.plot(h, L_inf_error, 's-', 'LineWidth', 2, label='L^inf error')
    plt.plot(h, L2error, 's-', markersize=1.5, linewidth=2, label='L^2 error with Weight')
    plt.ylabel('L2 err',fontsize=14, color='blue')
    plt.xlabel('h',fontsize=14, color='blue')
    plt.legend(loc='lower right')
    plt.savefig(os.path.join(output_dir, "L2_convergence.png"))


