import KratosMultiphysics
import KratosMultiphysics.IgaApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
import os

import matplotlib.tri as mtri

from matplotlib.ticker import FormatStrFormatter
import sympy as sp
import numpy as np

import matplotlib.pyplot as plt
import json
import matplotlib.patches as patches

from matplotlib.widgets import Slider

# import imageio

###########################################################################33
##############################################################################
def PlotVariable(ax, fig, n_objects, n_time_steps, variable_to_print, deformed_triangulations):
    lines = [[[] for _ in range(len(variable_to_print)) ] for _ in range(n_objects)]
    for i_body in range(n_objects):
        for i_step in range(n_time_steps):
            new_triangulation = deformed_triangulations[i_body][i_step]
            

            for i_var in range(0, len(variable_to_print)):
                # Crea i livelli di contorno
                num_levels = 100
                minmin = min(np.min(variable_to_print[i_var][0][i_step]), np.min(variable_to_print[i_var][1][i_step]))
                maxmax = max(np.max(variable_to_print[i_var][0][i_step]), np.max(variable_to_print[i_var][1][i_step]))

                # minmin = np.min(variable_to_print[i_var][1][i_step])
                # maxmax = np.max(variable_to_print[i_var][1][i_step])

                # for jj in range(0, len(variable_to_print[i_var][0][i_step])):
                #     if variable_to_print[i_var][0][i_step][jj] == maxmax:
                #        variable_to_print[i_var][0][i_step][jj] = 0.0

                # minmin = min(np.min(variable_to_print[i_var][0][i_step]), np.min(variable_to_print[i_var][1][i_step]))
                # maxmax = max(np.max(variable_to_print[i_var][0][i_step]), np.max(variable_to_print[i_var][1][i_step]))


                # print("minmin: %f maxmax: %f \n" % (minmin, maxmax))

                if (minmin > maxmax - 1e-6): maxmax = minmin + 1
                levels = np.linspace(minmin, maxmax, num_levels)
                
                i_plot_x = i_var%2
                i_plot_y = int(i_var/2)
                contour_set = ax[i_var].tricontourf(new_triangulation, variable_to_print[i_var][i_body][i_step], levels, alpha = 1.0, cmap='jet')
                # Nascondi tutte le collezioni inizialmente
                for collection in contour_set.collections:
                    collection.set_visible(False)
                lines[i_body][i_var].append(contour_set)



    # Creare la barra dei colori per il primo step
    cbar_list = []
    for i_var in range(0, len(variable_to_print)):
        i_plot_x = i_var%2
        i_plot_y = int(i_var/2)
        
        minVar = min(np.min(variable_to_print[i_var][0][0]), np.min(variable_to_print[i_var][1][0]))
        maxVar = max(np.max(variable_to_print[i_var][0][0]), np.max(variable_to_print[i_var][1][0]))
        cbar = fig.colorbar(lines[1][i_var][0], ax=ax[i_var])
        cbar.set_ticks(np.linspace(minVar, maxVar, 10))
        cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.5f'))

        cbar_list.append(cbar)
    
    return cbar_list, lines

# Funzione per aggiornare il grafico con lo slider
def update_plot(index, variable_to_print, cbar_list, lines, fig, ax):
    # Nascondi tutte le collezioni inizialmente
    for i_body in range(n_objects):
        for i_var in range(len(variable_to_print)):
            for i, line in enumerate(lines[i_body][i_var]):
                is_visible = (i == index)
                for collection in line.collections:
                    collection.set_visible(is_visible)

    # Per ogni variabile e plot, aggiorna i valori di contorno e la barra dei colori
    for i_var in range(len(variable_to_print)):
        i_plot_x = i_var % 2
        i_plot_y = int(i_var / 2)

        # Ottieni i valori correnti
        current_values = variable_to_print[i_var][0][index], variable_to_print[i_var][1][index]
        minVar = min(np.min(current_values[0]), np.min(current_values[1]))
        maxVar = max(np.max(current_values[0]), np.max(current_values[1]))

        # Aggiorna la barra dei colori
        cbar_list[i_var].remove()  # Rimuovi la barra dei colori esistente
        cbar_list[i_var] = fig.colorbar(lines[i_body][i_var][index], ax=ax[i_var])
        cbar_list[i_var].set_ticks(np.linspace(minVar, maxVar, 10))
        cbar_list[i_var].ax.yaxis.set_major_formatter(FormatStrFormatter('%.5f'))

    # Aggiorna il titolo
    fig.suptitle(f"Grafico TIME: {simulation.time_step[-n_time_steps + index]}")
    fig.canvas.draw_idle()

    # # Imposta i limiti del primo grafico (ax[0])
    # ax[0].set_xlim(-1, 17)    # Limiti per l'asse x
    # ax[0].set_ylim(-1, 17)    # Limiti per l'asse y

    # # Imposta i limiti del secondo grafico (ax[1])
    # ax[1].set_xlim(-1, 17)    # Limiti per l'asse x
    # ax[1].set_ylim(-1, 17)  # Limiti per l'asse y




class CustomAnalysisStage(StructuralMechanicsAnalysis):
    
    time_step = []
    nSTEP = 0
    reaction_force_body1 = [[], []]
    reaction_force_body2 = [[], []]
    
    contact_force_body1 = [[],[]]
    contact_force_body2 = [[],[]]

    sigma_x = [[],[]]
    sigma_y = [[],[]]
    sigma_xy = [[],[]]
    sigma_VM =[[],[]]

    old_time = 0

    def Initialize(self):
        super().Initialize()
        self.old_time = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]


    def InitializeSolutionStep(self):   
        # Qui puoi accedere al current time step
        current_time = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]

        model["skin_Body1"].ProcessInfo[KratosMultiphysics.TIME] = current_time
        model["skin_Body2"].ProcessInfo[KratosMultiphysics.TIME] = current_time

        super().InitializeSolutionStep()
        # self.time_step.append(current_time)
    
    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """
        while self.KeepAdvancingSolutionLoop():
            self.time = self._AdvanceTime()
            
            is_converged = False

            it = 0
            # while (it < 2):
            # while (is_converged == False) :
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            is_converged = self._GetSolver().SolveSolutionStep()
                # self._AnalysisStage__CheckIfSolveSolutionStepReturnsAValue(is_converged)
            # it += 1
            # is_converged = False

            self.FinalizeSolutionStep()


            
            self.OutputSolutionStep()
        
        
    def OutputSolutionStep(self):
        super().OutputSolutionStep()
        # Qui puoi accedere al numero di step attuale
        # mp = model["IgaModelPart.Body1.Support_outer"]
        current_time = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
        if (current_time == self.old_time): return
        self.old_time = current_time
        self.time_step.append(current_time)
        self.nSTEP += 1
        
        mp = model["IgaModelPart.Body1.StructuralAnalysisDomain"]
        sigma_x_current_time_step = []
        sigma_y_current_time_step = []
        sigma_VM_current_time_step = []
        sigma_xy_current_time_step = []

        s_max = 0
        ppp = 0
        for elem in mp.Elements:
            geom = elem.GetGeometry()
            N = geom.ShapeFunctionsValues()
            weight = elem.GetValue(KratosMultiphysics.INTEGRATION_WEIGHT)
            
            # GET GAUSS POINT COORD
            stress_xx = elem.GetValue(KratosMultiphysics.IgaApplication.CAUCHY_STRESS_XX)
            stress_yy = elem.GetValue(KratosMultiphysics.IgaApplication.CAUCHY_STRESS_YY)
            stress_xy = elem.GetValue(KratosMultiphysics.IgaApplication.CAUCHY_STRESS_XY)
            sigma_vm = elem.GetValue(KratosMultiphysics.YIELD_STRESS)
            sigma_x_current_time_step.append(stress_xx);  sigma_y_current_time_step.append(stress_yy)  
            sigma_VM_current_time_step.append(sigma_vm); sigma_xy_current_time_step.append(stress_xy)  
            if (sigma_vm > s_max): 
                gauss_pos = np.zeros(2)
                i = 0
                for n in geom:
                    gauss_pos[0] += N[0, i] * n.X
                    gauss_pos[1] += N[0, i] * n.Y
                    i += 1

                s_max = sigma_vm
                ppp = gauss_pos
        
        t_vec = []
        VM_vec = []
        s_xx_vec = []
        # for cond in mp.Conditions:
        #     geom = cond.GetGeometry()
        #     N = geom.ShapeFunctionsValues()
        #     weight = cond.GetValue(KratosMultiphysics.INTEGRATION_WEIGHT)
            
        #     # GET GAUSS POINT COORD
        #     stress_xx = cond.GetValue(KratosMultiphysics.IgaApplication.CAUCHY_STRESS_XX)
        #     stress_yy = cond.GetValue(KratosMultiphysics.IgaApplication.CAUCHY_STRESS_YY)
        #     stress_xy = cond.GetValue(KratosMultiphysics.IgaApplication.CAUCHY_STRESS_XY)
        #     sigma_vm = cond.GetValue(KratosMultiphysics.YIELD_STRESS)
        #     sigma_x_current_time_step.append(stress_xx);  sigma_y_current_time_step.append(stress_yy)   
        #     sigma_VM_current_time_step.append(sigma_vm); sigma_xy_current_time_step.append(stress_xy)   

        #     sigma_n = cond.GetValue(KratosMultiphysics.NORMAL_STRESS)

        #     gauss_pos = np.zeros(2)
        #     i = 0
        #     for n in geom:
        #         gauss_pos[0] += N[0, i] * n.X
        #         gauss_pos[1] += N[0, i] * n.Y
        #         i += 1

        #     if (sigma_vm > s_max): 
        #         s_max = sigma_vm
        #         ppp = gauss_pos

        self.sigma_x[0].append(sigma_x_current_time_step);  self.sigma_y[0].append(sigma_y_current_time_step)   
        self.sigma_VM[0].append(sigma_VM_current_time_step); self.sigma_xy[0].append(sigma_xy_current_time_step)  

        
        print("max point_1:---------------------------------------------")
        print(ppp)
        print(s_max)
        print("---------------------------------------------")  

        
        mp = model["IgaModelPart.Body2.StructuralAnalysisDomain"]
        s_max = 0
        sigma_x_current_time_step = []
        sigma_y_current_time_step = []
        sigma_VM_current_time_step = []
        sigma_xy_current_time_step = []
        for elem in mp.Elements:
            geom = elem.GetGeometry()
            N = geom.ShapeFunctionsValues()
            weight = elem.GetValue(KratosMultiphysics.INTEGRATION_WEIGHT)
            
            # GET GAUSS POINT COORD
            stress_xx = elem.GetValue(KratosMultiphysics.IgaApplication.CAUCHY_STRESS_XX)
            stress_yy = elem.GetValue(KratosMultiphysics.IgaApplication.CAUCHY_STRESS_YY)
            stress_xy = elem.GetValue(KratosMultiphysics.IgaApplication.CAUCHY_STRESS_XY)
            sigma_vm = elem.GetValue(KratosMultiphysics.YIELD_STRESS)
            sigma_x_current_time_step.append(stress_xx);  sigma_y_current_time_step.append(stress_yy)   
            sigma_VM_current_time_step.append(sigma_vm); sigma_xy_current_time_step.append(stress_xy)  

        # for cond in mp.Conditions:
        #     geom = cond.GetGeometry()
        #     N = geom.ShapeFunctionsValues()
        #     weight = cond.GetValue(KratosMultiphysics.INTEGRATION_WEIGHT)
            
        #     # GET GAUSS POINT COORD
        #     stress_xx = cond.GetValue(KratosMultiphysics.IgaApplication.CAUCHY_STRESS_XX)
        #     stress_yy = cond.GetValue(KratosMultiphysics.IgaApplication.CAUCHY_STRESS_YY)
        #     stress_xy = cond.GetValue(KratosMultiphysics.IgaApplication.CAUCHY_STRESS_XY)
        #     sigma_vm = cond.GetValue(KratosMultiphysics.YIELD_STRESS)
        #     sigma_x_current_time_step.append(stress_xx);  sigma_y_current_time_step.append(stress_yy)   
        #     sigma_VM_current_time_step.append(sigma_vm); sigma_xy_current_time_step.append(stress_xy)    

        #     if (sigma_vm > s_max): 
        #         gauss_pos = np.zeros(2)
        #         i = 0
        #         for n in geom:
        #             gauss_pos[0] += N[0, i] * n.X
        #             gauss_pos[1] += N[0, i] * n.Y
        #             i += 1

        #         s_max = sigma_vm
        #         ppp = gauss_pos

        print("max point_2:---------------------------------------------")
        print(ppp)
        print(s_max)
        print("---------------------------------------------")

        self.sigma_x[1].append(sigma_x_current_time_step);  self.sigma_y[1].append(sigma_y_current_time_step)   
        self.sigma_VM[1].append(sigma_VM_current_time_step);  self.sigma_xy[1].append(sigma_xy_current_time_step)     

        #%++++++++++++++++++
        # mp = model["IgaModelPart.Body1.Load_Support_outer_1"]
        mp = model["IgaModelPart.Body1.SBM_Support_outer.top"]
        step_reaction_force_x = 0.0
        step_reaction_force_y = 0.0
        for cond in mp.Conditions:
            weight = cond.GetValue(KratosMultiphysics.INTEGRATION_WEIGHT)
            sigma_n = cond.GetValue(KratosMultiphysics.NORMAL_STRESS)
            
            step_reaction_force_x += weight * sigma_n[0]
            step_reaction_force_y += weight * sigma_n[1]
            
        self.reaction_force_body1[0].append(step_reaction_force_x)
        self.reaction_force_body1[1].append(step_reaction_force_y)
        
        #%++++++++++++++++++
        mp = model["IgaModelPart.Body2.SBM_Support_outer.bottom"]
        
        step_reaction_force_x = 0.0
        step_reaction_force_y = 0.0
        for cond in mp.Conditions:
            weight = cond.GetValue(KratosMultiphysics.INTEGRATION_WEIGHT)
            sigma_n = cond.GetValue(KratosMultiphysics.NORMAL_STRESS)
            step_reaction_force_x += weight * sigma_n[0]
            step_reaction_force_y += weight * sigma_n[1]

        self.reaction_force_body2[0].append(step_reaction_force_x)
        self.reaction_force_body2[1].append(step_reaction_force_y)

        
        mp = model["IgaModelPart.Body1.SBM_Support_outer"]
        
        step_contact_force_x_body1 = 0.0
        step_contact_force_y_body1 = 0.0
        #
        step_contact_force_x_body2 = 0.0
        step_contact_force_y_body2 = 0.0
        for cond in mp.Conditions:
            geom = cond.GetGeometry()
            N = geom.ShapeFunctionsValues()

            weight = cond.GetValue(KratosMultiphysics.INTEGRATION_WEIGHT)
            sigma_n = cond.GetValue(KratosMultiphysics.NORMAL_STRESS)
            
            step_contact_force_x_body1 += weight * sigma_n[0]
            step_contact_force_y_body1 += weight * sigma_n[1]
            gauss_pos = np.zeros(2)
            i = 0
            for n in geom:
                gauss_pos[0] += N[0, i] * n.X
                gauss_pos[1] += N[0, i] * n.Y
                i += 1

            t_vec.append(gauss_pos[1]); 
            VM_vec.append(cond.GetValue(KratosMultiphysics.YIELD_STRESS)); s_xx_vec.append(sigma_n[0])

        # fig, ax = plt.subplots()
        # ax.plot(t_vec, VM_vec, "-b") 
        # ax.plot(t_vec, s_xx_vec, "-g") 
        # ax.legend("VM_stress", "S_xx")
        # plt.show()

            
        self.contact_force_body1[0].append(step_contact_force_x_body1)
        self.contact_force_body1[1].append(step_contact_force_y_body1)
        
        ## BODY 2 CONTACT INTERFACE
        
        mp = model["IgaModelPart.Body2.SBM_Support_outer"]
        
        step_contact_force_x_body2 = 0.0
        step_contact_force_y_body2 = 0.0
        for cond in mp.Conditions:
            weight = cond.GetValue(KratosMultiphysics.INTEGRATION_WEIGHT)
            sigma_n = cond.GetValue(KratosMultiphysics.NORMAL_STRESS)
            
            step_contact_force_x_body2 += weight * sigma_n[0]
            step_contact_force_y_body2 += weight * sigma_n[1]
        
        self.contact_force_body2[0].append(step_contact_force_x_body2)
        self.contact_force_body2[1].append(step_contact_force_y_body2)
        
        # mp = model["IgaModelPart.ContactInterface.Contact1"]
        
        # step_contact_force_x_body1 = 0.0
        # step_contact_force_y_body1 = 0.0
        # step_contact_force_x_body2 = 0.0
        # step_contact_force_y_body2 = 0.0
        # for cond in mp.Conditions:
        #     weight = cond.GetValue(KratosMultiphysics.INTEGRATION_WEIGHT)

        #     sigma_master = cond.GetValue(KratosMultiphysics.IgaApplication.STRESS_MASTER)
        #     sigma_slave = cond.GetValue(KratosMultiphysics.IgaApplication.STRESS_SLAVE)

        #     normal_master = cond.GetValue(KratosMultiphysics.IgaApplication.NORMAL_MASTER)
        #     normal_slave = cond.GetValue(KratosMultiphysics.IgaApplication.NORMAL_SLAVE)

        #     ###
        #     stress_master = np.zeros(2)
        #     stress_master[0] = sigma_master[0]*normal_master[0] + sigma_master[2]*normal_master[1]
        #     stress_master[1] = sigma_master[2]*normal_master[0] + sigma_master[1]*normal_master[1]

        #     stress_slave = np.zeros(2)
        #     stress_slave[0] = sigma_slave[0]*normal_slave[0] + sigma_slave[2]*normal_slave[1]
        #     stress_slave[1] = sigma_slave[2]*normal_slave[0] + sigma_slave[1]*normal_slave[1]
        #     ###
            
        #     step_contact_force_x_body1 += weight * stress_master[0]
        #     step_contact_force_y_body1 += weight * stress_master[1]

        #     step_contact_force_x_body2 += weight * stress_slave[0]
        #     step_contact_force_y_body2 += weight * stress_slave[1]
        
        # self.contact_force_body1[0].append(step_contact_force_x_body1)
        # self.contact_force_body1[1].append(step_contact_force_y_body1)

        # self.contact_force_body2[0].append(step_contact_force_x_body2)
        # self.contact_force_body2[1].append(step_contact_force_y_body2)





def GetModelInfo(mp, i_body):
    for elem in mp.Elements:
        geom = elem.GetGeometry()

        N = geom.ShapeFunctionsValues()

        gauss_pos = np.zeros(2)
        
        weight = elem.GetValue(KratosMultiphysics.INTEGRATION_WEIGHT)
        

        # GET GAUSS POINT COORD
        i = 0
        for n in geom:
            gauss_pos[0] += N[0, i] * n.X
            gauss_pos[1] += N[0, i] * n.Y
            i += 1

        x_coord[i_body].append(gauss_pos[0])
        y_coord[i_body].append(gauss_pos[1])
        
        #-------------------
        # GET SOLUTION STEP VALUES
        for i_step in range(0, n_time_steps):
            curr_disp_x = 0
            curr_disp_y = 0
            i = 0
            for n in geom:
                curr_disp_x += N[0, i] * n.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, n_time_steps-1-i_step)
                curr_disp_y += N[0, i] * n.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, n_time_steps-1-i_step)
                i += 1
            
            disp_x[i_body][i_step].append(curr_disp_x)
            disp_y[i_body][i_step].append(curr_disp_y)

            
    # INCLUDE SUPPORT
    
    for cond in mp.Conditions:
        geom = cond.GetGeometry()

        if (len(geom) == 0): continue

        N = geom.ShapeFunctionsValues()

        gauss_pos = np.zeros(2)

        # GET GAUSS POINT COORD
        i = 0
        for n in geom:
            gauss_pos[0] += N[0, i] * n.X
            gauss_pos[1] += N[0, i] * n.Y
            i += 1

        x_coord[i_body].append(gauss_pos[0])
        y_coord[i_body].append(gauss_pos[1])
        
        #-------------------
        # GET SOLUTION STEP VALUES
        for i_step in range(0, n_time_steps):
            curr_disp_x = 0
            curr_disp_y = 0
            i = 0
            for n in geom:
                curr_disp_x += N[0, i] * n.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, n_time_steps-1-i_step)
                curr_disp_y += N[0, i] * n.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, n_time_steps-1-i_step)
                i += 1
            
            disp_x[i_body][i_step].append(curr_disp_x)
            disp_y[i_body][i_step].append(curr_disp_y)
        
        plt.plot(gauss_pos[0], gauss_pos[1], "k*", markersize = 1)
            
        
###########################################################################33
##############################################################################


if __name__ == "__main__":

    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    
    if os.path.exists("physical_space.txt"):
        os.remove("physical_space.txt")
    if os.path.exists("DISPLACEMENT_X.txt"):
        os.remove("DISPLACEMENT_X.txt")
    if os.path.exists("txt_files/output_results_GPs.txt"):
        os.remove("txt_files/output_results_GPs.txt")
        
        
    if os.path.exists("txt_files/output_results_GPs_master.txt"):
        os.remove("txt_files/output_results_GPs_master.txt")
    if os.path.exists("txt_files/output_results_GPs_slave.txt"):
        os.remove("txt_files/output_results_GPs_slave.txt")
    # if os.path.exists("txt_files/gap_output.txt"):
    #     os.remove("txt_files/gap_output.txt")
    if os.path.exists("txt_files/projection_contact_slave.txt"):
        os.remove("txt_files/projection_contact_slave.txt")
    if os.path.exists("txt_files/projection_contact_master.txt"):
        os.remove("txt_files/projection_contact_master.txt")
    if os.path.exists("txt_files/error_projection_contact_slave.txt"):
        os.remove("txt_files/error_projection_contact_slave.txt")
    if os.path.exists("txt_files/error_projection_contact_master.txt"):
        os.remove("txt_files/error_projection_contact_master.txt")
        
    if os.path.exists("txt_files/gap_output_master.txt"):
        os.remove("txt_files/gap_output_master.txt")
    if os.path.exists("txt_files/gap_output_slave.txt"):
        os.remove("txt_files/gap_output_slave.txt")
    if os.path.exists("txt_files/contact_displacement.txt"):
        os.remove("txt_files/contact_displacement.txt")
        
        
    if os.path.exists("txt_files/new_normal.txt"):
        os.remove("txt_files/new_normal.txt")
        
                
    if os.path.exists("txt_files/Gauss_Point_coordinates.txt"):
        os.remove("txt_files/Gauss_Point_coordinates.txt")
    if os.path.exists("txt_files/Gauss_Point_Contact_coordinates.txt"):
        os.remove("txt_files/Gauss_Point_Contact_coordinates.txt")
    if os.path.exists("txt_files/Id_active_control_points.txt"):
        os.remove("txt_files/Id_active_control_points.txt")
    if os.path.exists("txt_files/Id_active_control_points_condition.txt"):
        os.remove("txt_files/Id_active_control_points_condition.txt")
    if os.path.exists("txt_files/Id_active_control_points_condition_neumann.txt"):
        os.remove("txt_files/Id_active_control_points_condition_neumann.txt")
    if os.path.exists("txt_files/Projection_Coordinates.txt"):
        os.remove("txt_files/Projection_Coordinates.txt")
    if os.path.exists("txt_files/Contact_Projection_Coordinates.txt"):
        os.remove("txt_files/Contact_Projection_Coordinates.txt")
    if os.path.exists("txt_files/Snake_coordinates.txt"):
        os.remove("txt_files/Snake_coordinates.txt")
    if os.path.exists("txt_files/Snake_coordinates_outer.txt"):
        os.remove("txt_files/Snake_coordinates_outer.txt")
    if os.path.exists("txt_files/boundary_GPs.txt"):
        os.remove("txt_files/boundary_GPs.txt")

    model = KratosMultiphysics.Model()
    # simulation = StructuralMechanicsAnalysis(model, parameters)
    
    simulation = CustomAnalysisStage(model, parameters)
    simulation.Run()

    
    fig, ax = plt.subplots()
    index_param_space= 2
    # body 1 parameter space
    n_knot_span_u = parameters["modelers"][index_param_space]["Parameters"]["number_of_knot_spans"][0].GetInt()
    n_knot_span_v = parameters["modelers"][index_param_space]["Parameters"]["number_of_knot_spans"][1].GetInt()
    initial_u = parameters["modelers"][index_param_space]["Parameters"]["lower_point_uvw"][0].GetDouble()
    initial_v = parameters["modelers"][index_param_space]["Parameters"]["lower_point_uvw"][1].GetDouble()
    final_u = parameters["modelers"][index_param_space]["Parameters"]["upper_point_uvw"][0].GetDouble()
    final_v = parameters["modelers"][index_param_space]["Parameters"]["upper_point_uvw"][1].GetDouble()

    knots_u_1 = [initial_u]
    knots_v_1 = [initial_v]

    for j in range(1, n_knot_span_u):
        knots_u_1.append(initial_u + (final_u-initial_u) / (n_knot_span_u) * j)
    for j in range(1, n_knot_span_v):
        knots_v_1.append(initial_v + (final_v-initial_v) / (n_knot_span_v) * j)
    knots_u_1.append(final_u)
    knots_v_1.append(final_v)

    # print(knots_u_1)

    # Drawing lines for knots
    # for u in knots_u_1:
    #     ax.axvline(x=u, color='red', linestyle='--', linewidth=1.0, zorder=4)
    # for v in knots_v_1:
    #     ax.axhline(y=v, color='red', linestyle='--', linewidth=1.0, zorder=4)

    # # Enhancing grid visibility
    # ax.grid(True, linestyle='--', linewidth=0.5)

    # Setting the aspect ratio
    ax.set_aspect('equal', adjustable='box')
    #--------------------------------------------------------------------
    # body 2 parameter space
    index_param_space += 1
    n_knot_span_u = parameters["modelers"][index_param_space]["Parameters"]["number_of_knot_spans"][0].GetInt()
    n_knot_span_v = parameters["modelers"][index_param_space]["Parameters"]["number_of_knot_spans"][1].GetInt()
    initial_u = parameters["modelers"][index_param_space]["Parameters"]["lower_point_uvw"][0].GetDouble()
    initial_v = parameters["modelers"][index_param_space]["Parameters"]["lower_point_uvw"][1].GetDouble()
    final_u = parameters["modelers"][index_param_space]["Parameters"]["upper_point_uvw"][0].GetDouble()
    final_v = parameters["modelers"][index_param_space]["Parameters"]["upper_point_uvw"][1].GetDouble()

    knots_u_2 = [initial_u]
    knots_v_2 = [initial_v]

    for j in range(1, n_knot_span_u):
        knots_u_2.append(initial_u + (final_u-initial_u) / (n_knot_span_u) * j)
    for j in range(1, n_knot_span_v):
        knots_v_2.append(initial_v + (final_v-initial_v) / (n_knot_span_v) * j)
    knots_u_2.append(final_u)
    knots_v_2.append(final_v)

    # Drawing lines for knots
    # for u in knots_u_2:
    #     ax.axvline(x=u, color='red', linestyle='--', linewidth=1.0, zorder=4)
    # for v in knots_v_2:
    #     ax.axhline(y=v, color='red', linestyle='--', linewidth=1.0, zorder=4)

    # Enhancing grid visibility
    ax.grid(True, linestyle='--', linewidth=0.5)

    # Setting the aspect ratio
    ax.set_aspect('equal', adjustable='box')

    #--------------------------------------
    ax.set_xlim(min(min(knots_u_1), min(knots_u_2)), max(max(knots_u_1), max(knots_u_2)))
    ax.set_ylim(min(min(knots_v_1), min(knots_v_2)), max(max(knots_v_1), max(knots_v_2)))
    #-------------------------------------------------------------------------------------
    # DRAW SURROGATE BOUNDARY 
    # BODY_1
    mp = model["IgaModelPart.Body1.surrogate_outer"]
    
    x_surrogate_vec_out = []
    y_surrogate_vec_out = []
    count = 0
    for cond in mp.Conditions:
        geom = cond.GetGeometry()
        
        if count == 0:
            x_surrogate_vec_out.append(geom[0].X)
            y_surrogate_vec_out.append(geom[0].Y)
        
        x_surrogate_vec_out.append(geom[1].X)
        y_surrogate_vec_out.append(geom[1].Y)

    x_surrogate_vec_out.append(x_surrogate_vec_out[0])
    y_surrogate_vec_out.append(y_surrogate_vec_out[0])

    plt.plot(x_surrogate_vec_out, y_surrogate_vec_out, "b-")
    # BODY_2
    mp = model["IgaModelPart.Body2.surrogate_outer"]
    
    x_surrogate_vec_out = []
    y_surrogate_vec_out = []
    count = 0
    for cond in mp.Conditions:
        geom = cond.GetGeometry()
        
        if count == 0:
            x_surrogate_vec_out.append(geom[0].X)
            y_surrogate_vec_out.append(geom[0].Y)
        
        x_surrogate_vec_out.append(geom[1].X)
        y_surrogate_vec_out.append(geom[1].Y)

    x_surrogate_vec_out.append(x_surrogate_vec_out[0])
    y_surrogate_vec_out.append(y_surrogate_vec_out[0])

    plt.plot(x_surrogate_vec_out, y_surrogate_vec_out, "b-")
    
    #--------------------------------------------------------------
    # DRAW SKIN BOUNDARY 
    # BODY_1
    mp = model["skin_Body1.outer"]
    
    x_skin_vec_out = []
    y_skin_vec_out = []
    count = 0
    for cond in mp.Conditions:
        geom = cond.GetGeometry()
        
        if count == 0:
            x_skin_vec_out.append(geom[0].X)
            y_skin_vec_out.append(geom[0].Y)
        
        x_skin_vec_out.append(geom[1].X)
        y_skin_vec_out.append(geom[1].Y)

    x_skin_vec_out.append(x_skin_vec_out[0])
    y_skin_vec_out.append(y_skin_vec_out[0])

    plt.plot(x_skin_vec_out, y_skin_vec_out, "g-", markersize=2)
    # BODY_2
    mp = model["skin_Body2.outer"]
    
    x_skin_vec_out = []
    y_skin_vec_out = []
    count = 0
    for cond in mp.Conditions:
        geom = cond.GetGeometry()
        
        if count == 0:
            x_skin_vec_out.append(geom[0].X)
            y_skin_vec_out.append(geom[0].Y)
        
        x_skin_vec_out.append(geom[1].X)
        y_skin_vec_out.append(geom[1].Y)

    x_skin_vec_out.append(x_skin_vec_out[0])
    y_skin_vec_out.append(y_skin_vec_out[0])

    plt.plot(x_skin_vec_out, y_skin_vec_out, "g-", markersize=2)

    # with open("txt_files/Contact_Projection_Coordinates.txt", "r") as file:
    #     lines = file.readlines()
    #     for line in lines[:]:
    #         values = line.strip().split()

    #         x_coord_master_slave = [float(values[0]), float(values[2])]
    #         y_coord_master_slave = [float(values[1]), float(values[3])]
    #         plt.plot(x_coord_master_slave, y_coord_master_slave, "b--", markersize=2)
    
    # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    # # fig = plt.figure()
    
    # # fig3, ax3 = plt.subplots()

    n_objects = 2

    disp_x = []  # Using lists instead of np.zeros
    disp_y = []
    x_coord = []
    y_coord = []
    
    gp_contact_x = []
    gp_contact_y = []
    gp_contact_deformed_x = []
    gp_contact_deformed_y = []
    
    #++++++++++++++++++++++++++++++++
    #++++++++++++++++++++++++++++++++
    for i_body in range(0, n_objects):
        disp_x.append([])
        disp_y.append([])
        x_coord.append([])
        y_coord.append([])
        
        gp_contact_x.append([])
        gp_contact_y.append([])
        gp_contact_deformed_x.append([])
        gp_contact_deformed_y.append([])
    
    # #++++++++++++++++++++++++++++++++
    # #++++++++++++++++++++++++++++++++
    
    
    x = sp.symbols('x')
    y = sp.symbols('y')
    t = sp.symbols('t')
    # internal_load_x = -1000*(1+0.3)/(1-0.3*0.3)*sp.cos(x)*sp.sinh(y)
    # internal_load_y = -1000*(1+0.3)/(1-0.3*0.3)*sp.sin(x)*sp.cosh(y)
    
    # internal_load_x = 0.08*t
    # internal_load_y = 0
    
    internal_load_x = 0
    internal_load_y = 0
    internal_load_handle = [sp.lambdify((x, y,t), internal_load_x), sp.lambdify((x, y,t), internal_load_y)]

    n_time_steps = simulation.nSTEP
    
    external_force = [np.zeros(n_time_steps), np.zeros(n_time_steps)]
    
    
    for i_step in range(n_time_steps-1, -1, -1):
        
        for i_body in range(0, n_objects):
            disp_x[i_body].append([])
            disp_y[i_body].append([])
        
    
    on_deformed = True
    mp_0 = model["IgaModelPart.Body1.StructuralAnalysisDomain"]
    GetModelInfo(mp_0, 0)

    # mp_0 = model["IgaModelPart.Body1.SBM_Support_outer"]
    # GetModelInfo(mp_0, 0)

    mp_1 = model["IgaModelPart.Body2.StructuralAnalysisDomain"]
    GetModelInfo(mp_1, 1)

    # mp_1 = model["IgaModelPart.Body2.SBM_Support_outer"]
    # GetModelInfo(mp_1, 1)

    # ax.scatter(x_coord[0], y_coord[0], marker='o', color='blue', s=2, zorder=2)
    # ax.scatter(x_coord[1], y_coord[1], marker='o', color='blue', s=2, zorder=2)


    if os.path.exists("txt_files/Projection_Coordinates.txt"):
        with open('txt_files/Projection_Coordinates.txt', 'r') as file:
            lines = file.readlines()
        x_values = []
        y_values = []
        x_gauss = []
        y_gauss = []
        for line in lines:
            parts = line.split()  # Dividi la riga in parti utilizzando lo spazio come separatore
            if len(parts) == 4:  # Assicurati che ci siano due valori su ogni riga
                x_values.append(float(parts[0]))
                y_values.append(float(parts[1]))
                x_gauss.append(float(parts[2]))
                y_gauss.append(float(parts[3]))
        for i in range(len(x_values)):
            plt.plot([x_values[i], x_gauss[i]], [y_values[i], y_gauss[i]], color='green')


        plt.scatter(x_values, y_values, marker='o', color='green')
    
    #---------------------------------
    # CONTACT
    #*********************************
    # mp = model["IgaModelPart.ContactInterface"]
    
    # points_for_gp_connections = []
    # points_for_gp_connections_deformed = []
    # for cond in mp.Conditions:
        
    #     temp_vec = []
    #     temp_vec_deformed = []
    #     for i in range(0,2):
    #         geom = cond.GetGeometry().GetGeometryPart(i)

    #         if (len(geom) == 0): continue

    #         N = geom.ShapeFunctionsValues()

    #         gauss_pos = np.zeros(2)

    #         # GET GAUSS POINT COORD
    #         i = 0
    #         for n in geom:
    #             gauss_pos[0] += N[0, i] * n.X
    #             gauss_pos[1] += N[0, i] * n.Y
    #             i += 1
                
    #         # print("GAUSS_POS: ", gauss_pos)
                

    #         if gauss_pos[0] > 4.98 + 1e-13:
    #             x_coord[1].append(gauss_pos[0])
    #             y_coord[1].append(gauss_pos[1])
                
    #             #-------------------
    #             # GET SOLUTION STEP VALUES
    #             for i_step in range(0, n_time_steps):
    #                 curr_disp_x = 0
    #                 curr_disp_y = 0
    #                 i = 0
    #                 for n in geom:
    #                     curr_disp_x += N[0, i] * n.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, n_time_steps-1-i_step)
    #                     curr_disp_y += N[0, i] * n.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, n_time_steps-1-i_step)
    #                     i += 1
                    
    #                 disp_x[1][i_step].append(curr_disp_x)
    #                 disp_y[1][i_step].append(curr_disp_y)
                
    #             gp_contact_x[1].append(gauss_pos[0])
    #             gp_contact_y[1].append(gauss_pos[1])
    #             gp_contact_deformed_x[1].append(gauss_pos[0]+curr_disp_x)
    #             gp_contact_deformed_y[1].append(gauss_pos[1]+curr_disp_y)
                        
    #         else:
    #             x_coord[0].append(gauss_pos[0])
    #             y_coord[0].append(gauss_pos[1])
                
    #             #-------------------
    #             # GET SOLUTION STEP VALUES
    #             for i_step in range(0, n_time_steps):
    #                 curr_disp_x = 0
    #                 curr_disp_y = 0
    #                 i = 0
    #                 for n in geom:
    #                     curr_disp_x += N[0, i] * n.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, n_time_steps-1-i_step)
    #                     curr_disp_y += N[0, i] * n.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, n_time_steps-1-i_step)
    #                     i += 1
                    
    #                 disp_x[0][i_step].append(curr_disp_x)
    #                 disp_y[0][i_step].append(curr_disp_y)
                
    #             gp_contact_x[0].append(gauss_pos[0])
    #             gp_contact_y[0].append(gauss_pos[1])
    #             gp_contact_deformed_x[0].append(gauss_pos[0]+curr_disp_x)
    #             gp_contact_deformed_y[0].append(gauss_pos[1]+curr_disp_y)
            
    #         # temp_vec.append(gauss_pos[0])
    #         # temp_vec.append(gauss_pos[1])
    #         # temp_vec_deformed.append(gauss_pos[0]+curr_disp_x)
    #         # temp_vec_deformed.append(gauss_pos[1]+curr_disp_y)
            
    #     points_for_gp_connections.append(temp_vec)
    #     points_for_gp_connections_deformed.append(temp_vec_deformed)
    #-------------------------------------------------
            
    # Labeling axes
    # ax.set_xlabel("x")
    # ax.set_ylabel("y")

    # # Setting limits based on the knots
    # ax.set_xlim(-5.2, 7.0)
    # ax.set_ylim(-5.2, 6.2)
    
    # # Setting the aspect ratio
    # ax.set_aspect('equal', adjustable='box')
   
        
    ################################################################################33

    # mp = model["IgaModelPart.Body1.Load_Support_outer"]
    # sum_tot = 0
    # for cond in mp.Conditions:
            
    #     geom = cond.GetGeometry()
    #     N = geom.ShapeFunctionsValues()

    #     gauss_pos = np.zeros(2)
    #     i = 0
    #     for n in geom:
    #         gauss_pos[0] += N[0, i] * n.X
    #         gauss_pos[1] += N[0, i] * n.Y
    #         i += 1
        
    #     print(gauss_pos)

    #     weight = cond.GetValue(KratosMultiphysics.INTEGRATION_WEIGHT)

    #     sum_tot += weight
    # print(sum_tot)
    # emomf
    ##################################################################################
    triangulation_1 = mtri.Triangulation(x_coord[0], y_coord[0])
    triangulation_2 = mtri.Triangulation(x_coord[1], y_coord[1])
    
    triangulations = [triangulation_1, triangulation_2]

    # Funzione per ottenere la triangolazione
    def GetTriangulation(curr_disp_x, curr_disp_y, i_body):
        if on_deformed:
            x_coord_deformed = list(map(lambda a, b: a + b, x_coord[i_body], curr_disp_x))
            y_coord_deformed = list(map(lambda a, b: a + b, y_coord[i_body], curr_disp_y))
            # new_triangulation = mtri.Triangulation(x_coord_deformed, y_coord_deformed)
            
            areas = []
            for triangle in triangulations[i_body].triangles:
                x1, y1 = x_coord_deformed[triangle[0]], y_coord_deformed[triangle[0]]
                x2, y2 = x_coord_deformed[triangle[1]], y_coord_deformed[triangle[1]]
                x3, y3 = x_coord_deformed[triangle[2]], y_coord_deformed[triangle[2]]
                # Calculate the area of the triangle using the shoelace formula
                area = 0.5 * abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
                areas.append(area)
            average_area = sum(areas) / len(areas)
            new_triangles = []
            for i in range(len(areas)):
                if True: #areas[i] <= average_area*2.5:
                    new_triangles.append(triangulations[i_body].triangles[i])
            
            final_triangulation = mtri.Triangulation(x_coord_deformed, y_coord_deformed, new_triangles)
            return final_triangulation
        else:
            return triangulations[i_body]
        

    deformed_triangulations = [[], []]
    for i_step in range(n_time_steps):
         
         for i_body in range(0,2):
            new_triangulation = GetTriangulation(disp_x[i_body][i_step], disp_y[i_body][i_step], i_body)
            deformed_triangulations[i_body].append(new_triangulation)


    variable_to_print = [disp_x, disp_y]
    
    fig, ax = plt.subplots(1,2 ,figsize=(10, 8))
    plt.subplots_adjust(bottom=0.25)  # Lascia spazio per lo slider    
    cbar_list, lines = PlotVariable(ax, fig, n_objects, n_time_steps, variable_to_print, deformed_triangulations)

    # Creare lo slider
    ax_slider = plt.axes([0.25, 0.1, 0.5, 0.03])  # Posizionare lo slider
    slider = Slider(ax_slider, 'Indice', 0, n_time_steps - 1, valinit=0, valstep=1)
    slider.on_changed(lambda val: update_plot(int(val), variable_to_print, cbar_list, lines, fig, ax))

    update_plot(0, variable_to_print, cbar_list, lines, fig, ax)

    #################
    variable_to_print1 = [simulation.sigma_x, simulation.sigma_y]
    
    fig1, ax1 = plt.subplots(1,2 ,figsize=(10, 8))
    plt.subplots_adjust(bottom=0.25)  # Lascia spazio per lo slider    
    cbar_list1, lines1 = PlotVariable(ax1, fig1, n_objects, n_time_steps, variable_to_print1, deformed_triangulations)

    # Creare lo slider
    ax_slider1 = plt.axes([0.25, 0.1, 0.5, 0.03])  # Posizionare lo slider
    slider1 = Slider(ax_slider1, 'Indice', 0, n_time_steps - 1, valinit=0, valstep=1)
    slider1.on_changed(lambda val: update_plot(int(val), variable_to_print1, cbar_list1, lines1, fig1, ax1))
    

    update_plot(0, variable_to_print1, cbar_list1, lines1, fig1, ax1)

    ################
    variable_to_print2 = [simulation.sigma_xy, simulation.sigma_VM]
    
    fig2, ax2 = plt.subplots(1,2 ,figsize=(10, 8))
    plt.subplots_adjust(bottom=0.25)  # Lascia spazio per lo slider    
    cbar_list2, lines2 = PlotVariable(ax2, fig2, n_objects, n_time_steps, variable_to_print2, deformed_triangulations)

    # Creare lo slider
    ax_slider2 = plt.axes([0.25, 0.1, 0.5, 0.03])  # Posizionare lo slider
    slider2 = Slider(ax_slider2, 'Indice', 0, n_time_steps - 1, valinit=0, valstep=1)
    slider2.on_changed(lambda val: update_plot(int(val), variable_to_print2, cbar_list2, lines2, fig2, ax2))
    

    update_plot(0, variable_to_print2, cbar_list2, lines2, fig2, ax2)

    
    print(external_force)
    print(simulation.reaction_force_body1[1])
    reaction_1 = simulation.reaction_force_body1[1]
    reaction_2 = simulation.reaction_force_body2[1]
    
    contact_force_1 = simulation.contact_force_body1[1]
    contact_force_2 = simulation.contact_force_body2[1]
    
    print(contact_force_2)
    
    print(list(map(lambda a, b: a + b, simulation.reaction_force_body1[1], simulation.reaction_force_body2[1])))
    
    
    for i in range(0, n_time_steps):
        print("STEP %d -- TIME %f" % (i, simulation.time_step[-n_time_steps + i]))
        print("     EXT_FORCE: %7.4f    REACTION_1: %7.4f   REACTION_2: %7.4f" % (external_force[0][i], reaction_1[i], reaction_2[i]))
        
        print("     CONTACT FORCE_1: %7.4e    CONTACT FORCE_2: %7.4e " % (contact_force_1[i], contact_force_2[i]))
        
        print(" --> ERROR: %e" % (np.abs(external_force[0][i] + (reaction_1[i] + reaction_2[i]))))
        print(" --> ERROR_2: %e" % (np.abs(external_force[0][i] + (reaction_1[i] + contact_force_1[i]))))
        
        print("--------------------------------------")


    
    

    # # Directory per salvare le immagini temporanee
    # temp_dir = "temp_images"
    # os.makedirs(temp_dir, exist_ok=True)

    # # Lista per salvare i percorsi dei file delle immagini temporanee
    # image_files = []

    # # Salva ogni frame come immagine
    # for i_step in range(0, n_time_steps):
    #     update_plot(i_step)
        
    #     # Salva l'immagine per lo step corrente
    #     image_path = f"{temp_dir}/frame_{i_step}.png"
    #     plt.savefig(image_path)
    #     image_files.append(image_path)

    # # Creazione della GIF
    # with imageio.get_writer('Contact Circle-Wall.gif', mode='I', duration=0.2) as writer:
    #     for filename in image_files:
    #         image = imageio.imread(filename)
    #         writer.append_data(image)

    # # Opzionalmente, rimuovi le immagini temporanee
    # for filename in image_files:
    #     os.remove(filename)

    # # Rimuovi la cartella temporanea
    # os.rmdir(temp_dir)


    plt.show()
