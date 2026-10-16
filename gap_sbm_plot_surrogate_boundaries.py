import KratosMultiphysics
import KratosMultiphysics.IgaApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as patches

from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import os

if __name__ == "__main__":

    if os.path.exists("centers.txt"):
        os.remove("centers.txt")
    if os.path.exists("centers_load.txt"):
        os.remove("centers_load.txt")
    if os.path.exists("txt_files/Id_active_control_points.txt"):
        os.remove("txt_files/Id_active_control_points.txt")
    if os.path.exists("txt_files/Id_active_control_points_condition.txt"):
        os.remove("txt_files/Id_active_control_points_condition.txt")

    with open("ProjectParametersSolidGap.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    
    simulation = StructuralMechanicsAnalysis(model, parameters)
    simulation.Initialize()
    simulation.InitializeSolutionStep()

    index_param_space= 0
    while not (parameters["modelers"][index_param_space]["Parameters"].Has("number_of_knot_spans")):
        index_param_space += 1
    
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

    # Create figure and axes
    # fig, ax = plt.subplots()
    fig = plt.gcf()      # get current figure (se non c’è, ne crea una)
    ax  = fig.gca() 
    ax.set_aspect('equal', adjustable='box')
    ax.grid(False)  
    label_used = set()

    # === PARAMETER SPACE GRID: BODY 1 ===
    for u in knots_u_1:
        ax.plot([u, u], [knots_v_1[0], knots_v_1[-1]], linestyle=':', color='gray', linewidth=1.5, label='Parameter space' if u == knots_u_1[0] else "")
    for v in knots_v_1:
        ax.plot([knots_u_1[0], knots_u_1[-1]], [v, v], linestyle=':', color='gray', linewidth=1.5)
    
    label_used.add("Parameter space")
    
    for side in ["outer", "inner"]:
        surrogate_part_name = f"IgaModelPart.surrogate_{side}"
        skin_part_name = f"skin_model_part.{side}"
                
        if model.HasModelPart(surrogate_part_name):
            mp = model.GetModelPart(surrogate_part_name)
            
            mskin = model.GetModelPart(skin_part_name)
            
            if len(mp.Conditions) > 0:
                for cond in mp.Conditions:
                    x_surr, y_surr = [], []
                    geom = cond.GetGeometry()
                    x_surr.append(geom[0].X)
                    y_surr.append(geom[0].Y)
                    x_surr.append(geom[1].X)
                    y_surr.append(geom[1].Y)
                    ax.plot(x_surr, y_surr, color='darkred', linewidth=5, label="Surrogate boundary")
                                    
                    projection_node_id = geom[0].GetValue(KratosMultiphysics.IgaApplication.PROJECTION_NODE_ID)
                    x1 = geom[0].X
                    y1 = geom[0].Y
    
                    x_projection1 = mskin.GetNode(projection_node_id).X
                    y_projection1 = mskin.GetNode(projection_node_id).Y
                    ax.plot([x1, x_projection1], [y1, y_projection1], "k--", linewidth=2)
                    
                    projection_node_id = geom[1].GetValue(KratosMultiphysics.IgaApplication.PROJECTION_NODE_ID)
                    x2 = geom[1].X
                    y2 = geom[1].Y
    
                    x_projection2 = mskin.GetNode(projection_node_id).X
                    y_projection2 = mskin.GetNode(projection_node_id).Y
                    ax.plot([x2, x_projection2], [y2, y_projection2], "k--", linewidth=2)
                    
                    # ax.plot([x_projection1, x_projection2], [y_projection1, y_projection2], "b--", linewidth=1)
                        
                        
                        
                label_used.add("Surrogate boundary")                                                
        
        if model.HasModelPart(skin_part_name):
            mp = model.GetModelPart(skin_part_name)
            if len(mp.Conditions) > 0:
                x_skin, y_skin = [], []
                for i, cond in enumerate(mp.Conditions):
                    geom = cond.GetGeometry()
                    if i == 0:
                        x_skin.append(geom[0].X)
                        y_skin.append(geom[0].Y)
                    x_skin.append(geom[1].X)
                    y_skin.append(geom[1].Y)
                x_skin.append(x_skin[0])
                y_skin.append(y_skin[0])
                ax.plot(x_skin, y_skin, color='deepskyblue', linewidth=2, label="True boundary")
                label_used.add("True boundary")

        # gp_part_name = f"IgaModelPart.SBM_Support_{side}"
        # if model.HasModelPart(gp_part_name):
        #     mp = model.GetModelPart(gp_part_name)
        #     if len(mp.Conditions) > 0:
        #         for cond in mp.Conditions:
        #             geom = cond.GetGeometry()
        #             N = geom.ShapeFunctionsValues()
        #             gp = np.zeros(2)
        #             for i, node in enumerate(geom):
        #                 gp[0] += N[0, i] * node.X
        #                 gp[1] += N[0, i] * node.Y
                    # plot gp, arrows, etc...

    try:
        is_plot_all = int(input("Plot all the active knot spans and gauss points? 1: yes, 0: no "))
    except ValueError:
        is_plot_all = 0  # defaultis_plot_all = input("plot all the active knot spans and gauss points ? 1: yes. 0: no")
        
    area = 0

    if is_plot_all == 1:
        # === ACTIVE GAUSS POINTS: FROM ELEMENT CENTERS ===
        mp = model["IgaModelPart"]
        active_gauss_points = []
        for cond in mp.Conditions:
            center = cond.GetGeometry().Center()
            x = center.X
            y = center.Y
            ax.plot(x, y, "go", markersize=5, label="ESBM Gauss Points" if "ESBM Gauss Points" not in label_used else "")
            label_used.add("ESBM Gauss Points")


        mp = model["IgaModelPart.StructuralAnalysisDomain"]
        active_gauss_points = []

        for elem in mp.Elements:
            center = elem.GetGeometry().Center()
            x = center.X
            y = center.Y
            active_gauss_points.append((x, y))
            ax.plot(x, y, "ro", markersize=1, label="Elemental Gauss Points" if "Elemental Gauss Points" not in label_used else "")
            label_used.add("Elemental Gauss Points")

            weight = elem.GetValue(KratosMultiphysics.INTEGRATION_WEIGHT)
            area+= weight

    
    mp = model["IgaModelPart.GapElements"]
    for elem in mp.Elements:
        center = elem.GetGeometry().Center()
        x = center.X
        y = center.Y
        ax.plot(x, y, "ro", markersize=2, label="Elemental Gauss Points" if "Elemental Gauss Points" not in label_used else "")
        weight = elem.GetValue(KratosMultiphysics.INTEGRATION_WEIGHT)
        # print(weight)
        area += weight
        
    print("Total area of the model:", area)


    mp = model["IgaModelPart.Layer0"]
    mp_skin = model["skin_model_part.inner"]
    perimeter = 0.0
    for cond in mp.Conditions:
        center = cond.GetGeometry().Center()
        x = center.X
        y = center.Y
        ax.plot(x, y, "bo", markersize=3, label="ESBM Gauss Points" if "ESBM Gauss Points" not in label_used else "")
        weight = cond.GetValue(KratosMultiphysics.INTEGRATION_WEIGHT)
        normal = cond.GetValue(KratosMultiphysics.NORMAL)
        # ax.arrow(x, y, normal[0] * 0.02, normal[1] * 0.02, head_width=0.001, head_length=0.02, fc='blue', ec='blue')
        
        projection = cond.GetValue(KratosMultiphysics.IgaApplication.PROJECTION_NODE_COORDINATES)
        
        if (len(projection) == 0): continue
        # ax.plot([x, projection[0]], [y, projection[1]], "g-o", markersize=3)
        
        interpolation_node_ids = cond.GetValue(KratosMultiphysics.IgaApplication.TEMP_INTERPOLATION_NODES_ID)
        for id in interpolation_node_ids:
            id =int(id)
            node = mp_skin.GetNode(id)
            
            ax.plot(node.X, node.Y, "ko", markersize=3)
        perimeter += weight

    mp = model["IgaModelPart.GapInterfaces"]
    # perimeter = 0.0
    for cond in mp.Conditions:
        center = cond.GetGeometry().Center()
        x = center.X
        y = center.Y
        ax.plot(x, y, "ko", markersize=2, label="ESBM Gauss Points" if "ESBM Gauss Points" not in label_used else "")
        weight = cond.GetValue(KratosMultiphysics.INTEGRATION_WEIGHT)
        normal = cond.GetValue(KratosMultiphysics.NORMAL)
        # ax.arrow(x, y, normal[0] * 0.1, normal[1] * 0.1, head_width=0.05, head_length=0.1, fc='blue', ec='blue')
        # perimeter += weight
        

    true_perimeter = 2*np.pi*0.17
    true_area = 1-np.pi*0.17*0.17
    print("Total perimeter of the model:", perimeter)
    # true_area = 16-np.pi #34
    # true_perimeter = 2*np.pi #4*np.sqrt(2)
    error_area = abs(area - true_area) / true_area
    error_perimeter = abs(perimeter - true_perimeter) / true_perimeter
    print("Area error:", error_area)
    print("Perimeter error:", error_perimeter)

    # FIXME:
    # with open("centers.txt", "r") as file:
    #     lines = file.readlines()
    #     for line in lines[:]:
    #         values = line.strip().split()

    #         x_coord_true_surrogate = [float(values[0]), float(values[3])]
    #         y_coord_true_surrogate = [float(values[1]), float(values[4])]
    #         ax.plot(x_coord_true_surrogate, y_coord_true_surrogate, "g--", markersize=2, label="SBM projections" if "SBM projections" not in label_used else "")
    # #         label_used.add("SBM projections")

    # === FIND ACTIVE CELLS AND COLOR THEM LIGHT YELLOW ===
    # active_cells = set()
    # for (x_gp, y_gp) in active_gauss_points:
    #     for i in range(len(knots_u_1) - 1):
    #         for j in range(len(knots_v_1) - 1):
    #             u_min, u_max = knots_u_1[i], knots_u_1[i+1]
    #             v_min, v_max = knots_v_1[j], knots_v_1[j+1]
    #             if u_min <= x_gp <= u_max and v_min <= y_gp <= v_max:
    #                 active_cells.add((i, j))

    # for (i, j) in active_cells:
    #     u_min, u_max = knots_u_1[i], knots_u_1[i+1]
    #     v_min, v_max = knots_v_1[j], knots_v_1[j+1]
    #     rect = patches.Rectangle((u_min, v_min), u_max - u_min, v_max - v_min,
    #                             linewidth=0, edgecolor=None, facecolor="#FFFACD", zorder=0,
    #                             label="Active elements" if "Active elements" not in label_used else "")
    #     ax.add_patch(rect)
    #     label_used.add("Active elements")

    # # Custom handles
    # custom_handles = []

    # if "Parameter space" in label_used:
    #     custom_handles.append(Line2D([0], [0], color='gray', linestyle=':', label='Parameter space'))

    # if "Surrogate boundary" in label_used:
    #     custom_handles.append(Line2D([0], [0], color='darkred', linewidth=2, label='Surrogate boundary'))

    # if "True boundary" in label_used:
    #     custom_handles.append(Line2D([0], [0], color='deepskyblue', linewidth=2, label='True boundary'))

    # if "Elemental Gauss Points" in label_used:
    #     custom_handles.append(Line2D([0], [0], marker='o', color='red', linestyle='None', markersize=5, label='Elemental Gauss Points'))

    # if "Boundary Gauss Points" in label_used:
    #     custom_handles.append(Line2D([0], [0], marker='o', color='red', linestyle='None', markersize=5, label='Boundary Gauss Points'))

    # if "Active elements" in label_used:
    #     custom_handles.append(Patch(facecolor='#FFFACD', edgecolor='goldenrod', linestyle=':', linewidth=1.5, label='Active elements'))
    
    # if "SBM projections" in label_used:
    #     custom_handles.append(Line2D([0], [0], color='green', linewidth=2, label='SBM projections'))

    # Clean layout for publication
    ax.set_aspect('equal')
    ax.set_xlabel("")  # or add label if meaningful
    ax.set_ylabel("")
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)  # hide ticks
    plt.tight_layout()
    plt.subplots_adjust(right=0.8)  # extra space for legend

    # ax.legend(handles=custom_handles, loc='center left', bbox_to_anchor=(1, 0.5))

    fig.tight_layout()
    fig.savefig("Gap_Sbm.pdf", format="pdf", dpi=1000, bbox_inches="tight")
    plt.show()


    
