"""Rigenera le configurazioni commentate da un'unica definizione analitica."""
import json
from pathlib import Path
import sympy as sp
from analytic_fields import FIELDS, PHYSICS, HARDENING, expressions, expression_string

ROOT = Path(__file__).resolve().parent

def entry(part, variable, values):
    return {'iga_model_part': 'IgaModelPart.'+part,
            'variables': [{'variable_name': variable, 'value': values}]}

def build(method, field, physics):
    u, _, _, body = expressions(field, physics)
    displacement = [expression_string(v) for v in u]+['0.0']
    force = [expression_string(sp.simplify(v) if physics == 'plastic_hardening' else v)
             for v in body]
    lower, upper = ([0.,0.,0.], [4.,2.,0.]) if method == 'body_fitted' else ([-.5,-.3,0.], [4.2,2.1,0.])
    geometry = dict(echo_level=0, model_part_name='IgaModelPart', lower_point_xyz=lower,
                    upper_point_xyz=upper,
                    polynomial_order=[2,2], number_of_knot_spans=[20,10])
    modelers = []
    if method != 'body_fitted':
        cut = 'Gap' if method == 'gap_sbm' else ''
        modelers.append({'modeler_name':'ImportNurbsSbmModeler', 'Parameters':{
            'input_filename':'../geometry/rectangle.nurbs.json', 'model_part_name':'initial_skin_model_part_out',
            'link_layer_to_condition_name':[{'layer_name':side, 'condition_name':cut+'SbmSolidCondition'} for side in ('bottom','right','top','left')]}})
        geometry.update(number_of_inner_loops=0,
                        number_initial_points_if_importing_nurbs=100,
                        skin_model_part_outer_initial_name='initial_skin_model_part_out', skin_model_part_name='skin_model_part')
        if method == 'sbm':
            geometry['lambda_outer'] = 1.
    if method == 'gap_sbm':
        geometry.update(number_internal_divisions=1, gap_approximation_order=1,
                        gap_sbm_type='interpolation', gap_element_name='GapSbmSolidElement',
                        gap_interface_condition_name='GapSbmSolidInterfaceCondition')
    modelers.append({'modeler_name':'NurbsGeometryModelerGapSbm' if method == 'gap_sbm' else 'NurbsGeometryModelerSbm', 'Parameters':geometry})
    elements = [{'geometry_type':'GeometrySurface', 'iga_model_part':'StructuralAnalysisDomain',
                 'type':'element', 'name':'SolidElement', 'shape_function_derivatives_order':2}]
    if method != 'gap_sbm':
        support = {'geometry_type':'SurfaceEdge', 'iga_model_part':'SBM_Support_outer',
                   'type':'condition', 'name':'SupportSolidCondition' if method == 'body_fitted' else 'SbmCondition',
                   'shape_function_derivatives_order':7}
        if method == 'body_fitted':
            support['brep_ids'] = [2,3,4,5]
        else:
            support['sbm_parameters'] = {'is_inner':False, 'projection_type':'linealized'}
        elements.append(support)
    modelers.append({'modeler_name':'IgaModelerSbm', 'Parameters':{'echo_level':0,
        'skin_model_part_name':'skin_model_part', 'analysis_model_part_name':'IgaModelPart', 'element_condition_list':elements}})
    assignments = [entry('StructuralAnalysisDomain', 'BODY_FORCE', force)]
    if method == 'gap_sbm':
        assignments += [entry('GapElements', 'BODY_FORCE', force)]
        assignments += [entry(side, 'DISPLACEMENT', displacement) for side in ('bottom','right','top','left')]
    else:
        assignments += [entry('SBM_Support_outer', 'DISPLACEMENT', displacement)]
    solver = {'model_part_name':'IgaModelPart', 'domain_size':2, 'echo_level':0, 'buffer_size':2,
              'analysis_type':'non_linear', 'model_import_settings':{'input_type':'use_input_model_part'},
              'material_import_settings':{'materials_filename':'../materials_'+physics+'.json'},
              'time_stepping':{'time_step':.0002}, 'rotation_dofs':False, 'reform_dofs_at_each_step':False,
              'compute_reactions':False, 'clear_storage':False, 'move_mesh_flag':False,
              'solving_strategy_settings':{'type':'newton_raphson','advanced_settings':{}},
              'builder_and_solver_settings':{'type':'block','advanced_settings':{}},
              'convergence_criterion':'and_criterion', 'displacement_relative_tolerance':1e-10,
              'displacement_absolute_tolerance':1e-12, 'residual_relative_tolerance':1e-9,
              'residual_absolute_tolerance':1e-10, 'max_iteration':100, 'solver_type':'Static',
              'linear_solver_settings':{'solver_type':'LinearSolversApplication.sparse_lu'},
              'auxiliary_variables_list':[], 'auxiliary_dofs_list':[], 'auxiliary_reaction_list':[]}
    if method == 'sbm' and physics == 'plastic_hardening':
        # Il passo 2e-4 non converge sulla mesh 10x10 oltre lo snervamento.
        solver['time_stepping']['time_step'] = .0001
    return {'problem_data':{'problem_name':method+'_'+field+'_'+physics,'echo_level':0,'parallel_type':'OpenMP',
                            'start_time':0.,'end_time':.002},
            'solver_settings':solver, 'modelers':modelers,
            'processes':{'additional_processes':[{'python_module':'assign_iga_external_conditions_process',
                'kratos_module':'KratosMultiphysics.IgaApplication','process_name':'AssignIgaExternalConditionsProcess',
                'Parameters':{'echo_level':0,'element_condition_list':assignments}}],
                'dirichlet_process_list':[], 'constraints_process_list':[]},
            'output_processes':{'output_process_list':[]}}

DESCRIPTIONS = {
 'linear':'Campo lineare nello spazio: u=(x*t,0). Patch test; forza di volume nulla.',
 'polynomial':'Campo polinomiale storico: u=(x*(x-4)*y*(y-2)*t,0). Non e un campo lineare.',
 'cosine_middle':'u=(-2*cos(pi*x/4)*sin(pi*y/2)*t,4*sin(pi*x/4)*cos(pi*y/2)*t). Massimo stress trial a y=1: snervamento nel mezzo.',
 'cosine_boundary':'Stesso campo con y sostituito da y+1. Massimo stress trial a y=0,2: snervamento al boundary.'}

def write(method, field, physics):
    params = build(method,field,physics)
    header = ('// '+method+'; campo='+field+'; materiale='+physics+'\n// '+DESCRIPTIONS[field]+'\n'
              '// Materiale distinto dal campo: "linear" indica u(x), non la legge costitutiva.\n'
              '// Dominio fisico [0,4] x [0,2]. Dirichlet esatto su TUTTO il contorno.\n'
              '// BODY_FORCE=-div(sigma), densita=1; stessi campi per bulk e GapElements.\n'
              '// plastic: J2 perfettamente plastico, E=67000, nu=0, sigma_y=201; t>=0 monotono.\n'
              '// elastic: plane strain lineare, E=67000, nu=0; nessuno snervamento.\n'
              + (f'// plastic_hardening: J2 isotropo lineare, H={HARDENING}, sigma_y=201; t>=0 monotono.\n'
                 if physics == 'plastic_hardening' else '') +
              '// Per cambiare campo/materiale scegliere un altro file completo; non decommentare formule isolate.\n')
    text = header+json.dumps(params,indent=2)+'\n'
    (ROOT/method/f'ProjectParameters_{field}_{physics}.json').write_text(text)
    if field == 'linear' and physics == 'elastic':
        (ROOT/method/'ProjectParameters.json').write_text(text)

def write_material(physics):
    variables = dict(DENSITY=1.,YOUNG_MODULUS=67000.,POISSON_RATIO=0.,THICKNESS=1.,PENALTY_FACTOR=-1.)
    law = 'LinearElasticPlaneStrain2DLaw'
    if physics == 'plastic':
        law = 'SmallStrainIsotropicPlasticityPlaneStrainVonMisesVonMises'
        variables.update(YIELD_STRESS=201., FRACTURE_ENERGY=1e10, HARDENING_CURVE=3, TANGENT_OPERATOR_ESTIMATION=2)
    elif physics == 'plastic_hardening':
        law = 'SmallStrainJ2PlasticityPlaneStrain2DThreeComponentsLaw'
        variables.update(YIELD_STRESS=201., ISOTROPIC_HARDENING_MODULUS=float(HARDENING),
                         EXPONENTIAL_SATURATION_YIELD_STRESS=201., HARDENING_EXPONENT=0.)
    # Gap e condizioni ricevono lo stesso property ID del bulk, incluso il materiale.
    materials = {'properties':[{'model_part_name':'IgaModelPart','properties_id':2,
                 'Material':{'constitutive_law':{'name':law}, 'Variables':variables, 'Tables':{}}}]}
    (ROOT/f'materials_{physics}.json').write_text(json.dumps(materials,indent=2)+'\n')

if __name__ == '__main__':
    source = json.loads((ROOT.parent/'nurbs_files/simplestBox.json').read_text())
    for line, side in zip(source['Lines'], ('bottom','right','top','left')):
        line['Layer']=side
    (ROOT/'geometry/rectangle.nurbs.json').write_text(json.dumps(source,indent=2)+'\n')
    for physics in PHYSICS:
        write_material(physics)
    for method in ('body_fitted','sbm','gap_sbm'):
        for field in FIELDS:
            for physics in PHYSICS:
                if physics == 'plastic_hardening' and field != 'cosine_boundary':
                    continue
                write(method,field,physics)
