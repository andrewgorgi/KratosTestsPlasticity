"""Generate only the annular benchmark; leave the rectangular cases untouched."""
import json
import sys
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent/'AnalyticalTestRevised'))
spec = importlib.util.spec_from_file_location('rectangle_builder', ROOT.parent/'AnalyticalTestRevised/build_parameters.py')
rectangle = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rectangle)
# Import local constants explicitly after loading the rectangular helper.
from analytic import E, NU, YIELD, H, END_TIME


def build(method, boundary='dirichlet'):
    params = rectangle.build(method, 'cosine_boundary', 'plastic_hardening')
    params['problem_data']['problem_name'] = f'{method}_annulus_{boundary}_plastic_hardening'
    params['problem_data']['end_time'] = END_TIME
    solver = params['solver_settings']
    solver['material_import_settings']['materials_filename'] = '../materials.json'
    solver['time_stepping']['time_step'] = END_TIME / 20
    if method == 'sbm' or (method == 'gap_sbm' and boundary == 'dirichlet'):
        solver['solving_strategy_settings']['type'] = 'line_search'
        solver['time_stepping']['time_step'] = END_TIME / 2
    params['processes']['additional_processes'] = []
    geometry = params['modelers'][1]['Parameters']
    geometry.update(lower_point_xyz=[-2.2,-2.2,0.], upper_point_xyz=[2.2,2.2,0.],
                    number_of_knot_spans=[20,20], number_of_inner_loops=1,
                    number_initial_points_if_importing_nurbs=1000,
                    skin_model_part_inner_initial_name='initial_skin_inner',
                    skin_model_part_outer_initial_name='initial_skin_outer')
    if method == 'sbm':
        # Keep the Dirichlet surrogate closer to the curved boundary.
        # lambda=1 failed on the p=1, 80x80 annulus even with line search.
        geometry['lambda_inner'] = .5 if boundary == 'dirichlet' else .5
        geometry['lambda_outer'] = .5 if boundary == 'dirichlet' else .5
    imports = []
    for side in ('outer','inner'):
        condition = ('GapSbmLoadSolidCondition' if method == 'gap_sbm' and side == 'inner' and boundary == 'mixed' else
                     'GapSbmSolidCondition' if method == 'gap_sbm' else
                     'SbmLoadSolidCondition' if side == 'inner' and boundary == 'mixed' else
                     'SbmSolidCondition')
        imports.append({'modeler_name':'ImportNurbsSbmModeler','Parameters':{
            'input_filename':f'../../AnalyticalTestRevised/geometry/{side}_circle.json',
            'model_part_name':f'initial_skin_{side}',
            'link_layer_to_condition_name':[{'layer_name':'Layer0','condition_name':condition}]}})
    iga = params['modelers'][2]
    entries = iga['Parameters']['element_condition_list']
    entries[0]['number_of_integration_points_per_span'] = 4
    if method == 'sbm':
        inner = json.loads(json.dumps(entries[1]))
        inner['iga_model_part'] = 'SBM_Inner'
        inner['sbm_parameters']['is_inner'] = True
        entries.append(inner)
    params['modelers'] = imports + [params['modelers'][1], iga]
    return params


def main():
    material = {'properties':[{'model_part_name':'IgaModelPart','properties_id':2,'Material':{
        'constitutive_law':{'name':'SmallStrainJ2PlasticityPlaneStrain2DThreeComponentsLaw'},
        'Variables':{'DENSITY':1.,'YOUNG_MODULUS':E,'POISSON_RATIO':NU,'THICKNESS':1.,
                     'PENALTY_FACTOR':100*E,'YIELD_STRESS':YIELD,'ISOTROPIC_HARDENING_MODULUS':H,
                     'EXPONENTIAL_SATURATION_YIELD_STRESS':YIELD,'HARDENING_EXPONENT':0.},
        'Tables':{}}}]}
    (ROOT/'materials.json').write_text(json.dumps(material,indent=2)+'\n')
    for method, boundary in [('sbm','dirichlet'),('sbm','mixed'),('gap_sbm','dirichlet'),('gap_sbm','mixed')]:
        path = ROOT/method/f'ProjectParameters_{boundary}.json'
        path.write_text(json.dumps(build(method,boundary),indent=2)+'\n')


if __name__ == '__main__':
    main()
