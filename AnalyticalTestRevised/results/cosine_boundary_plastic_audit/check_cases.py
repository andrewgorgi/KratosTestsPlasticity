import json,sys,tempfile
from pathlib import Path
from types import SimpleNamespace
from collections import Counter
import numpy as np
import sympy as sp
root=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(root))
import study
from analytic_fields import expressions,x,y,t
study.KM.ParallelUtilities.SetNumThreads(4)
configs={method:json.loads(study.KM.Parameters((root/method/'ProjectParameters_cosine_boundary_plastic.json').read_text()).WriteJsonString()) for method in ('body_fitted','sbm','gap_sbm')}
reference=configs['body_fitted']
assignments={}
summary={}
for method,config in configs.items():
 assert config['solver_settings']==reference['solver_settings']
 assert config['problem_data']['start_time']==reference['problem_data']['start_time']
 assert config['problem_data']['end_time']==reference['problem_data']['end_time']
 entries=config['processes']['additional_processes'][0]['Parameters']['element_condition_list']
 assignments[method]={}
 for item in entries:
  for variable in item['variables']:
   assignments[method].setdefault(variable['variable_name'],[]).append(variable['value'])
 for variable,values in assignments[method].items():
  assert all(v==assignments['body_fitted'][variable][0] for v in values)
 geometry=next(m['Parameters'] for m in config['modelers'] if 'number_of_knot_spans' in m['Parameters'])
 bulk=next(m['Parameters']['element_condition_list'][0] for m in config['modelers'] if m['modeler_name']=='IgaModelerSbm')
 summary[method]=dict(geometry=geometry,bulk_quadrature_per_direction=bulk.get('number_of_integration_points_per_span','degree+1'),solver_settings_equal=True,formula_strings_equal=True)

body=sp.lambdify((x,y,t),list(expressions('cosine_boundary','plastic')[3]),'numpy')
u=study.evaluator('cosine_boundary','plastic')
formula_checks=[]
points=[(.4,.2),(2.,0.),(2.,.7),(2.,1.),(3.4,2.)]
for time in (.00001,.00125,.0025):
 for point in points:
  for variable,formula in assignments['body_fitted'].items():
   actual=np.array([study.KM.GenericFunctionUtility(s).CallFunction(*point,0.,time,*point,0.) for s in formula[0]])
   with np.errstate(divide='ignore',invalid='ignore'):
    expected=np.array(body(*point,time),dtype=float) if variable=='BODY_FORCE' else np.r_[u(np.array([point]),time)[0][0],0.]
   assert np.allclose(actual,expected,rtol=1e-10,atol=1e-10),(variable,time,point,actual,expected)
   formula_checks.append(dict(variable=variable,time=time,point=point,maximum_absolute_difference=float(max(abs(actual-expected)))))

with tempfile.TemporaryDirectory(prefix='cosine_case_audit_',dir='/tmp') as temporary:
 workspace=Path(temporary)
 (workspace/'materials_plastic.json').symlink_to(root/'materials_plastic.json')
 (workspace/'geometry').symlink_to(root/'geometry',target_is_directory=True)
 for method in configs:
  folder=workspace/method;folder.mkdir();(folder/'txt_files').mkdir()
  parameters=study.load_parameters(root/method,SimpleNamespace(parameters=None,field='cosine_boundary',physics='plastic',spans=[8,4],degree=None))
  with study.working_directory(folder):
   model=study.KM.Model();analysis=study.Analysis(model,parameters);analysis.Initialize()
   solver=analysis._GetSolver();analysis.time=solver.AdvanceInTime(analysis.time);analysis.InitializeSolutionStep()
   counts=Counter();traction_checks=[]
   for condition in model['IgaModelPart'].Conditions:
    name=str(condition).split(' #')[0].strip('"');counts[name]+=1
    if 'LoadSolidCondition' not in name:continue
    center=condition.GetGeometry().Center()
    point=np.array([center.X,center.Y])
    if name=='SbmLoadSolidCondition':point=np.array(condition.GetValue(study.IGA.PROJECTION_NODE_COORDINATES))[:2]
    assert np.isclose(point[1],0.) or np.isclose(point[1],2.),(method,name,point)
    # Escludere le proiezioni esattamente al vertice, dove la normale è ambigua.
    if np.isclose(point[0],0.) or np.isclose(point[0],4.):continue
    normal=np.array(condition.GetValue(study.KM.NORMAL))[:2]
    true_normal=np.array([0.,-1. if point[1]<1. else 1.])
    assert float(normal@true_normal)>.99,(method,name,point,normal)
    rhs=study.KM.Vector();condition.CalculateRightHandSide(rhs,model['IgaModelPart'].ProcessInfo)
    actual=np.array(rhs).reshape(-1,2).sum(axis=0)/condition.GetValue(study.KM.INTEGRATION_WEIGHT)
    stress=u(point[None,:],analysis.time)[2][0]
    expected=np.array([stress[0]*true_normal[0]+stress[2]*true_normal[1],stress[2]*true_normal[0]+stress[1]*true_normal[1]])
    assert np.allclose(actual,expected,rtol=1e-10,atol=1e-9),(method,name,point,actual,expected)
    traction_checks.append(dict(point=point.tolist(),traction=actual.tolist(),maximum_absolute_difference=float(max(abs(actual-expected)))))
   summary[method].update(runtime_mesh=[8,4],condition_counts=dict(counts),neumann_checks=len(traction_checks),maximum_traction_absolute_difference=max(v['maximum_absolute_difference'] for v in traction_checks),neumann_checks_details=traction_checks)
   print('CASE_AUDIT',method,json.dumps({k:v for k,v in summary[method].items() if k not in ('geometry','neumann_checks_details')}),flush=True)

report=dict(methods=summary,formula_checks=formula_checks)
(Path(__file__).resolve().parent / 'case_audit.json').write_text(json.dumps(report,indent=2)+'\n')
