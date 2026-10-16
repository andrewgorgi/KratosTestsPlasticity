"""Evaluate the assembled RHS at exactly representable affine solutions.

The residual is normalized by the norm of the load-condition RHS vector.
No Newton solve or material finalization is performed.
"""
import sys,json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import annulus_study as s
from diagnose_mixed import reference,analysis_type
s.KM.ParallelUtilities.SetNumThreads(1)
from tempfile import TemporaryDirectory
temporary=TemporaryDirectory()
rows=[]
for case in ('patch_elastic','patch_hardening','patch_perfect'):
 for degree in (1,2):
  for n in (10,20):
   folder=s.ROOT/'sbm';args=SimpleNamespace(parameters=None,boundary='mixed',degree=degree,dt=None)
   p=s.load_parameters(folder,args,[n,n])
   material_data=json.loads(s.KM.Parameters((s.ROOT/'materials.json').read_text()).PrettyPrintJsonString())
   if case=='patch_perfect':
    material_data['properties'][0]['Material']['Variables']['ISOTROPIC_HARDENING_MODULUS']=0.
   material=Path(temporary.name)/f'{case}_{degree}_{n}.json'
   material.write_text(json.dumps(material_data))
   p['solver_settings']['material_import_settings']['materials_filename'].SetString(str(material))
   evaluate,body=reference(case)
   with s.shared.working_directory(folder):
    m=s.KM.Model();a=analysis_type(evaluate,body)(m,p);a.Initialize()
    part=m['IgaModelPart'];a.time=.002;part.CloneTimeStep(.002)
    nodes=list(part.Nodes);u=evaluate([[node.X0,node.Y0] for node in nodes],.002)[0]
    for node,value in zip(nodes,u):
     node.SetSolutionStepValue(s.KM.DISPLACEMENT_X,float(value[0]));node.SetSolutionStepValue(s.KM.DISPLACEMENT_Y,float(value[1]))
    a.InitializeSolutionStep()
    total={};loads={};sum_norm=0.
    entities=[e for e in part.Elements if e.GetValue(s.KM.ACTIVATION_LEVEL)!=2]
    entities += [c for c in part.Conditions if s.condition_name(c).strip('"') in ('SbmLoadSolidCondition','SbmSolidCondition')]
    for entity in entities:
     rhs=s.KM.Vector();entity.CalculateRightHandSide(rhs,part.ProcessInfo)
     sum_norm+=np.linalg.norm(np.array(rhs))
     is_load='SbmLoadSolidCondition' in s.condition_name(entity) if isinstance(entity,s.KM.Condition) else False
     for dof,value in zip(entity.GetDofList(part.ProcessInfo),rhs):
      key=(dof.Id(),dof.GetVariable().Name());total[key]=total.get(key,0.)+value
      if is_load:loads[key]=loads.get(key,0.)+value
    error=np.linalg.norm(list(total.values()))/max(np.linalg.norm(list(loads.values())),1e-30)
    row=dict(case=case,degree=degree,spans=n,relative_exact_patch_residual=float(error),residual_norm=float(np.linalg.norm(list(total.values()))))
    rows.append(row);print('PATCH_RESIDUAL',json.dumps(row),flush=True)
output=s.ROOT/'diagnostics/exact_patch_residual.json'
output.parent.mkdir(parents=True,exist_ok=True)
output.write_text(json.dumps(rows,indent=2)+'\n')
temporary.cleanup()
