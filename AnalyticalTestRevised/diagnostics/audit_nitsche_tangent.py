"""Audit della tangente locale Nitsche. Eseguire con PYTHONPATH della build KratosPR."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import study
import numpy as np
K=study.KM
p=K.Parameters((study.ROOT/'gap_sbm/ProjectParameters_cosine_boundary_plastic.json').read_text())
for i in range(2):study.geometry_settings(p)['number_of_knot_spans'][i].SetInt(10)
K.ParallelUtilities.SetNumThreads(1)
with study.working_directory(study.ROOT/'gap_sbm'):
 model=K.Model();a=study.Analysis(model,p);a.Initialize()
 mp=model['IgaModelPart']; info=mp.ProcessInfo
 print('PARTS',[(part.FullName(),part.NumberOfConditions()) for part in mp.SubModelParts],flush=True)
 conditions=sorted([c for c in mp.Conditions if 'CutSbmSolidCondition' in str(c) and abs(c.GetValue(K.NORMAL)[1])>0.9],key=lambda c:100*abs(c.GetGeometry().Center().X-2)+abs(c.GetGeometry().Center().Y))[:3]
 conditions += sorted([c for c in mp.Conditions if 'CutSbmSolidInterfaceCondition' in str(c)],key=lambda c:100*abs(c.GetGeometry().Center().X-2)+abs(c.GetGeometry().Center().Y))[:3]
 print('SELECTED',len(conditions),flush=True)
 print('BOUNDARY_COUNTS',[(str(c).splitlines()[0],c.GetGeometry().Center().X,c.GetGeometry().Center().Y) for c in conditions],flush=True)
 for time in [0.0001,0.00125]:
  a.time=time;info[K.TIME]=time;info[K.STEP]=2; info[K.NL_ITERATION_NUMBER]=2
  a.InitializeSolutionStep()
  nodes=list(mp.Nodes);vals=study.evaluator('cosine_boundary','plastic')([[n.X0,n.Y0] for n in nodes],time)[0]
  for n,u in zip(nodes,vals):
   n.SetSolutionStepValue(K.DISPLACEMENT_X,float(u[0]));n.SetSolutionStepValue(K.DISPLACEMENT_Y,float(u[1]))
  for c in conditions:
   dofs=c.GetDofList(info);lhs=K.Matrix();rhs=K.Vector();c.CalculateLocalSystem(lhs,rhs,info);base=np.array(lhs)
   keys=[(d.Id,str(d.GetVariable())) for d in dofs]
   expected=np.column_stack([base[:,[k for k,key in enumerate(keys) if key==keys[j]]].sum(axis=1) for j in range(len(dofs))])
   for h in [1e-7,1e-8]:
    numerical=np.zeros_like(base)
    for j,d in enumerate(dofs):
     old=d.GetSolutionStepValue();d.SetSolutionStepValue(old+h);rp=K.Vector();c.CalculateRightHandSide(rp,info)
     d.SetSolutionStepValue(old-h);rm=K.Vector();c.CalculateRightHandSide(rm,info)
     d.SetSolutionStepValue(old);numerical[:,j]=-(np.array(rp)-np.array(rm))/(2*h)
    print('AUDIT',time,c.Id,'x',c.GetGeometry().Center().X,'h',h,'relative_error',np.linalg.norm(expected-numerical)/np.linalg.norm(numerical),flush=True)
