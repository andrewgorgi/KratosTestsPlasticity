"""Local finite-difference Jacobian checks at a plastic affine trial state.

No material history is committed during the central-difference probes.
Uses the current mixed JSON/material without modifying them.
"""
import sys,json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import annulus_study as s
s.KM.ParallelUtilities.SetNumThreads(1)
rows=[]
for degree in (1,2):
 args=SimpleNamespace(parameters=None,boundary='mixed',degree=degree,dt=None)
 folder=s.ROOT/'sbm';p=s.load_parameters(folder,args,[10,10])
 with s.shared.working_directory(folder):
  model=s.KM.Model();a=s.Analysis(model,p);a.Initialize()
  a.time=.001;model['IgaModelPart'].CloneTimeStep(.001);a.InitializeSolutionStep()
  part=model['IgaModelPart']
  for node in part.Nodes:
   node.SetSolutionStepValue(s.KM.DISPLACEMENT_X,.01*node.X0+.003*node.Y0)
   node.SetSolutionStepValue(s.KM.DISPLACEMENT_Y,.002*node.X0-.004*node.Y0)
  used={}
  for cond in part.Conditions:
   name=s.condition_name(cond).strip('"')
   if name not in ('SbmLoadSolidCondition','SbmSolidCondition') or used.get(name,0)>=3:continue
   used[name]=used.get(name,0)+1
   dofs=cond.GetDofList(part.ProcessInfo)
   variables=[dof.GetVariable() for dof in dofs]
   nodes=[part.GetNode(dof.Id()) for dof in dofs]
   base=np.array([node.GetSolutionStepValue(var) for node,var in zip(nodes,variables)])
   rng=np.random.default_rng(45);direction=rng.normal(size=len(base));direction/=np.linalg.norm(direction)
   K=s.KM.Matrix();R=s.KM.Vector();cond.CalculateLocalSystem(K,R,part.ProcessInfo)
   expected=np.array(K)@direction
   for h in (1e-7,1e-8):
    values=[]
    for sign in (1,-1):
     for node,var,value in zip(nodes,variables,base+sign*h*direction):node.SetSolutionStepValue(var,float(value))
     rhs=s.KM.Vector();cond.CalculateRightHandSide(rhs,part.ProcessInfo);values.append(np.array(rhs))
    for node,var,value in zip(nodes,variables,base):node.SetSolutionStepValue(var,float(value))
    fd=-(values[0]-values[1])/(2*h)
    error=np.linalg.norm(fd-expected)/max(np.linalg.norm(fd),np.linalg.norm(expected),1e-30)
    row=dict(degree=degree,condition=name,id=cond.Id,h=h,relative_tangent_error=float(error),fd_norm=float(np.linalg.norm(fd)))
    rows.append(row);print('TANGENT_AUDIT',json.dumps(row),flush=True)
output=s.ROOT/'diagnostics/local_tangent.json'
output.parent.mkdir(parents=True,exist_ok=True)
output.write_text(json.dumps(rows,indent=2)+'\n')
