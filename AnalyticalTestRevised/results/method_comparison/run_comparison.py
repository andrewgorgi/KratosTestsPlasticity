import sys,json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
root=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(root))
import study
study.KM.ParallelUtilities.SetNumThreads(4)
mode=sys.argv[1];method='body_fitted' if mode=='body_fitted' else 'sbm';folder=root/method
p=study.load_parameters(folder,SimpleNamespace(parameters=None,field='cosine_middle',physics='plastic',spans=[10,10],degree=2))
p['modelers'][p['modelers'].size()-1]['Parameters']['element_condition_list'][0]['number_of_integration_points_per_span'].SetInt(10)
if mode=='sbm_aligned':
 g=study.geometry_settings(p)
 for key,values in [('lower_point_xyz',[0.,0.,0.]),('lower_point_uvw',[0.,0.,0.]),('upper_point_xyz',[4.,2.,0.]),('upper_point_uvw',[4.,2.,0.])]:g[key].SetVector(study.KM.Vector(values))
(Path(__file__).resolve().parent/(mode+'_parameters.json')).write_text(p.PrettyPrintJsonString()+'\n')
iterations=[]
class Checked(study.Analysis):
 def Initialize(self):
  super().Initialize();solver=self._GetSolver();old=solver.SolveSolutionStep
  def checked():
   ok=old();iterations.append(dict(time=self.time,iterations=solver.GetComputingModelPart().ProcessInfo[study.KM.NL_ITERATION_NUMBER],converged=bool(ok)))
   if not ok:raise RuntimeError('Newton not converged')
   return ok
  solver.SolveSolutionStep=checked
study.Analysis=Checked
s,m=study.solve(folder,p)
r=dict(mode=mode,**study.errors(s,'cosine_middle','plastic'),iterations=iterations)
(Path(__file__).resolve().parent/(mode+'_solve.json')).write_text(json.dumps(r,indent=2));print('SUMMARY',json.dumps(r),flush=True)
