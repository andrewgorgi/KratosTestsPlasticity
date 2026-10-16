import json, sys
from pathlib import Path
import numpy as np
import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication
import KratosMultiphysics.ConstitutiveLawsApplication as CLA

root = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(root))
from analytic_fields import evaluator
model=KM.Model();part=model.CreateModelPart('PointCheck');part.ProcessInfo[KM.DOMAIN_SIZE]=2
nodes=[part.CreateNewNode(i+1,x,y,0) for i,(x,y) in enumerate(((0.,0.),(1.,0.),(0.,1.)))]
geometry=KM.Triangle2D3(*nodes)
prop=part.CreateNewProperties(2)
material=json.loads((root/'materials_plastic.json').read_text())['properties'][0]['Material']
for name,value in material['Variables'].items():prop[KM.KratosGlobals.GetVariable(name)]=value
law=CLA.SmallStrainIsotropicPlasticityPlaneStrainVonMisesVonMises()
law.InitializeMaterial(prop,geometry,KM.Vector([1/3]*3))
parameters=KM.ConstitutiveLawParameters(geometry,prop,part.ProcessInfo)
flags=KM.Flags()
for flag in (KM.ConstitutiveLaw.USE_ELEMENT_PROVIDED_STRAIN,KM.ConstitutiveLaw.COMPUTE_STRESS,KM.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR):flags.Set(flag,True)
parameters.SetOptions(flags)
strain=KM.Vector(3);stress=KM.Vector(3);matrix=KM.Matrix(3,3)
parameters.SetStrainVector(strain);parameters.SetStressVector(stress);parameters.SetConstitutiveMatrix(matrix)
rows=[]
for time in (.00125,.0025):
 part.ProcessInfo[KM.TIME]=time
 _,grad,expected=evaluator('cosine_boundary','plastic')(np.array([[2.,0.]]),time)
 strain[0]=grad[0,0];strain[1]=grad[0,3];strain[2]=grad[0,1]+grad[0,2]
 law.InitializeMaterialResponseCauchy(parameters)
 law.CalculateMaterialResponseCauchy(parameters)
 reference=expected[0,:3]
 row=dict(time=time,strain=list(strain),stress_incremental=list(stress),stress_analytic=reference.tolist(),relative_difference=float(np.linalg.norm(np.array(stress)-reference)/np.linalg.norm(reference)))
 rows.append(row)
 law.FinalizeMaterialResponseCauchy(parameters)
print('MATERIAL_CHECK',json.dumps(rows),flush=True)
(Path(__file__).resolve().parent / 'material_point_check.json').write_text(json.dumps(rows,indent=2)+'\n')
