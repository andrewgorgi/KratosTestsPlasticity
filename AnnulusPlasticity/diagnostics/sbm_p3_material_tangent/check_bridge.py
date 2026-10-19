"""Verify wrapper stress, history and Clone against the native law over a cycle."""
import json,sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]))
import check_projection_history as check
KM=check.KM
import KratosMultiphysics.ConstitutiveLawsApplication as CLA
sys.path.insert(0,json.loads((HERE/'bridge.json').read_text())['directory'])
import sbm_tangent_bridge as bridge
model=KM.Model();part=model.CreateModelPart('material_test')
part.CreateNewNode(1,0.,0.,0.);part.CreateNewNode(2,1.,0.,0.)
props=part.CreateNewProperties(1)
params=KM.Parameters((HERE/'local/materials_input.json').read_text())
variables=params['properties'][0]['Material']['Variables']
for name in ('YOUNG_MODULUS','POISSON_RATIO','YIELD_STRESS','ISOTROPIC_HARDENING_MODULUS','EXPONENTIAL_SATURATION_YIELD_STRESS','HARDENING_EXPONENT'):
    props.SetValue(KM.KratosGlobals.GetVariable(name),variables[name].GetDouble())
props.SetValue(KM.CONSTITUTIVE_LAW,KM.KratosGlobals.GetConstitutiveLaw('SmallStrainJ2PlasticityPlaneStrain2DThreeComponentsLaw').Clone())
e=part.CreateNewCondition('LineCondition2D2N',1,[1,2],props)
points=[check.MaterialPoint(e,part.ProcessInfo) for _ in range(3)]
shape=KM.Vector(np.array(e.GetGeometry().ShapeFunctionsValues())[0])
for point,order in zip(points[1:],(1,2)):
    point.law=bridge.NumericalJ2(order,-1.);point.law.InitializeMaterial(props,e.GetGeometry(),shape)
max_stress=0.;max_alpha=0.;max_clone_stress=0.
for eps in ([1e-5,0.,0.],[.005,-.002,.003],[.006,-.004,.007],[.004,-.001,.002],[-.006,.003,-.008],[0.,0.,0.],[.008,-.002,.01]):
    eps=np.array(eps);reference=points[0].response(eps)[0]
    for point in points:
        max_stress=max(max_stress,float(np.linalg.norm(point.response(eps)[0]-reference)))
        point.commit(eps)
    alphas=[point.law.GetValue(CLA.EQUIVALENT_PLASTIC_STRAIN) for point in points]
    max_alpha=max(max_alpha,max(abs(a-alphas[0]) for a in alphas))
    for point in points:
        clone=point.law.Clone();stress=point.response(eps)[0];clone.CalculateMaterialResponseCauchy(point.values)
        max_clone_stress=max(max_clone_stress,float(np.linalg.norm(np.array(point.stress)-stress)))
assert max_stress<1e-10 and max_alpha<1e-12 and max_clone_stress<1e-10
result=dict(max_stress_difference=max_stress,max_committed_alpha_difference=max_alpha,max_clone_stress_difference=max_clone_stress)
(HERE/'bridge_validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
